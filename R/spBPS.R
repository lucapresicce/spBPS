utils::globalVariables("i")

#' Unified spatial BPS workflow (multivariate path, works for q = 1)
#'
#' @description
#' Orchestrates subsetting, local stacking weight estimation, global stacking
#' combination, and optional posterior or predictive simulation using the
#' multivariate Student-t spatial model. Works for both multivariate outcomes
#' and the univariate case via `q = 1`.
#'
#' @param data List with matrices `Y` (response) and `X` (covariates).
#' @param priors List of priors for the multivariate model (`mu_B`, `V_r`,
#'   `Psi`, `nu`).
#' @param coords Matrix of observation coordinates.
#' @param hyperpar List with elements `alpha` and `phi` (vectors allowed).
#' @param subset_size Target subset size when `K` is not provided. Default 500.
#' @param K Optional number of subsets. When `NULL`, computed as
#'   `ceiling(nrow(Y) / subset_size)` and lower-bounded at 1.
#' @param cv_folds Number of folds for local cross-validation (default 5).
#' @param rp Fraction of rows used when recomputing global stacking weights
#'   (passed to `BPS_combine`). Ignored when `combine_method = "pseudoBMA"`.
#' @param combine_method Choose between Bayesian Predictive Stacking (`"bps"`)
#'   or pseudo-BMA (`"pseudoBMA"`) for combining subsets.
#' @param draws Number of joint posterior/predictive draws to return (0 to
#'   skip). When positive, `newdata` must be supplied because draws are obtained
#'   via `BPS_post_MvT` which jointly samples posterior and predictive.
#' @param newdata Optional list with `X` and `coords` for prediction locations;
#'   required when either draw count is positive.
#' @param include_latent Logical; if `TRUE`, posterior draws include latent
#'   processes.
#' @param cores Optional integer; when >1 a parallel backend is registered
#'   internally via `doParallel::registerDoParallel(cores)` for the fit and draw
#'   loops. When `NULL`, the existing foreach backend (if any) is used.
#'
#' @return List with components `subsets`, `weights_global`, `weights_local`,
#'   `epd`, and optional `posterior` and `predictive` draws.
#'
#' @examples
#' \donttest{
#' n <- 1000
#' p <- 2
#' q <- 1
#'
#' Y <- matrix(rnorm(n*q), ncol = q)
#' X <- matrix(rnorm(n*p), ncol = p)
#' coords <- matrix(runif(n*2), ncol = 2)
#'
#' data <- list(Y = Y, X = X)
#' priors <- list(mu_B = matrix(0, nrow = p, ncol = q),
#'                              V_r = diag(10, p),
#'                              Psi = diag(1, q),
#'                              nu = 3)
#' hyperpar <- list(alpha = 0.5, phi = 1)
#' subset_size <- 200
#'
#' res <- spBPS(data, priors, coords, hyperpar, subset_size = subset_size)
#'
#' }
#'
#' @export
spBPS <- function(data,
                  priors,
                  coords,
                  hyperpar,
                  subset_size = 500L,
                  K = NULL,
                  cv_folds = 5L,
                  rp = 1,
                  combine_method = c("bps", "pseudoBMA"),
                  draws = 0L,
                  newdata = NULL,
                  include_latent = FALSE,
                  cores = NULL) {

  cat("\n====================================================\n")
  cat("         Welcome to spBPS Bayesian Engine\n")
  cat("====================================================\n\n")

  combine_method <- match.arg(combine_method)

  map_parallel <- function(indices, fn, exports = NULL, packages = NULL) {
    if (requireNamespace("foreach", quietly = TRUE)) {
      `%dopar%` <- foreach::`%dopar%`
      workers <- tryCatch(foreach::getDoParWorkers(), error = function(e) 1L)
      if (!is.null(workers) && workers > 1L) {
        export_vars <- unique(c(exports, ls(envir = parent.frame()), ls(envir = environment(fn))))
        return(foreach::foreach(i = indices,
                                 .export = export_vars,
                                 .packages = packages) %dopar% fn(i))
      }
    }
    lapply(indices, fn)
  }

  # optional internal backend
  if (!is.null(cores) && cores > 1L) {
    cat("Preparing parallel backend with", cores, "cores... \n\n")
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is required when 'cores' > 1.")
    }
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
  }

  Y <- data$Y
  if (is.null(dim(Y))) {
    Y <- matrix(Y, ncol = 1L)
  }
  X <- data$X
  n <- nrow(Y)

  if (is.null(K)) {
    K <- max(1L, min(n, as.integer(ceiling(n / subset_size))))
  }

  # Subset data once using the provided coordinates
  cat("Pritioning data into K =", K, "subsets ... \n\n")
  subsets <- spBPS::subset_data(data = list(Y = Y, X = X, crd = coords), K = K)
  # Local models: compute weights and epd per subset
  J <- nrow(expand_grid_cpp(hyperpar$alpha, hyperpar$phi))
  cat("Computing local stacking weights over J =", J,"models ...\n")
  fit_list <- map_parallel(seq_len(K), function(i) {
    Yi <- subsets$Y_list[[i]]
    Xi <- subsets$X_list[[i]]
    crd_i <- subsets$crd_list[[i]]
    p <- ncol(Xi)
    q <- ncol(Yi)

    local <- suppressWarnings(BPS_weights_MvT(
      data = list(Y = Yi, X = Xi),
      priors = priors,
      coords = crd_i,
      hyperpar = hyperpar,
      K = cv_folds
    ))

    list(local$epd, local$W)
  }, exports = c("subsets", "priors", "hyperpar", "cv_folds"), packages = "spBPS")
  cat("Local weights computed.\n\n")

  # Global combination
  cat("Computing global stacking weights over K =", K,"partitions ...\n")
  if (combine_method == "bps") {
    comb <- suppressWarnings(BPS_combine(fit_list = fit_list, K = K, rp = rp))
  } else {
    comb <- BPS_PseudoBMA(fit_list = fit_list)
  }

  W_global <- comb$W
  W_local <- comb$W_list
  cat("Global weights computed.\n\n")

  out <- list(
    subsets = subsets,
    weights_global = W_global,
    weights_local = W_local,
    epd = lapply(fit_list, `[[`, 1)
  )

  # Joint posterior + predictive draws via BPS_post_MvT (requires newdata when draws requested)
  draw_total <- draws
  if (draw_total > 0) {

    if (is.null(newdata)) {

      cat("Posterior Sampling (R =",draws, "draws) ...\n")
      subset_ind <- sample.int(K, draw_total, replace = TRUE, prob = W_global)
      draw_list <- map_parallel(seq_len(draw_total), function(r) {
        idx <- subset_ind[r]
        Yi <- subsets$Y_list[[idx]]
        Xi <- subsets$X_list[[idx]]
        crd_i <- subsets$crd_list[[idx]]
        res <- BPS_postdraws_MvT(
          data = list(Y = Yi, X = Xi),
          priors = priors,
          coords = crd_i,
          hyperpar = hyperpar,
          W = W_local[[idx]],
          R = 1,
          par = !isTRUE(include_latent)
        )
        list(post = res[[1]][[1]])
      }, exports = c("subsets", "priors", "hyperpar", "W_local", "newdata"), packages = "spBPS")

      out$posterior <- lapply(draw_list, `[[`, "post")
      cat("Posterior sampling completed.\n\n")

    } else {

    cat("Posterior & Predictive Sampling (R =",draws, "draws) ...\n")
    subset_ind <- sample.int(K, draw_total, replace = TRUE, prob = W_global)
    draw_list <- map_parallel(seq_len(draw_total), function(r) {
      idx <- subset_ind[r]
      Yi <- subsets$Y_list[[idx]]
      Xi <- subsets$X_list[[idx]]
      crd_i <- subsets$crd_list[[idx]]
      res <- BPS_post_MvT(
        data = list(Y = Yi, X = Xi),
        X_u = newdata$X,
        priors = priors,
        coords = crd_i,
        crd_u = newdata$coords,
        hyperpar = hyperpar,
        W = W_local[[idx]],
        R = 1
      )
      list(post = res$Post[[1]], pred = res$Pred[[1]])
    }, exports = c("subsets", "priors", "hyperpar", "W_local", "newdata"), packages = "spBPS")

    out$posterior <- lapply(draw_list, `[[`, "post")
    out$predictive <- lapply(draw_list, `[[`, "pred")
    cat("Posterior & Predictive sampling completed.\n\n")
  }

  }

  cat("====================================================\n")
  cat("     spBPS pipeline completed successfully!\n")
  cat("====================================================\n\n")

  return(out)
}
