#' Unified spatial BPS workflow (multivariate path, works for q = 1)
#'
#' @description
#' Orchestrates subsetting, local stacking weight estimation, global stacking combination, and optional posterior or predictive simulation using Double Bayesian Predictive Stacking for latent spatial regression. Works for both multivariate outcomes and the univariate case via `q = 1`.
#'
#' @param data List with matrices `Y` (response, n x q) and `X` (covariates, n x p).
#' @param priors List of priors for the multivariate model: `mu_B` (p x q mean
#'   matrix), `V_r` (p x p covariance), `Psi` (q x q scale), `nu` (degrees of
#'   freedom).
#' @param coords Matrix of observation coordinates (n x d).
#' @param hyperpar List with elements `alpha` and `phi` (vectors allowed); defines
#'   the grid of models over which stacking weights are computed.
#' @param subset_size Target subset size when `K` is not provided. Default 500.
#' @param K Optional number of subsets. When `NULL`, computed as
#'   `ceiling(nrow(Y) / subset_size)` and lower-bounded at 1.
#' @param cv_folds Number of folds for local cross-validation. Default 5.
#' @param rp Fraction of rows used when recomputing global stacking weights
#'   (passed to `BPS_combine`). Ignored when `combine_method = "pseudoBMA"`.
#' @param combine_method Global combination method: Bayesian Predictive Stacking
#'   (`"bps"`) or pseudo-BMA (`"pseudoBMA"`).
#' @param draws Number of posterior/predictive draws to return (0 to skip
#'   sampling). When positive and `newdata` is supplied, joint posterior and
#'   predictive draws are returned. When positive and `newdata` is `NULL`,
#'   only posterior draws (beta, sigma) are returned.
#' @param newdata Optional list with `X` (u x p) and `coords` (u x d) for
#'   prediction locations. Required when `draws > 0` and predictions are desired.
#'   For large `u`, use `pred_batch_size` to control memory usage.
#' @param include_latent Unused; kept for compatibility.
#' @param n_cores Number of cores for parallel computation. Controls both the
#'   local weight estimation (parallel over K subsets) and predictive sampling.
#'   Default 1.
#' @param pred_batch_size Batch size for streaming prediction. Controls the
#'   maximum number of prediction sites processed at once, hence the peak memory
#'   usage. Default 200. Rule of thumb: peak RAM (MB) ? batch_size x q x draws x 8 / 1e6.
#'   Set `NULL` for automatic selection (min(200, u)).
#'
#' @return Object of class `"spBPS"` ? a list with components:
#'   \describe{
#'     \item{subsets}{Partition information: Y_list, X_list, crd_list.}
#'     \item{weights_global}{K-vector of global stacking weights.}
#'     \item{weights_local}{K-list of local stacking weight vectors.}
#'     \item{epd}{K-list of n_k x J log-density matrices.}
#'     \item{priors, hyperpar}{Stored for use by \code{\link{predict.spBPS}}.}
#'     \item{timings}{Named numeric vector: fitting, combination, sampling (seconds).}
#'     \item{posterior}{(if draws > 0) List with three-dimensional arrays:
#'       \code{beta} (p x q x R), \code{sigma} (q x q x R), \code{model} (R).}
#'     \item{predictive}{(if draws > 0 and newdata supplied) List with
#'       three-dimensional arrays: \code{Wu}, \code{Yu}, \code{MY} (each u x q x R).}
#'   }
#'
#' @seealso \code{\link{predict.spBPS}} for generating predictions on new data
#'   after fitting.
#'
#' @examples
#' \donttest{
#' n <- 1000
#' p <- 2
#' q <- 1
#'
#' Y <- matrix(rnorm(n * q), ncol = q)
#' X <- matrix(rnorm(n * p), ncol = p)
#' coords <- matrix(runif(n * 2), ncol = 2)
#'
#' data    <- list(Y = Y, X = X)
#' priors  <- list(mu_B = matrix(0, nrow = p, ncol = q),
#'                 V_r  = diag(10, p),
#'                 Psi  = diag(1, q),
#'                 nu   = 3)
#' hyperpar <- list(alpha = 0.5, phi = 1)
#'
#' res <- spBPS(data, priors, coords, hyperpar, subset_size = 200)
#' }
#'
#' @export
spBPS <- function(data,
                  priors,
                  coords,
                  hyperpar,
                  subset_size     = 500L,
                  K               = NULL,
                  cv_folds        = 5L,
                  rp              = 1.0,
                  combine_method  = c("bps", "pseudoBMA"),
                  draws           = 0L,
                  newdata         = NULL,
                  include_latent  = FALSE,
                  n_cores         = 1L,
                  pred_batch_size = 200L) {

  cat("\n====================================================\n")
  cat("         Welcome to spBPS Bayesian Engine\n")
  cat("====================================================\n\n")

  combine_method <- match.arg(combine_method)
  timings <- c(fitting = NA_real_, combination = NA_real_, sampling = NA_real_)

  # -- Normalise inputs ------------------------------------------------------
  data$Y <- as.matrix(data$Y)
  data$X <- as.matrix(data$X)
  coords  <- as.matrix(coords)
  n <- nrow(data$Y)

  if (is.null(K))
    K <- max(1L, min(n, as.integer(ceiling(n / subset_size))))

  Grid <- expand_grid_cpp(hyperpar$alpha, hyperpar$phi)
  J    <- nrow(Grid)

  # -- Step 1: Partition -----------------------------------------------------
  cat("Partitioning data into K =", K, "subsets ...\n\n")
  subsets <- subset_data(list(Y = data$Y, X = data$X, crd = coords), K = K)

  # -- Step 2: Local weight estimation --------------------------------------
  cat("Computing local stacking weights over J =", J, "models ...\n")
  t_fit <- system.time({
    fit_list <- compute_all_local_weights_v2(
      Y_list   = subsets$Y_list,
      X_list   = subsets$X_list,
      crd_list = subsets$crd_list,
      priors   = priors,
      hyperpar = hyperpar,
      cv_folds = cv_folds,
      n_cores  = n_cores)
  })
  timings["fitting"] <- t_fit["elapsed"]
  cat(sprintf("Local weights computed in %.1f s.\n\n", t_fit["elapsed"]))

  # -- Step 3: Global combination --------------------------------------------
  cat("Computing global stacking weights over K =", K, "partitions ...\n")
  t_comb <- system.time({
    if (combine_method == "bps") {
      comb <- BPS_combine_v2(fit_list = fit_list, K = K, rp = rp)
    } else {
      comb <- BPS_PseudoBMA_v2(fit_list = fit_list, n_cores = n_cores)
    }
  })
  timings["combination"] <- t_comb["elapsed"]
  W_global <- as.numeric(comb$W)
  W_local  <- comb$W_list
  cat(sprintf("Global weights computed in %.1f s.\n\n", t_comb["elapsed"]))

  out <- list(
    subsets        = subsets,
    weights_global = W_global,
    weights_local  = W_local,
    epd            = lapply(fit_list, `[[`, 1L),
    priors         = priors,
    hyperpar       = hyperpar,
    timings        = timings
  )
  class(out) <- "spBPS"

  # -- Step 4: Posterior / predictive sampling -------------------------------
  if (draws > 0L) {

    K_loc      <- length(subsets$Y_list)
    W_loc_vecs <- lapply(W_local, as.numeric)
    subset_ind <- sample.int(K_loc, draws, replace = TRUE, prob = W_global) - 1L

    if (is.null(newdata)) {
      # Posterior-only draws
      cat("Posterior sampling (R =", draws, "draws) ...\n")
      t_samp <- system.time({
        res_post <- BPS_post_all_subsets_v2(
          Y_list       = subsets$Y_list,
          X_list       = subsets$X_list,
          crd_list     = subsets$crd_list,
          W_local_list = W_loc_vecs,
          subset_ind   = as.integer(subset_ind),
          priors       = priors,
          hyperpar     = hyperpar,
          n_cores      = n_cores)
      })
      timings["sampling"] <- t_samp["elapsed"]
      out$posterior <- list(beta  = res_post$beta,
                            sigma = res_post$sigma,
                            model = res_post$model)
      cat(sprintf("Posterior sampling completed in %.1f s.\n\n", t_samp["elapsed"]))

    } else {
      # Posterior + predictive draws
      newdata$X      <- as.matrix(newdata$X)
      newdata$coords <- as.matrix(newdata$coords)
      X_u   <- newdata$X
      crd_u <- newdata$coords
      u     <- nrow(X_u)

      B_eff <- if (is.null(pred_batch_size) || pred_batch_size <= 0L)
                 min(200L, u)
               else
                 min(as.integer(pred_batch_size), u)

      cat(sprintf("Posterior & Predictive Sampling (R = %d draws) ...\n", draws))
      t_samp <- system.time({
        res_all <- BPS_pred_streaming_v2(
          Y_list          = subsets$Y_list,
          X_list          = subsets$X_list,
          crd_list        = subsets$crd_list,
          W_local_list    = W_loc_vecs,
          subset_ind      = as.integer(subset_ind),
          crd_u           = crd_u,
          X_u             = X_u,
          priors          = priors,
          hyperpar        = hyperpar,
          pred_batch_size = B_eff,
          n_cores         = as.integer(n_cores),
          return_summary  = FALSE,
          prob_lo         = 0.025,
          prob_hi         = 0.975)

        res_post <- BPS_post_all_subsets_v2(
          Y_list       = subsets$Y_list,
          X_list       = subsets$X_list,
          crd_list     = subsets$crd_list,
          W_local_list = W_loc_vecs,
          subset_ind   = as.integer(subset_ind),
          priors       = priors,
          hyperpar     = hyperpar,
          n_cores      = n_cores)
      })
      timings["sampling"] <- t_samp["elapsed"]

      out$posterior  <- list(beta  = res_post$beta,
                             sigma = res_post$sigma,
                             model = res_post$model)
      out$predictive <- list(Wu = res_all$Wu,
                             Yu = res_all$Yu,
                             MY = res_all$MY)
      cat(sprintf("Posterior & Predictive sampling completed in %.1f s.\n\n",
                  t_samp["elapsed"]))
    }
  }

  out$timings <- timings

  cat("====================================================\n")
  cat("     spBPS pipeline completed successfully!\n")
  cat("====================================================\n")
  cat(sprintf("  Fitting:     %.1f s\n", timings["fitting"]))
  cat(sprintf("  Combination: %.1f s\n", timings["combination"]))
  if (!is.na(timings["sampling"]))
    cat(sprintf("  Sampling:    %.1f s\n", timings["sampling"]))
  cat(sprintf("  TOTAL:       %.1f s\n\n",
              sum(timings, na.rm = TRUE)))
  return(out)
}

