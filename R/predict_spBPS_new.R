#' Predict at new locations using a fitted spBPS model
#'
#' Generates posterior predictive samples at new spatial locations using a
#' fitted \code{spBPS} object. Memory usage is controlled via
#' \code{pred_batch_size}
#'
#' @param object Object of class \code{"spBPS"} returned by \code{spBPS()}.
#' @param newdata List with elements \code{X} (u x p covariate matrix) and \code{coords} (u x d coordinate matrix).
#' @param draws Integer. Number of posterior predictive draws. Default 200.
#' @param pred_batch_size Integer. Number of prediction sites processed per batch. Default 200. Smaller values use less RAM; larger values may be faster.
#' @param n_cores Integer. Number of OMP threads. Default 1.
#' @param return_summary Logical. If \code{FALSE} (default) return full draw arrays. If \code{TRUE} return summary statistics.
#' @param probs Numeric vector of length 2. Quantile probabilities for credible intervals when \code{return_summary = TRUE}. Default \code{c(0.025, 0.975)}.
#' @param ... Currently unused.
#'
#' @return When \code{return_summary = FALSE}, a list with arrays \code{Wu},
#'   \code{Yu}, \code{MY} of dimension u x q x R.
#'   When \code{return_summary = TRUE}, a list with u x q matrices
#'   \code{mu_Yu}, \code{sd_Yu}, \code{q_lo}, \code{q_hi}, \code{mu_Wu}.
#'
#' @seealso \code{\link{spBPS}}
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
#'
#' u <- 500
#' p <- 2
#' q <- 1
#'
#' X_u <- matrix(rnorm(u * p), ncol = p)
#' crd_u <- matrix(runif(u * 2), ncol = 2)
#'
#' # Predictive sampling
#' pred <- predict(res, newdata=list(X=X_u, coords=crd_u),
#'                 draws=200L, pred_batch_size=500L)
#'
#' # Summary statistics only (memory-efficient for large u)
#' pred_sum <- predict(res, newdata=list(X=X_u, coords=crd_u),
#'                     draws=500L, return_summary=TRUE,
#'                     probs=c(0.025, 0.975))
#' }
#'
#'
#' @method predict spBPS
#' @export
predict.spBPS <- function(object,
                          newdata,
                          draws           = 200L,
                          pred_batch_size = 200L,
                          n_cores         = 1L,
                          return_summary  = FALSE,
                          probs           = c(0.025, 0.975),
                          ...) {

  if (!inherits(object, "spBPS"))
    stop("object must be of class 'spBPS' (output of spBPS()).")
  if (is.null(object$subsets) || is.null(object$priors) || is.null(object$hyperpar))
    stop("object is missing $subsets, $priors, or $hyperpar. Re-run spBPS().")
  if (!is.list(newdata) || is.null(newdata$X) || is.null(newdata$coords))
    stop("newdata must be a list with elements $X and $coords.")
  if (length(probs) != 2L || probs[1L] >= probs[2L] ||
      any(probs < 0) || any(probs > 1))
    stop("probs must be a length-2 vector with 0 <= probs[1] < probs[2] <= 1.")

  newdata$X      <- as.matrix(newdata$X)
  newdata$coords <- as.matrix(newdata$coords)
  X_u   <- newdata$X
  crd_u <- newdata$coords
  u     <- nrow(X_u)
  q     <- ncol(object$subsets$Y_list[[1L]])
  K     <- length(object$subsets$Y_list)

  B_eff <- if (is.null(pred_batch_size) || pred_batch_size <= 0L)
              min(200L, u) else min(as.integer(pred_batch_size), u)

  W_global   <- object$weights_global
  W_loc_vecs <- lapply(object$weights_local, as.numeric)
  subset_ind <- sample.int(K, draws, replace = TRUE, prob = W_global) - 1L

  cat(sprintf("[predict.spBPS]  u = %d  draws = %d  cores = %d\n",
              u, draws, n_cores))

  t_pred <- system.time({
    res <- BPS_pred_streaming_v2(
      Y_list          = object$subsets$Y_list,
      X_list          = object$subsets$X_list,
      crd_list        = object$subsets$crd_list,
      W_local_list    = W_loc_vecs,
      subset_ind      = as.integer(subset_ind),
      crd_u           = crd_u,
      X_u             = X_u,
      priors          = object$priors,
      hyperpar        = object$hyperpar,
      pred_batch_size = B_eff,
      n_cores         = as.integer(n_cores),
      return_summary  = return_summary,
      prob_lo         = probs[1L],
      prob_hi         = probs[2L])
  })

  cat(sprintf("  Completed in %.1f s.\n", t_pred["elapsed"]))

  if (!return_summary) {
    structure(list(Wu = res$Wu, Yu = res$Yu, MY = res$MY,
                   u = u, draws = draws, batch_size = B_eff,
                   elapsed = t_pred["elapsed"]),
              class = "spBPS_pred")
  } else {
    structure(list(mu_Yu = res$mu_Yu, sd_Yu = res$sd_Yu,
                   q_lo  = res$q_lo,  q_hi  = res$q_hi,
                   mu_Wu = res$mu_Wu, probs = probs,
                   u = u, draws = draws, batch_size = B_eff,
                   elapsed = t_pred["elapsed"]),
              class = "spBPS_pred_summary")
  }
}
