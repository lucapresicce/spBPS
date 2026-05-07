#' @keywords internal
#' @export
print.spBPS <- function(x, ...) {
  cat("spBPS object\n")
  cat(sprintf("  K subsets:   %d\n", length(x$subsets$Y_list)))
  cat(sprintf("  J models:    %d\n", length(x$weights_global)))
  if (!is.null(x$posterior))
    cat(sprintf("  Posterior draws:  %d\n", dim(x$posterior$beta)[3L]))
  if (!is.null(x$predictive))
    cat(sprintf("  Predictive draws: %d  (u = %d)\n",
                dim(x$predictive$Yu)[3L], dim(x$predictive$Yu)[1L]))
  invisible(x)
}

#' @keywords internal
#' @export
print.spBPS_pred <- function(x, ...) {
  cat(sprintf("spBPS_pred  [u=%d  q=%d  R=%d  %.1f s]\n",
              x$u, dim(x$Yu)[2L], x$draws, x$elapsed))
  invisible(x)
}

#' @keywords internal
#' @export
print.spBPS_pred_summary <- function(x, ...) {
  cat(sprintf("spBPS_pred_summary  [u=%d  q=%d  R=%d  probs=(%.3f,%.3f)  %.1f s]\n",
              x$u, ncol(x$mu_Yu), x$draws,
              x$probs[1L], x$probs[2L], x$elapsed))
  invisible(x)
}
