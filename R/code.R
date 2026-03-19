#' Solver for Bayesian Predictive Stacking of Predictive densities convex optimization problem
#'
#' @param scores [matrix] \eqn{N \times K} of expected predictive density evaluations for the K models considered
#'
#' @return W [matrix] of Bayesian Predictive Stacking weights for the K models considered
#'
#' @importFrom CVXR Variable Maximize Problem psolve status value
#'
conv_opt <- function(scores) {
  # library(CVXR, quietly = T)

  # set up minimization problem and solve it
  weights <- Variable( ncol(scores) )
  constraints <- list(weights >= 0, sum(weights) == 1)

  # the constraint for sum up to 1 with positive weights
  f <- Maximize( mean( log( scores %*% weights ) ) )
  problem <- Problem(f, constraints)

  # set the solver
  solver_to_use <- if (requireNamespace("ECOSolveR", quietly = TRUE)) {
    "ECOS_BB"
  } else {
    "SCS"  # CLARABEL, ECOS, OSQP
  }
  result <- psolve(problem, solver = solver_to_use)

  # return the weights
  W <- if(status( problem ) == "solver_error") {
    matrix(rep(1/ncol(scores), ncol(scores)))
  } else {
    value( weights )
  }
  return(W)
}
