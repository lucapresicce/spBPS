#ifndef SPBPS_STACKING_SOLVER_H
#define SPBPS_STACKING_SOLVER_H

#include <RcppArmadillo.h>
#include <string>

namespace spBPS {
namespace stacking {

enum class SolverMethod {
  MirrorDescent,
  ProjectedGradient,
  Adam,
  CVXR
};

struct SolverConfig {
  double tolerance = 1e-8;
  int max_iterations = 1000;
  double step_size = 0.5;
  double min_step_size = 1e-6;
  double smoothing = 1e-12;
  double backtracking = 0.5;
  double learning_rate = 0.05;
  double beta1 = 0.9;
  double beta2 = 0.999;
  double epsilon = 1e-8;
  SolverMethod method = SolverMethod::MirrorDescent;
};

struct SolverResult {
  arma::vec weights;
  double objective;
  int iterations;
  bool converged;
  double gradient_norm;
  std::string solver;
};

SolverResult solve_weights(const arma::mat& scores,
                           const SolverConfig& config = SolverConfig());

} // namespace stacking
} // namespace spBPS

#endif // SPBPS_STACKING_SOLVER_H
