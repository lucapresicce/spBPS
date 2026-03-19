// [[Rcpp::plugins(openmp)]]

#include "stacking_solver.h"
#include "code.h"

#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>

namespace spBPS {
namespace stacking {

namespace {

double projected_kkt_norm(const arma::vec& gradient,
                          const arma::vec& weights,
                          double zero_threshold) {
  double lambda = arma::dot(gradient, weights);
  arma::vec residual = gradient - lambda;
  for (arma::uword k = 0; k < residual.n_elem; ++k) {
    if (weights(k) <= zero_threshold && residual(k) < 0.0) {
      residual(k) = 0.0;
    }
  }
  return arma::norm(residual, "inf");
}

arma::vec compute_density(const arma::mat& scores,
                          const arma::vec& weights,
                          double smoothing) {
  arma::uword N = scores.n_rows;
  arma::uword K = scores.n_cols;
  arma::vec dens(N, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::uword i = 0; i < N; ++i) {
    double accum = 0.0;
    for (arma::uword k = 0; k < K; ++k) {
      double val = scores(i, k) > smoothing ? scores(i, k) : smoothing;
      accum += val * weights(k);
    }
    dens(i) = accum > smoothing ? accum : smoothing;
  }

  return dens;
}

double compute_objective(const arma::vec& dens) {
  double acc = 0.0;
  arma::uword N = dens.n_elem;
  for (arma::uword i = 0; i < N; ++i) {
    acc += std::log(dens(i));
  }
  return acc / static_cast<double>(N);
}

arma::vec compute_gradient(const arma::mat& scores,
                           const arma::vec& dens,
                           double smoothing) {
  arma::uword N = scores.n_rows;
  arma::uword K = scores.n_cols;
  arma::vec grad(K, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel
  {
    arma::vec local_grad(K, arma::fill::zeros);
#pragma omp for nowait schedule(static)
    for (arma::uword i = 0; i < N; ++i) {
      double inv_d = 1.0 / dens(i);
      for (arma::uword k = 0; k < K; ++k) {
        double val = scores(i, k) > smoothing ? scores(i, k) : smoothing;
        local_grad(k) += val * inv_d;
      }
    }
#pragma omp critical
    {
      grad += local_grad;
    }
  }
#else
  for (arma::uword i = 0; i < N; ++i) {
    double inv_d = 1.0 / dens(i);
    for (arma::uword k = 0; k < K; ++k) {
      double val = scores(i, k) > smoothing ? scores(i, k) : smoothing;
      grad(k) += val * inv_d;
    }
  }
#endif

  grad /= static_cast<double>(N);
  return grad;
}

arma::vec project_simplex(const arma::vec& v) {
  arma::uword K = v.n_elem;
  arma::vec u = arma::sort(v, "descend");
  double cumulative = 0.0;
  arma::uword rho = 0;
  for (arma::uword j = 0; j < K; ++j) {
    cumulative += u(j);
    double t = (cumulative - 1.0) / static_cast<double>(j + 1);
    if (u(j) > t) {
      rho = j + 1;
    } else {
      break;
    }
  }
  arma::vec w(K, arma::fill::zeros);
  if (rho == 0) {
    w.fill(1.0 / static_cast<double>(K));
    return w;
  }
  double tau = (arma::accu(u.subvec(0, rho - 1)) - 1.0) / static_cast<double>(rho);
  for (arma::uword i = 0; i < K; ++i) {
    double val = v(i) - tau;
    w(i) = val > 0.0 ? val : 0.0;
  }
  double sum_w = arma::sum(w);
  if (sum_w <= 0.0) {
    w.fill(1.0 / static_cast<double>(K));
  } else {
    w /= sum_w;
  }
  return w;
}

SolverResult solve_mirror(const arma::mat& scores,
                          const SolverConfig& config) {
  arma::uword K = scores.n_cols;
  SolverResult result;
  result.weights = arma::ones<arma::vec>(K) / static_cast<double>(K);
  result.iterations = 0;
  result.objective = -std::numeric_limits<double>::infinity();
  result.converged = false;
  result.solver = "mirror";
  result.gradient_norm = 0.0;

  double zero_threshold = std::max(1e-10, std::sqrt(std::max(config.tolerance, 1e-16)));
  arma::vec dens = compute_density(scores, result.weights, config.smoothing);
  double objective = compute_objective(dens);
  double step = config.step_size;

  for (int iter = 0; iter < config.max_iterations; ++iter) {
    arma::vec gradient = compute_gradient(scores, dens, config.smoothing);
    double kkt_norm = projected_kkt_norm(gradient, result.weights, zero_threshold);
    result.gradient_norm = kkt_norm;
    if (kkt_norm < config.tolerance) {
      result.converged = true;
      result.iterations = iter;
      break;
    }
    arma::vec previous_weights = result.weights;
    double previous_objective = objective;

    bool accepted = false;
    arma::vec candidate_weights(K);
    arma::vec candidate_dens;
    double candidate_objective = previous_objective;

    int backtrack = 0;
    while (!accepted && backtrack < 50) {
      arma::vec scaled_grad = gradient * step;
      double max_entry = scaled_grad.max();
      arma::vec exponent = arma::exp(scaled_grad - max_entry);
      candidate_weights = previous_weights % exponent;
      double sum_weights = arma::sum(candidate_weights);
      if (sum_weights <= 0.0 || !std::isfinite(sum_weights)) {
        step *= 0.5;
        ++backtrack;
        continue;
      }
      candidate_weights /= sum_weights;

      candidate_weights.elem(arma::find(candidate_weights < config.smoothing)).fill(0.0);
      double norm_sum = arma::sum(candidate_weights);
      if (norm_sum <= 0.0) {
        candidate_weights.fill(1.0 / static_cast<double>(K));
      } else {
        candidate_weights /= norm_sum;
      }

      candidate_dens = compute_density(scores, candidate_weights, config.smoothing);
      candidate_objective = compute_objective(candidate_dens);

      if (candidate_objective + 1e-10 >= previous_objective || step <= config.min_step_size) {
        accepted = true;
      } else {
        step *= 0.5;
        ++backtrack;
      }
    }

    result.weights = candidate_weights;
    dens = candidate_dens;
    objective = candidate_objective;
    result.iterations = iter + 1;

    double weight_change = arma::norm(result.weights - previous_weights, 1);
    double grad_norm = projected_kkt_norm(gradient, previous_weights, zero_threshold);
    result.gradient_norm = grad_norm;

    if (weight_change < config.tolerance || grad_norm < config.tolerance) {
      result.converged = true;
      break;
    }

    if (step < config.min_step_size) {
      break;
    }
  }

  result.objective = objective;
  arma::vec dens_final = compute_density(scores, result.weights, config.smoothing);
  arma::vec grad_final = compute_gradient(scores, dens_final, config.smoothing);
  result.gradient_norm = projected_kkt_norm(grad_final, result.weights, zero_threshold);
  if (!result.converged && result.gradient_norm < 10 * config.tolerance) {
    result.converged = true;
  }

  return result;
}

SolverResult solve_projected(const arma::mat& scores,
                             const SolverConfig& config) {
  arma::uword K = scores.n_cols;
  SolverResult result;
  result.weights = arma::ones<arma::vec>(K) / static_cast<double>(K);
  result.iterations = 0;
  result.objective = -std::numeric_limits<double>::infinity();
  result.converged = false;
  result.solver = "projected";
  result.gradient_norm = 0.0;

  double zero_threshold = std::max(1e-10, std::sqrt(std::max(config.tolerance, 1e-16)));
  arma::vec dens = compute_density(scores, result.weights, config.smoothing);
  double objective = compute_objective(dens);

  double base_lr = config.learning_rate;
  if (base_lr <= 0.0 || !std::isfinite(base_lr)) {
    base_lr = 0.5;
  }
  double current_lr = base_lr;
  double min_lr = config.min_step_size;
  if (min_lr <= 0.0 || !std::isfinite(min_lr)) {
    min_lr = 1e-6;
  }
  double shrink = config.backtracking > 0.0 && config.backtracking < 1.0
    ? config.backtracking
    : 0.5;

  for (int iter = 0; iter < config.max_iterations; ++iter) {
    arma::vec gradient = compute_gradient(scores, dens, config.smoothing);
    double kkt_norm_current = projected_kkt_norm(gradient, result.weights, zero_threshold);
    result.gradient_norm = kkt_norm_current;
    if (kkt_norm_current < config.tolerance) {
      result.converged = true;
      result.iterations = iter;
      break;
    }
    arma::vec previous_weights = result.weights;
    double previous_objective = objective;

    double step = current_lr;
    arma::vec candidate_weights = project_simplex(previous_weights + step * gradient);
    arma::vec candidate_dens = compute_density(scores, candidate_weights, config.smoothing);
    double candidate_objective = compute_objective(candidate_dens);

    int backtrack = 0;
    while (candidate_objective + 1e-10 < previous_objective && step > min_lr && backtrack < 50) {
      step *= shrink;
      candidate_weights = project_simplex(previous_weights + step * gradient);
      candidate_dens = compute_density(scores, candidate_weights, config.smoothing);
      candidate_objective = compute_objective(candidate_dens);
      ++backtrack;
    }

    if (candidate_objective + 1e-10 < previous_objective && step <= min_lr) {
      step = min_lr;
      candidate_weights = project_simplex(previous_weights + step * gradient);
      candidate_dens = compute_density(scores, candidate_weights, config.smoothing);
      candidate_objective = compute_objective(candidate_dens);
    }

    result.weights = candidate_weights;
    dens = candidate_dens;
    objective = candidate_objective;
    result.iterations = iter + 1;

    double weight_change = arma::norm(result.weights - previous_weights, 1);
    double grad_norm = projected_kkt_norm(gradient, previous_weights, zero_threshold);
    result.gradient_norm = grad_norm;

    if (weight_change < config.tolerance || grad_norm < config.tolerance) {
      result.converged = true;
      break;
    }

    if (step <= min_lr && backtrack > 0 && candidate_objective + 1e-10 < previous_objective) {
      break;
    }

    if (backtrack == 0) {
      step = std::min(base_lr, step * 1.05);
    }
    current_lr = step;
  }

  arma::vec dens_final = compute_density(scores, result.weights, config.smoothing);
  arma::vec grad_final = compute_gradient(scores, dens_final, config.smoothing);
  result.gradient_norm = projected_kkt_norm(grad_final, result.weights, zero_threshold);
  result.objective = compute_objective(dens_final);
  if (!result.converged && result.gradient_norm < 10 * config.tolerance) {
    result.converged = true;
  }

  return result;
}

SolverResult solve_adam(const arma::mat& scores,
                        const SolverConfig& config) {
  arma::uword K = scores.n_cols;
  SolverResult result;
  result.weights = arma::ones<arma::vec>(K) / static_cast<double>(K);
  result.iterations = 0;
  result.objective = -std::numeric_limits<double>::infinity();
  result.converged = false;
  result.solver = "adam";
  result.gradient_norm = 0.0;

  arma::vec m(K, arma::fill::zeros);
  arma::vec v(K, arma::fill::zeros);
  double beta1 = config.beta1;
  double beta2 = config.beta2;
  double epsilon = config.epsilon;
  double lr = config.learning_rate;

  arma::vec dens = compute_density(scores, result.weights, config.smoothing);
  double objective = compute_objective(dens);

  for (int t = 1; t <= config.max_iterations; ++t) {
    arma::vec gradient = compute_gradient(scores, dens, config.smoothing);
    arma::vec previous_weights = result.weights;
    double previous_objective = objective;

    m = beta1 * m + (1.0 - beta1) * gradient;
    v = beta2 * v + (1.0 - beta2) * arma::square(gradient);
    arma::vec m_hat = m / (1.0 - std::pow(beta1, t));
    arma::vec v_hat = v / (1.0 - std::pow(beta2, t));

    arma::vec candidate_weights = project_simplex(previous_weights + lr * m_hat / (arma::sqrt(v_hat) + epsilon));
    arma::vec candidate_dens = compute_density(scores, candidate_weights, config.smoothing);
    double candidate_objective = compute_objective(candidate_dens);

    int backtrack = 0;
    while (candidate_objective + 1e-10 < previous_objective && backtrack < 20) {
      lr *= 0.5;
      candidate_weights = project_simplex(previous_weights + lr * m_hat / (arma::sqrt(v_hat) + epsilon));
      candidate_dens = compute_density(scores, candidate_weights, config.smoothing);
      candidate_objective = compute_objective(candidate_dens);
      ++backtrack;
    }

    result.weights = candidate_weights;
    dens = candidate_dens;
    objective = candidate_objective;
    result.iterations = t;

    double weight_change = arma::norm(result.weights - previous_weights, 1);
    double grad_norm = arma::norm(gradient, "inf");
    result.gradient_norm = grad_norm;

    if (weight_change < config.tolerance || grad_norm * lr < config.tolerance) {
      result.converged = true;
      break;
    }
  }

  arma::vec dens_final = compute_density(scores, result.weights, config.smoothing);
  arma::vec grad_final = compute_gradient(scores, dens_final, config.smoothing);
  result.gradient_norm = arma::norm(grad_final, "inf");
  result.objective = compute_objective(dens_final);
  if (!result.converged && result.gradient_norm < 10 * config.tolerance) {
    result.converged = true;
  }

  return result;
}

SolverResult solve_cvxr(const arma::mat& scores,
                        const SolverConfig& config) {
  SolverResult result;
  result.iterations = 0;
  result.objective = -std::numeric_limits<double>::infinity();
  result.converged = true;
  result.solver = "cvxr";
  result.gradient_norm = 0.0;

  Rcpp::NumericMatrix solver_out = CVXR_opt(scores);
  arma::mat weights_mat = Rcpp::as<arma::mat>(solver_out);
  if (weights_mat.n_elem == 0) {
    Rcpp::stop("CVXR solver returned empty weights");
  }
  arma::vec weights = arma::vectorise(weights_mat);
  weights.elem(arma::find(weights < 0)).zeros();
  double sum_w = arma::sum(weights);
  if (sum_w <= 0.0) {
    weights.fill(1.0 / static_cast<double>(weights.n_elem));
  } else {
    weights /= sum_w;
  }
  result.weights = weights;
  arma::vec dens = compute_density(scores, result.weights, config.smoothing);
  result.objective = compute_objective(dens);
  arma::vec grad = compute_gradient(scores, dens, config.smoothing);
  result.gradient_norm = arma::norm(grad, "inf");
  return result;
}

} // anonymous namespace

SolverResult solve_weights(const arma::mat& scores,
                           const SolverConfig& config) {
  if (scores.n_rows == 0 || scores.n_cols == 0) {
    Rcpp::stop("Scores matrix must have positive dimensions");
  }

  switch (config.method) {
  case SolverMethod::MirrorDescent:
    return solve_mirror(scores, config);
  case SolverMethod::ProjectedGradient: {
    SolverResult res = solve_projected(scores, config);
    if (!res.converged) {
      SolverConfig backup = config;
      backup.method = SolverMethod::MirrorDescent;
      SolverResult fallback = solve_mirror(scores, backup);
      if (fallback.objective >= res.objective) {
        res = fallback;
        res.solver = std::string("projected+mirror");
      }
    }
    return res;
  }
  case SolverMethod::Adam:
    return solve_adam(scores, config);
  case SolverMethod::CVXR:
    return solve_cvxr(scores, config);
  default:
    Rcpp::stop("Unknown solver method");
  }
}

} // namespace stacking
} // namespace spBPS
