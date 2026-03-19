// [[Rcpp::plugins(openmp)]]

#include "random_utils.h"

#include <limits>
#include <cmath>
#include <vector>

namespace spBPS {
namespace random {

namespace {

struct EnginePool {
  std::atomic<uint64_t> base_seed{5489u};
  std::atomic<uint64_t> counter{0u};
  std::atomic<uint64_t> epoch{0u};

  static EnginePool& instance() {
    static EnginePool pool;
    return pool;
  }

  void reseed(uint64_t seed) {
    if (seed == 0u) {
      seed = 1u;
    }
    base_seed.store(seed, std::memory_order_relaxed);
    counter.store(0u, std::memory_order_relaxed);
    epoch.fetch_add(1u, std::memory_order_acq_rel);
  }

  uint64_t next_seed() {
    uint64_t id = counter.fetch_add(1u, std::memory_order_relaxed) + 1u;
    // mix the counter with the base seed using golden ratio constant
    constexpr uint64_t phi = 0x9E3779B97F4A7C15ULL;
    return base_seed.load(std::memory_order_relaxed) + phi * id;
  }

  uint64_t current_epoch() const {
    return epoch.load(std::memory_order_acquire);
  }
};

inline uint64_t draw_seed_from_R() {
  double u1 = R::runif(0.0, 1.0);
  double u2 = R::runif(0.0, 1.0);
  uint64_t s1 = static_cast<uint64_t>(std::floor(u1 * std::numeric_limits<uint32_t>::max()));
  uint64_t s2 = static_cast<uint64_t>(std::floor(u2 * std::numeric_limits<uint32_t>::max()));
  return (s1 << 32) ^ s2 ^ 0x9E3779B97F4A7C15ULL;
}

inline std::mt19937_64& thread_engine() {
  thread_local std::mt19937_64 eng(EnginePool::instance().next_seed());
  thread_local uint64_t eng_epoch = EnginePool::instance().current_epoch();

  uint64_t pool_epoch = EnginePool::instance().current_epoch();
  if (eng_epoch != pool_epoch) {
    eng.seed(EnginePool::instance().next_seed());
    eng_epoch = pool_epoch;
  }
  return eng;
}

inline std::normal_distribution<double>& normal_dist() {
  thread_local std::normal_distribution<double> dist(0.0, 1.0);
  return dist;
}

inline std::gamma_distribution<double> gamma_dist(double shape, double scale) {
  return std::gamma_distribution<double>(shape, scale);
}

inline arma::mat standard_normal_mat(arma::uword n_rows, arma::uword n_cols) {
  arma::mat Z(n_rows, n_cols);
  for (arma::uword j = 0; j < n_cols; ++j) {
    for (arma::uword i = 0; i < n_rows; ++i) {
      Z(i, j) = standard_normal();
    }
  }
  return Z;
}

} // anonymous namespace

void seed_from_R() {
  EnginePool::instance().reseed(draw_seed_from_R());
}

std::mt19937_64& engine() {
  return thread_engine();
}

double standard_normal() {
  return normal_dist()(thread_engine());
}

double chi_square(double df) {
  if (df <= 0.0) {
    Rcpp::stop("Degrees of freedom must be positive in chi_square");
  }
  std::gamma_distribution<double> dist(0.5 * df, 2.0);
  return dist(thread_engine());
}

double gamma(double shape, double scale) {
  if (shape <= 0.0 || scale <= 0.0) {
    Rcpp::stop("Invalid parameters for gamma distribution");
  }
  return gamma_dist(shape, scale)(thread_engine());
}

double beta(double alpha, double beta_par) {
  double x = gamma(alpha, 1.0);
  double y = gamma(beta_par, 1.0);
  return x / (x + y);
}

arma::vec standard_normal_vec(std::size_t n) {
  arma::vec out(n);
  for (std::size_t i = 0; i < n; ++i) {
    out(i) = standard_normal();
  }
  return out;
}

arma::vec gamma_vec(std::size_t n, double shape, double scale) {
  arma::vec out(n);
  for (std::size_t i = 0; i < n; ++i) {
    out(i) = gamma(shape, scale);
  }
  return out;
}

arma::vec dirichlet_sample(const arma::vec& alpha) {
  if (alpha.min() <= 0.0) {
    Rcpp::stop("Dirichlet parameters must be positive");
  }
  arma::vec draws(alpha.n_elem);
  double total = 0.0;
  for (arma::uword i = 0; i < alpha.n_elem; ++i) {
    double val = gamma(alpha(i), 1.0);
    draws(i) = val;
    total += val;
  }
  if (total == 0.0) {
    return arma::ones<arma::vec>(alpha.n_elem) / static_cast<double>(alpha.n_elem);
  }
  return draws / total;
}

arma::uvec sample_indices(std::size_t size, std::size_t n_draws, const arma::vec& prob, bool replace) {
  if (!replace) {
    Rcpp::stop("sample_indices currently supports replace = true only");
  }
  if (prob.n_elem != size) {
    Rcpp::stop("Probability vector length does not match size");
  }
  std::vector<double> weights(size);
  double sum_weights = 0.0;
  for (std::size_t i = 0; i < size; ++i) {
    double w = prob(i);
    if (w < 0.0) {
      Rcpp::stop("Probabilities must be non-negative");
    }
    weights[i] = w;
    sum_weights += w;
  }
  if (sum_weights == 0.0) {
    std::fill(weights.begin(), weights.end(), 1.0);
  }
  std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
  arma::uvec out(n_draws);
  for (std::size_t i = 0; i < n_draws; ++i) {
    out(i) = dist(thread_engine());
  }
  return out;
}

arma::mat inverse_wishart(double df, const arma::mat& scale) {
  arma::uword p = scale.n_rows;
  if (scale.n_cols != p) {
    Rcpp::stop("Scale matrix for inverse Wishart must be square");
  }
  if (df <= static_cast<double>(p) - 1.0) {
    Rcpp::stop("Degrees of freedom too small for inverse Wishart");
  }

  arma::mat chol_scale = arma::chol(scale, "lower");
  arma::mat chol_inv_scale = arma::solve(arma::trimatl(chol_scale), arma::eye<arma::mat>(p, p));
  arma::mat A(p, p, arma::fill::zeros);

  for (arma::uword i = 0; i < p; ++i) {
    A(i, i) = std::sqrt(chi_square(df - static_cast<double>(i)));
    for (arma::uword j = i + 1; j < p; ++j) {
      A(j, i) = standard_normal();
    }
  }

  arma::mat L = chol_inv_scale * A;
  arma::mat wishart_sample = L * L.t();
  arma::mat result = arma::inv_sympd(wishart_sample);
  return result;
}

arma::mat matrix_normal(const arma::mat& mean,
                        const arma::mat& row_cov,
                        const arma::mat& col_cov,
                        bool row_cov_is_chol,
                        bool col_cov_is_chol) {
  arma::mat Lr = row_cov_is_chol ? row_cov : arma::chol(row_cov, "lower");
  arma::mat Lc = col_cov_is_chol ? col_cov : arma::chol(col_cov, "lower");

  arma::mat Z(mean.n_rows, mean.n_cols);
  for (arma::uword j = 0; j < mean.n_cols; ++j) {
    for (arma::uword i = 0; i < mean.n_rows; ++i) {
      Z(i, j) = standard_normal();
    }
  }

  arma::mat sample = Lr * Z * Lc.t();
  sample += mean;
  return sample;
}

arma::mat matrix_t(const arma::mat& mean,
                   const arma::mat& row_cov,
                   const arma::mat& col_cov,
                   double df) {
  arma::mat scale = inverse_wishart(df + static_cast<double>(row_cov.n_rows) - 1.0, row_cov);
  arma::mat centered = matrix_normal(arma::zeros<arma::mat>(mean.n_rows, mean.n_cols), scale, col_cov);
  centered += mean;
  return centered;
}

arma::cube matrix_normal_draws(std::size_t n_draws,
                               const arma::mat& mean,
                               const arma::mat& row_cov,
                               const arma::mat& col_cov) {
  arma::uword n_rows = mean.n_rows;
  arma::uword n_cols = mean.n_cols;
  arma::mat Lr = arma::chol(row_cov, "lower");
  arma::mat Lc = arma::chol(col_cov, "lower");

  arma::cube samples(n_rows, n_cols, n_draws);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::int64_t idx = 0; idx < static_cast<std::int64_t>(n_draws); ++idx) {
    arma::mat Z = standard_normal_mat(n_rows, n_cols);
    arma::mat sample = Lr * Z * Lc.t();
    sample += mean;
    samples.slice(static_cast<arma::uword>(idx)) = sample;
  }

  return samples;
}

arma::cube matrix_t_draws(std::size_t n_draws,
                          const arma::mat& mean,
                          const arma::mat& row_cov,
                          const arma::mat& col_cov,
                          double df,
                          int threads) {
  if (df <= static_cast<double>(col_cov.n_rows) - 1.0) {
    Rcpp::stop("Degrees of freedom too small for matrix t draws");
  }

  arma::uword n_rows = mean.n_rows;
  arma::uword n_cols = mean.n_cols;
  arma::mat Lr = arma::chol(row_cov, "lower");

  arma::cube samples(n_rows, n_cols, n_draws);

#ifdef _OPENMP
  int max_threads = (threads > 0) ? threads : omp_get_max_threads();
  if (max_threads <= 0) {
    max_threads = 1;
  }
#pragma omp parallel for schedule(static) num_threads(max_threads)
#endif
  for (std::int64_t idx = 0; idx < static_cast<std::int64_t>(n_draws); ++idx) {
    arma::mat sigma_c = inverse_wishart(df, col_cov);
    arma::mat Lc = arma::chol(sigma_c, "lower");
    arma::mat Z = standard_normal_mat(n_rows, n_cols);
    arma::mat sample = Lr * Z * Lc.t();
    sample += mean;
    samples.slice(static_cast<arma::uword>(idx)) = sample;
  }

  return samples;
}

arma::vec multivariate_t(const arma::vec& mean,
                         const arma::mat& scale,
                         double df) {
  arma::mat chol_scale = arma::chol(scale, "lower");
  arma::vec z(mean.n_rows);
  for (arma::uword i = 0; i < mean.n_rows; ++i) {
    z(i) = standard_normal();
  }
  double g = chi_square(df);
  double factor = std::sqrt(df / g);
  arma::vec sample = mean + (chol_scale * z) * factor;
  return sample;
}

double dmvt(const arma::vec& x,
            const arma::vec& mean,
            const arma::mat& scale,
            double df,
            bool log_density) {
  arma::vec diff = x - mean;
  arma::mat L = arma::chol(scale, "lower");
  double log_det = 2.0 * arma::sum(arma::log(L.diag()));
  arma::vec y = arma::solve(arma::trimatl(L), diff);
  double quad = arma::dot(y, y);
  double p = static_cast<double>(mean.n_elem);
  double log_norm = std::lgamma(0.5 * (df + p)) - std::lgamma(0.5 * df) - 0.5 * (p * std::log(df * M_PI) + log_det);
  double log_kernel = -0.5 * (df + p) * std::log1p(quad / df);
  double log_val = log_norm + log_kernel;
  if (log_density) {
    return log_val;
  }
  return std::exp(log_val);
}

double dmatrix_t(const arma::mat& X,
                 const arma::mat& mean,
                 const arma::mat& row_cov,
                 const arma::mat& col_cov,
                 double df,
                 bool log_density) {
  if (X.n_rows != mean.n_rows || X.n_cols != mean.n_cols) {
    Rcpp::stop("Matrix dimensions mismatch in dmatrix_t");
  }
  arma::uword n = mean.n_rows;
  arma::uword p = mean.n_cols;

  arma::mat Lr = arma::chol(row_cov, "lower");
  arma::mat Lc = arma::chol(col_cov, "lower");

  arma::mat diff = X - mean;
  arma::mat temp = arma::solve(arma::trimatu(Lc.t()), diff.t());
  arma::mat standardized = arma::solve(arma::trimatl(Lr), temp.t());
  double quad = arma::accu(arma::square(standardized));

  double log_det_row = 2.0 * arma::sum(arma::log(Lr.diag()));
  double log_det_col = 2.0 * arma::sum(arma::log(Lc.diag()));
  double dim = static_cast<double>(n) * static_cast<double>(p);

  double log_norm = std::lgamma(0.5 * (df + dim)) - std::lgamma(0.5 * df)
    - 0.5 * dim * std::log(df * M_PI)
    - 0.5 * static_cast<double>(p) * log_det_row
    - 0.5 * static_cast<double>(n) * log_det_col;

  double log_kernel = -0.5 * (df + dim) * std::log1p(quad / df);
  double log_val = log_norm + log_kernel;
  if (log_density) {
    return log_val;
  }
  return std::exp(log_val);
}

} // namespace random
} // namespace spBPS
