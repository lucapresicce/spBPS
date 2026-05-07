// [[Rcpp::plugins(openmp)]]
#include "random_utils.h"
#include <Rmath.h>
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
    if (seed == 0u) seed = 1u;
    base_seed.store(seed, std::memory_order_relaxed);
    counter.store(0u, std::memory_order_relaxed);
    epoch.fetch_add(1u, std::memory_order_acq_rel);
  }
  uint64_t next_seed() {
    uint64_t id = counter.fetch_add(1u, std::memory_order_relaxed) + 1u;
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
  uint64_t s1 = static_cast<uint64_t>(
      std::floor(u1 * static_cast<double>(std::numeric_limits<uint32_t>::max())));
  uint64_t s2 = static_cast<uint64_t>(
      std::floor(u2 * static_cast<double>(std::numeric_limits<uint32_t>::max())));
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

inline arma::mat std_normal_mat(arma::uword nr, arma::uword nc) {
  arma::mat Z(nr, nc);
  auto& eng  = thread_engine();
  auto& dist = normal_dist();
  for (arma::uword j = 0; j < nc; ++j)
    for (arma::uword i = 0; i < nr; ++i)
      Z(i, j) = dist(eng);
  return Z;
}

// log Gamma_q(alpha) — identical to mniw's logMultiGamma
inline double log_multi_gamma(double alpha, int q) {
  double lmg = 0.5 * static_cast<double>(q * (q - 1)) * M_LN_SQRT_PI;
  for (int i = 0; i < q; ++i)
    lmg += std::lgamma(alpha - 0.5 * static_cast<double>(i));
  return lmg;
}

} // anon namespace

// ---- public API -------------------------------------------------------

void seed_from_R() {
  EnginePool::instance().reseed(draw_seed_from_R());
}

std::mt19937_64& engine() { return thread_engine(); }

double standard_normal() { return normal_dist()(thread_engine()); }

double chi_square(double df) {
  if (df <= 0.0) Rcpp::stop("chi_square: df must be positive");
  std::gamma_distribution<double> dist(0.5 * df, 2.0);
  return dist(thread_engine());
}

double gamma(double shape, double scale) {
  if (shape <= 0.0 || scale <= 0.0)
    Rcpp::stop("gamma: invalid parameters");
  std::gamma_distribution<double> d(shape, scale);
  return d(thread_engine());
}

double beta(double alpha, double beta_par) {
  double x = gamma(alpha, 1.0);
  double y = gamma(beta_par, 1.0);
  return x / (x + y);
}

arma::vec standard_normal_vec(std::size_t n) {
  arma::vec out(n);
  auto& eng = thread_engine();
  auto& dist = normal_dist();
  for (std::size_t i = 0; i < n; ++i) out(i) = dist(eng);
  return out;
}

arma::vec gamma_vec(std::size_t n, double shape, double scale) {
  arma::vec out(n);
  for (std::size_t i = 0; i < n; ++i) out(i) = gamma(shape, scale);
  return out;
}

arma::vec dirichlet_sample(const arma::vec& alpha) {
  if (alpha.min() <= 0.0) Rcpp::stop("Dirichlet: all alpha must be positive");
  arma::vec draws(alpha.n_elem);
  double total = 0.0;
  for (arma::uword i = 0; i < alpha.n_elem; ++i) {
    double val = gamma(alpha(i), 1.0);
    draws(i) = val; total += val;
  }
  if (total == 0.0)
    return arma::ones<arma::vec>(alpha.n_elem) / static_cast<double>(alpha.n_elem);
  return draws / total;
}

arma::uvec sample_indices(std::size_t size, std::size_t n_draws,
                           const arma::vec& prob, bool replace) {
  if (!replace) Rcpp::stop("sample_indices: only replace=true supported");
  if (prob.n_elem != size) Rcpp::stop("sample_indices: prob length mismatch");
  std::vector<double> w(size);
  double sw = 0.0;
  for (std::size_t i = 0; i < size; ++i) {
    if (prob(i) < 0.0) Rcpp::stop("sample_indices: negative probability");
    w[i] = prob(i); sw += prob(i);
  }
  if (sw == 0.0) std::fill(w.begin(), w.end(), 1.0);
  std::discrete_distribution<std::size_t> dist(w.begin(), w.end());
  arma::uvec out(n_draws);
  auto& eng = thread_engine();
  for (std::size_t i = 0; i < n_draws; ++i) out(i) = dist(eng);
  return out;
}

arma::mat normal_matrix(arma::uword n_rows, arma::uword n_cols) {
  return std_normal_mat(n_rows, n_cols);
}

// ---- Inverse-Wishart ---------------------------------------------------
// X ~ IW(scale, df)  via Bartlett decomposition.
// Verified correct; only hardened with numerical fallbacks.
arma::mat inverse_wishart(double df, const arma::mat& scale) {
  arma::uword p = scale.n_rows;
  if (df <= static_cast<double>(p) - 1.0)
    Rcpp::stop("inverse_wishart: df too small");

  arma::mat Lpsi;
  if (!arma::chol(Lpsi, scale, "lower")) {
    arma::mat reg = scale;
    reg.diag() += 1e-10 * arma::trace(scale) / p;
    arma::chol(Lpsi, reg, "lower");
  }
  // L_Psi^{-1} (lower triangular)
  arma::mat Lpsi_inv = arma::solve(arma::trimatl(Lpsi),
                                    arma::eye<arma::mat>(p, p));

  // Bartlett matrix A (lower triangular)
  arma::mat A(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < p; ++i) {
    A(i, i) = std::sqrt(chi_square(df - static_cast<double>(i)));
    for (arma::uword j = i + 1; j < p; ++j)
      A(j, i) = standard_normal();
  }

  arma::mat L = Lpsi_inv * A;  // L*L' ~ W(Psi^{-1}, df)
  arma::mat W = L * L.t();
  arma::mat result;
  if (!arma::inv_sympd(result, W)) {
    W.diag() += 1e-10 * arma::trace(W) / p;
    arma::inv_sympd(result, W);
  }
  return result;
}

// ---- Matrix-Normal sampler ---------------------------------------------
// X ~ MN(mean, row_cov, col_cov)
// X = mean + Lr * Z * Lc'   where Lr/Lc are lower Cholesky factors.
// Pass row_cov_is_chol=true / col_cov_is_chol=true to reuse pre-computed
// factors and avoid redundant O(u^3) Cholesky decompositions.
arma::mat matrix_normal(const arma::mat& mean,
                        const arma::mat& row_cov,
                        const arma::mat& col_cov,
                        bool row_cov_is_chol,
                        bool col_cov_is_chol) {
  arma::mat Lr, Lc;

  if (row_cov_is_chol) {
    Lr = row_cov;
  } else {
    if (!arma::chol(Lr, row_cov, "lower")) {
      arma::mat reg = row_cov;
      reg.diag() += 1e-10 * arma::trace(row_cov) / row_cov.n_rows;
      arma::chol(Lr, reg, "lower");
    }
  }

  if (col_cov_is_chol) {
    Lc = col_cov;
  } else {
    if (!arma::chol(Lc, col_cov, "lower")) {
      arma::mat reg = col_cov;
      reg.diag() += 1e-10 * arma::trace(col_cov) / col_cov.n_rows;
      arma::chol(Lc, reg, "lower");
    }
  }

  arma::mat Z = std_normal_mat(mean.n_rows, mean.n_cols);
  return mean + Lr * Z * Lc.t();
}

// ---- Matrix-T sampler --------------------------------------------------
// X ~ MT(mean, row_cov, col_cov, df)
//
// BUG FIX: original had row_cov and col_cov swapped in the IW call,
// causing the large (u x u) matrix to be drawn from IW instead of the
// small (q x q) col_cov. Correct parameterization (matches mniw::rMT):
//   Sigma_c ~ IW(col_cov, df)          [q×q small matrix]
//   X | Sigma_c ~ MN(mean, row_cov, Sigma_c)
arma::mat matrix_t(const arma::mat& mean,
                   const arma::mat& row_cov,
                   const arma::mat& col_cov,
                   double df) {
  arma::mat sigma_c = inverse_wishart(df, col_cov);   // q×q only
  return matrix_normal(mean, row_cov, sigma_c, false, false);
}

// ---- Batch Matrix-Normal (OpenMP, chol computed once) ------------------
arma::cube matrix_normal_draws(std::size_t n_draws,
                               const arma::mat& mean,
                               const arma::mat& row_cov,
                               const arma::mat& col_cov) {
  arma::uword nr = mean.n_rows, nc = mean.n_cols;
  arma::mat Lr, Lc;
  arma::chol(Lr, row_cov, "lower");
  arma::chol(Lc, col_cov, "lower");

  arma::cube out(nr, nc, n_draws);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::int64_t idx = 0; idx < static_cast<std::int64_t>(n_draws); ++idx) {
    arma::mat Z = std_normal_mat(nr, nc);
    out.slice(static_cast<arma::uword>(idx)) = mean + Lr * Z * Lc.t();
  }
  return out;
}

// ---- Batch Matrix-T (OpenMP, chol(row_cov) computed once) --------------
// IW on col_cov per draw (small q×q). Row chol reused across all draws.
arma::cube matrix_t_draws(std::size_t n_draws,
                          const arma::mat& mean,
                          const arma::mat& row_cov,
                          const arma::mat& col_cov,
                          double df,
                          int threads) {
  if (df <= static_cast<double>(col_cov.n_rows) - 1.0)
    Rcpp::stop("matrix_t_draws: df too small");

  arma::uword nr = mean.n_rows, nc = mean.n_cols;
  arma::mat Lr;
  arma::chol(Lr, row_cov, "lower");   // computed ONCE — the big win

  arma::cube out(nr, nc, n_draws);
#ifdef _OPENMP
  int mt = (threads > 0) ? threads : omp_get_max_threads();
  mt = std::max(1, mt);
#pragma omp parallel for schedule(static) num_threads(mt)
#endif
  for (std::int64_t idx = 0; idx < static_cast<std::int64_t>(n_draws); ++idx) {
    arma::mat sigma_c = inverse_wishart(df, col_cov);    // q×q per draw
    arma::mat Lc;
    arma::chol(Lc, sigma_c, "lower");
    arma::mat Z = std_normal_mat(nr, nc);
    out.slice(static_cast<arma::uword>(idx)) = mean + Lr * Z * Lc.t();
  }
  return out;
}

// ---- Multivariate-t sampler (vector) -----------------------------------
arma::vec multivariate_t(const arma::vec& mean,
                         const arma::mat& scale,
                         double df) {
  arma::mat L;
  arma::chol(L, scale, "lower");
  arma::vec z = standard_normal_vec(mean.n_rows);
  double g = chi_square(df);
  return mean + (L * z) * std::sqrt(df / g);
}

// ---- Multivariate-t log-density (vector) -------------------------------
double dmvt(const arma::vec& x,
            const arma::vec& mean,
            const arma::mat& scale,
            double df,
            bool log_density) {
  arma::mat L;
  arma::chol(L, scale, "lower");
  double ldet = 2.0 * arma::sum(arma::log(L.diag()));
  arma::vec y  = arma::solve(arma::trimatl(L), x - mean);
  double quad  = arma::dot(y, y);
  double p     = static_cast<double>(mean.n_elem);
  double lnorm = std::lgamma(0.5*(df+p)) - std::lgamma(0.5*df)
                 - 0.5*(p*std::log(df*M_PI) + ldet);
  double lv    = lnorm - 0.5*(df+p)*std::log1p(quad/df);
  return log_density ? lv : std::exp(lv);
}

// ---- Matrix-T log-density ----------------------------------------------
// BUG FIX: original used a quadratic-form + log1p(quad/df) formula which
// is incorrect for the Matrix-T distribution. Correct formula uses the
// Sylvester identity for the log-determinant ratio, matching mniw exactly.
//
// For X ~ MT_p×q(Lambda, SigmaR, SigmaC, nu):
//   nuq  = nu + q - 1,   nupq = nuq + p
//
// Case p <= q:  B = SigmaR + diff * SigmaC^{-1} * diff' (p×p)
//               ldet_ratio = 0.5*(log|B| - log|SigmaR|)
// Case p > q:   B = SigmaC + diff' * SigmaR^{-1} * diff (q×q)
//               ldet_ratio = 0.5*(log|B| - log|SigmaC|)
//
// log p = -(nupq*ldet_ratio + q*ldSigmaR_half + p*ldSigmaC_half
//           + p*q*M_LN_SQRT_PI)
//         + logGamma_q(nupq/2) - logGamma_q(nuq/2)
double dmatrix_t(const arma::mat& X,
                 const arma::mat& mean,
                 const arma::mat& row_cov,
                 const arma::mat& col_cov,
                 double df,
                 bool log_density) {
  if (X.n_rows != mean.n_rows || X.n_cols != mean.n_cols)
    Rcpp::stop("dmatrix_t: dimension mismatch");

  arma::uword p = mean.n_rows;
  arma::uword q = mean.n_cols;
  double nuq  = df + static_cast<double>(q) - 1.0;
  double nupq = nuq + static_cast<double>(p);

  arma::mat Lr, Lc;
  arma::chol(Lr, row_cov, "lower");
  arma::chol(Lc, col_cov, "lower");

  // Half log-determinants (= sum(log(diag(L))) = 0.5*log|S|)
  double ldR = arma::sum(arma::log(Lr.diag()));
  double ldC = arma::sum(arma::log(Lc.diag()));

  arma::mat diff = X - mean;   // p×q

  double ldet_ratio;
  if (p <= q) {
    // A = Lc^{-1} * diff'  (q×p)
    // A'*A = diff * SigmaC^{-1} * diff'  (p×p)
    arma::mat A  = arma::solve(arma::trimatl(Lc), diff.t());
    arma::mat B  = row_cov + A.t() * A;
    arma::mat LB;
    arma::chol(LB, B, "lower");
    ldet_ratio = arma::sum(arma::log(LB.diag())) - ldR;
  } else {
    // A = Lr^{-1} * diff  (p×q)
    // A'*A = diff' * SigmaR^{-1} * diff  (q×q)
    arma::mat A  = arma::solve(arma::trimatl(Lr), diff);
    arma::mat B  = col_cov + A.t() * A;
    arma::mat LB;
    arma::chol(LB, B, "lower");
    ldet_ratio = arma::sum(arma::log(LB.diag())) - ldC;
  }

  double pq = static_cast<double>(p * q);
  double lp =
      -(nupq * ldet_ratio
        + static_cast<double>(q) * ldR
        + static_cast<double>(p) * ldC
        + pq * M_LN_SQRT_PI)
      + log_multi_gamma(0.5 * nupq, static_cast<int>(q))
      - log_multi_gamma(0.5 * nuq,  static_cast<int>(q));

  return log_density ? lp : std::exp(lp);
}

} // namespace random
} // namespace spBPS
