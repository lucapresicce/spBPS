#ifndef SPBPS_RANDOM_UTILS_H
#define SPBPS_RANDOM_UTILS_H

#include <RcppArmadillo.h>
#include <random>
#include <atomic>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace spBPS {
namespace random {

// Seed the internal engine pool using R's RNG state (requires an active Rcpp::RNGScope).
void seed_from_R();

// Access a thread-local engine that is consistent with the most recent seed_from_R() call.
std::mt19937_64& engine();

// Draw helpers -------------------------------------------------------------

double standard_normal();

double chi_square(double df);

double gamma(double shape, double scale);

double beta(double alpha, double beta);

arma::vec standard_normal_vec(std::size_t n);

arma::vec gamma_vec(std::size_t n, double shape, double scale);

arma::vec dirichlet_sample(const arma::vec& alpha);

arma::uvec sample_indices(std::size_t size, std::size_t n_draws, const arma::vec& prob, bool replace = true);

arma::mat inverse_wishart(double df, const arma::mat& scale);

arma::mat matrix_normal(const arma::mat& mean,
                        const arma::mat& row_cov,
                        const arma::mat& col_cov,
                        bool row_cov_is_chol = false,
                        bool col_cov_is_chol = false);

arma::mat matrix_t(const arma::mat& mean,
                   const arma::mat& row_cov,
                   const arma::mat& col_cov,
                   double df);

arma::cube matrix_normal_draws(std::size_t n_draws,
                               const arma::mat& mean,
                               const arma::mat& row_cov,
                               const arma::mat& col_cov);

arma::cube matrix_t_draws(std::size_t n_draws,
                          const arma::mat& mean,
                          const arma::mat& row_cov,
                          const arma::mat& col_cov,
                          double df,
                          int threads = 0);

arma::vec multivariate_t(const arma::vec& mean,
                         const arma::mat& scale,
                         double df);

double dmvt(const arma::vec& x,
            const arma::vec& mean,
            const arma::mat& scale,
            double df,
            bool log_density = false);

double dmatrix_t(const arma::mat& X,
                 const arma::mat& mean,
                 const arma::mat& row_cov,
                 const arma::mat& col_cov,
                 double df,
                 bool log_density = false);

} // namespace random
} // namespace spBPS

#endif // SPBPS_RANDOM_UTILS_H
