#ifndef CODE_NEW_H
#define CODE_NEW_H
#include <RcppArmadillo.h>

Rcpp::List fit_cpp_MvT_v2(const Rcpp::List&, const Rcpp::List&,
                            const arma::mat&, const Rcpp::List&,
                            const arma::mat&);
Rcpp::List post_draws_MvT_v2(const Rcpp::List&, int, bool, int);
double d_pred_cpp_MvT_v2(const Rcpp::List&, const arma::mat&,
                          const arma::mat&, const arma::mat&,
                          const Rcpp::List&, const Rcpp::List&,
                          const arma::mat&);
arma::vec dens_kcv_MvT_v2(const Rcpp::List&, const Rcpp::List&,
                            const arma::mat&, const Rcpp::List&, int);
arma::mat models_dens_MvT_v2(const Rcpp::List&, const Rcpp::List&,
                               const arma::mat&, const Rcpp::List&, int);
Rcpp::List BPS_weights_MvT_v2(const Rcpp::List&, const Rcpp::List&,
                                const arma::mat&, const Rcpp::List&, int);
Rcpp::List BPS_combine_v2(const Rcpp::List&, int, double);
Rcpp::List compute_all_local_weights_v2(const Rcpp::List&, const Rcpp::List&,
                                         const Rcpp::List&, const Rcpp::List&,
                                         const Rcpp::List&, int, int);
// Full prediction (array output)
Rcpp::List BPS_post_MvT_v2(const Rcpp::List&, const arma::mat&,
                             const Rcpp::List&, const arma::mat&,
                             const arma::mat&, const Rcpp::List&,
                             const arma::vec&, int,
                             const arma::mat&,   // d_u
                             const Rcpp::List&,  // Rphi_u_list
                             int);               // n_cores
// Batch prediction
Rcpp::List BPS_post_MvT_v2_batch(const Rcpp::List&, const arma::mat&,
                                   const Rcpp::List&, const arma::mat&,
                                   const arma::mat&, const Rcpp::List&,
                                   const arma::vec&, int,
                                   const arma::mat&,
                                   const Rcpp::List&,
                                   int, int);
Rcpp::List BPS_postdraws_MvT_v2(const Rcpp::List&, const Rcpp::List&,
                                  const arma::mat&, const Rcpp::List&,
                                  const arma::vec&, int, bool);
double dmatrix_t_test(const arma::mat&, const arma::mat&, const arma::mat&,
                       const arma::mat&, double, bool);
#endif
