// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include "code_new.h"
#include "utilsC.h"
#include "random_utils.h"
#include "stacking_solver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// ─────────────────────────────────────────────────────────────────────────────
// spBPS_par_for: sequential-bypass parallel dispatcher
//
// On Windows/GOMP, "#pragma omp parallel for" with num_threads(1) still
// invokes WaitForSingleObject barriers per iteration, adding ~5s overhead
// each time. With n_cores==1 we use a plain C++ for-loop (zero overhead).
// With n_cores>1 we use the OMP parallel for.
// ─────────────────────────────────────────────────────────────────────────────
template<class Fn>
static void spBPS_par_for(int K, int n_cores, Fn fn) {
    if (n_cores <= 1) {
        for (int k = 0; k < K; ++k) fn(k);
        return;
    }
#ifdef _OPENMP
    omp_set_max_active_levels(1);
#pragma omp parallel for schedule(dynamic) num_threads(n_cores)
    for (int k = 0; k < K; ++k) fn(k);
#else
    for (int k = 0; k < K; ++k) fn(k);
#endif
}


// ═══════════════════════════════════════════════════════════════════════════
// INTERNAL STRUCTURES
// ═══════════════════════════════════════════════════════════════════════════

struct FitResult_v2 {
    arma::mat V_star;
    arma::mat mu_star;
    arma::mat Psi_star;
    double    nu_star;
    arma::mat iRphi_s;
    // Cholesky factors are NOT stored here — computed lazily in BPS_post_MvT_v2
    // (once per model j) to avoid K×J×cv_folds unnecessary decompositions.
};

// ═══════════════════════════════════════════════════════════════════════════
// INTERNAL HELPERS  (no R API — safe inside OMP parallel for)
// ═══════════════════════════════════════════════════════════════════════════

// 1×n_tr row of distances from one test point to all training sites
static arma::rowvec cross_dist_row(const arma::rowvec& crd_test,
                                    const arma::mat&    crd_train) {
    arma::uword n = crd_train.n_rows;
    arma::rowvec d(n);
    for (arma::uword i = 0; i < n; ++i)
        d(i) = arma::norm(crd_test - crd_train.row(i), 2);
    return d;
}

// u×n_tr cross-distance matrix
static arma::mat cross_dist_mat(const arma::mat& crd_u,
                                const arma::mat& crd_train) {
    int u = crd_u.n_rows, n = crd_train.n_rows;
    arma::mat D(u, n);
    for (int i = 0; i < u; ++i)
        for (int j = 0; j < n; ++j)
            D(i,j) = arma::norm(crd_u.row(i) - crd_train.row(j), 2);
    return D;
}

// Symmetrize + Cholesky with jitter fallback
static arma::mat safe_chol(const arma::mat& M, const char* which = "lower") {
    arma::mat S = 0.5 * (M + M.t());
    arma::mat L;
    if (!arma::chol(L, S, which)) {
        S.diag() += 1e-8 * arma::trace(S) / S.n_rows;
        if (!arma::chol(L, S, which)) {
            S.diag() += 1e-6;
            arma::chol(L, S, which);
        }
    }
    return L;
}

// Pure fitting — no R API, OMP-safe.
// iV_r is passed pre-computed (shared across all K×J×cv calls).
static FitResult_v2 fit_pure_v2(const arma::mat& Y,
                                  const arma::mat& X,
                                  const arma::mat& coords,
                                  const arma::mat& mu_B,
                                  const arma::mat& iV_r,
                                  const arma::mat& Psi,
                                  double nu, double alpha, double phi) {
    int n = Y.n_rows;
    arma::mat Rphi_s = arma::exp(-phi * arma_dist(coords));

    // arma::inv (LU) matches original fit_cpp_MvT exactly.
    // For n_k < 500 LU is as fast or faster than inv_sympd.
    arma::mat iR_s    = arma::inv(Rphi_s);
    double    a       = alpha / (1.0 - alpha);
    arma::mat tX      = X.t();
    arma::mat V_BW    = a * tX;
    arma::mat iV_star = arma::join_vert(
        arma::join_horiz(a * tX * X + iV_r, V_BW),
        arma::join_horiz(V_BW.t(), iR_s + a * arma::eye<arma::mat>(n,n)));
    arma::mat V_star  = arma::inv(iV_star);
    arma::mat mu_star = V_star * arma::join_vert(a * tX * Y + iV_r * mu_B, a * Y);
    arma::mat Psi_star = Psi + a * Y.t() * Y
                             + mu_B.t() * iV_r * mu_B
                             - mu_star.t() * iV_star * mu_star;
    // Symmetrize once — prevents chol() warnings downstream
    V_star   = 0.5 * (V_star   + V_star.t());
    Psi_star = 0.5 * (Psi_star + Psi_star.t());

    FitResult_v2 r;
    r.V_star = V_star; r.mu_star = mu_star; r.Psi_star = Psi_star;
    r.nu_star = nu + n; r.iRphi_s = iR_s;
    return r;
}

// Predictive density at one test point — no R API, OMP-safe
static double d_pred_pure_v2(const arma::rowvec& Y_i,
                               const arma::rowvec& X_i,
                               const arma::rowvec& cross_d,
                               const FitResult_v2& fit,
                               double alpha, double phi) {
    arma::rowvec Rphi_us = arma::exp(-phi * cross_d);
    arma::rowvec M_u     = Rphi_us * fit.iRphi_s;
    arma::rowvec Mgamma  = arma::join_horiz(X_i, M_u);
    arma::mat    mu_til  = Mgamma * fit.mu_star;
    double V_wu = 1.0 - arma::as_scalar(M_u * Rphi_us.t());
    arma::mat M_til(1,1);
    M_til(0,0) = arma::as_scalar(Mgamma * fit.V_star * Mgamma.t())
                 + V_wu + (1.0/alpha - 1.0);
    arma::mat Y_m(1, Y_i.n_elem); Y_m.row(0) = Y_i;
    return spBPS::random::dmatrix_t(Y_m, mu_til, M_til, fit.Psi_star,
                                     fit.nu_star, false);
}

// K-fold CV densities — no R API, OMP-safe
static arma::vec dens_kcv_pure_v2(const arma::mat& Y, const arma::mat& X,
                                   const arma::mat& coords,
                                   const arma::mat& mu_B, const arma::mat& iV_r,
                                   const arma::mat& Psi, double nu,
                                   double alpha, double phi, int K_folds) {
    int n = Y.n_rows;
    arma::vec pred(n, arma::fill::zeros);
    arma::vec p_eq = arma::ones<arma::vec>(K_folds) / K_folds;
    arma::uvec fold = spBPS::random::sample_indices(K_folds, n, p_eq, true) + 1;

    for (int k = 1; k <= K_folds; ++k) {
        arma::uvec te = arma::find(fold == k), tr = arma::find(fold != k);
        FitResult_v2 fit = fit_pure_v2(Y.rows(tr), X.rows(tr), coords.rows(tr),
                                        mu_B, iV_r, Psi, nu, alpha, phi);
        arma::mat Yte = Y.rows(te), Xte = X.rows(te), cte = coords.rows(te);
        arma::mat ctr = coords.rows(tr);
        for (arma::uword i = 0; i < te.n_elem; ++i) {
            arma::rowvec cd = cross_dist_row(cte.row(i), ctr);
            pred(te(i)) = d_pred_pure_v2(Yte.row(i), Xte.row(i), cd,
                                          fit, alpha, phi);
        }
    }
    return pred;
}

// ═══════════════════════════════════════════════════════════════════════════
// CORE PREDICTION KERNEL — shared by v2 and batch variants
// Draws R_j samples for one model j, writes into pre-allocated cubes.
// Called from BPS_post_MvT_v2 (full u) and BPS_post_MvT_v2_batch (u_batch).
// ═══════════════════════════════════════════════════════════════════════════
static void draw_group_j(
        // model params
        const FitResult_v2& fit,
        double alpha_j, double a_j,
        // spatial matrices (already computed for this j)
        const arma::mat& M_u,       // u×n  kriging weights
        const arma::mat& L_Vu,      // lower chol of V_u (u×u)
        const arma::mat& L_Vstar_j, // lower chol of V_star
        const arma::mat& M_gamma,   // u×(p+n) [X_u | M_u]
        // X_u for mean
        const arma::mat& X_u,       // u×p
        // draw indices in global arrays
        const arma::uvec& draws_j,
        int p, int n_subset,
        // OMP thread count for this inner loop
        int n_cores,
        // output (pre-allocated, indexed by draws_j(ri))
        arma::cube& Wu_all, arma::cube& Yu_all, arma::cube& MY_all,
        arma::cube& beta_all, arma::cube& sigma_all) {

    arma::uword n_dj = draws_j.n_elem;
    arma::mat mu_bar = M_gamma * fit.mu_star;   // u×q — same for all draws in j

#ifdef _OPENMP
    if (n_cores > 1) omp_set_max_active_levels(1);
#pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
    for (std::int64_t ri = 0; ri < static_cast<std::int64_t>(n_dj); ++ri) {
        int r = static_cast<int>(draws_j(static_cast<arma::uword>(ri)));

        // Sigma ~ IW(Psi_star, nu_star)  [q×q]
        arma::mat sigma = spBPS::random::inverse_wishart(fit.nu_star, fit.Psi_star);
        arma::mat L_sig = safe_chol(sigma, "lower");

        // (beta, w) ~ MN(mu_star, V_star, sigma)
        arma::mat b    = spBPS::random::matrix_normal(
                             fit.mu_star, L_Vstar_j, L_sig, true, true);
        arma::mat beta = b.rows(0, p-1);
        arma::mat w    = b.rows(p, p + n_subset - 1);

        // W_u ~ MN(M_u*w, V_u, sigma)  — L_Vu reused, no extra Cholesky
        arma::mat Wu = spBPS::random::matrix_normal(
                           M_u * w, L_Vu, L_sig, true, true);

        // Y_u = X_u*beta + W_u + sqrt(a)*Z*L_sig'  [diagonal noise, O(u*q²)]
        arma::mat Yu = X_u * beta + Wu
                       + std::sqrt(a_j) * spBPS::random::normal_matrix(Wu.n_rows, sigma.n_rows)
                         * L_sig.t();

        arma::uword ur = static_cast<arma::uword>(r);
        Wu_all.slice(ur)    = Wu;
        Yu_all.slice(ur)    = Yu;
        MY_all.slice(ur)    = mu_bar;
        beta_all.slice(ur)  = beta;
        sigma_all.slice(ur) = sigma;
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// EXPORTED WRAPPER FUNCTIONS
// ═══════════════════════════════════════════════════════════════════════════


// ═══════════════════════════════════════════════════════════════════════════
// BPS_post_all_subsets_v2
// OMP-parallel posterior-only sampling over all K subsets.
// Replaces the R-level loop in spBPS_new() for the draws>0, newdata=NULL case.
// Running OMP ensures BLAS is called from worker threads, not the main R
// thread — avoids the global BLAS mutex that causes 500x slowdown at 1 core.
// ═══════════════════════════════════════════════════════════════════════════
//' @keywords internal
// [[Rcpp::export]]
List BPS_post_all_subsets_v2(
        const List& Y_list, const List& X_list, const List& crd_list,
        const List& W_local_list,
        const arma::uvec& subset_ind,   // R-length, 0-indexed: draw r → subset k
        const List& priors,
        const List& hyperpar,
        int n_cores) {

    int K = Y_list.size();
    int R = (int)subset_ind.n_elem;

    std::vector<arma::mat> Y_v(K), X_v(K), crd_v(K);
    std::vector<arma::vec> W_v(K);
    for (int k = 0; k < K; ++k) {
        Y_v[k]   = as<arma::mat>(Y_list[k]);
        X_v[k]   = as<arma::mat>(X_list[k]);
        crd_v[k] = as<arma::mat>(crd_list[k]);
        arma::vec wk = as<arma::vec>(W_local_list[k]);
        double sw = arma::accu(wk); if (sw > 0.0) wk /= sw;
        W_v[k] = wk;
    }

    arma::mat mu_B = as<arma::mat>(priors["mu_B"]);
    arma::mat V_r  = as<arma::mat>(priors["V_r"]);
    arma::mat Psi  = as<arma::mat>(priors["Psi"]);
    double nu      = as<double>(priors["nu"]);
    arma::mat Grid = expand_grid_cpp(as<arma::vec>(hyperpar["alpha"]),
                                      as<arma::vec>(hyperpar["phi"]));
    int J = Grid.n_rows;
    int p = X_v[0].n_cols;
    int q = Y_v[0].n_cols;

    arma::mat iV_r;
    if (!arma::inv_sympd(iV_r, V_r)) { V_r.diag() += 1e-10; arma::inv_sympd(iV_r, V_r); }
    spBPS::random::seed_from_R();

    // Build subset→draws mapping
    std::vector<std::vector<int>> subset_draws(K);
    for (int r = 0; r < R; ++r)
        subset_draws[(int)subset_ind(r)].push_back(r);

    // Output arrays (pre-allocated, written to non-overlapping draw slices)
    arma::cube beta_all(p, q, R), sigma_all(q, q, R);
    arma::ivec model_idx(R, arma::fill::zeros);

    spBPS_par_for(K, n_cores, [&](int k) {
        if (subset_draws[k].empty()) return;

        int R_k = (int)subset_draws[k].size();
        arma::uvec kmod = spBPS::random::sample_indices(J, R_k, W_v[k], true);

        // Local fit cache for this subset (thread-local, no shared state)
        std::vector<bool>         done(J, false);
        std::vector<FitResult_v2> cache(J);
        std::vector<arma::mat>    LVs(J);

        for (int ri = 0; ri < R_k; ++ri) {
            int r = subset_draws[k][ri];
            int j = (int)kmod(ri);
            model_idx((arma::uword)r) = j;

            if (!done[j]) {
                cache[j] = fit_pure_v2(Y_v[k], X_v[k], crd_v[k], mu_B, iV_r, Psi,
                                        nu, Grid(j,0), Grid(j,1));
                LVs[j]   = safe_chol(cache[j].V_star, "lower");
                done[j]  = true;
            }

            arma::mat sigma = spBPS::random::inverse_wishart(
                                  cache[j].nu_star, cache[j].Psi_star);
            arma::mat Lsig  = safe_chol(sigma, "lower");
            arma::mat b     = spBPS::random::matrix_normal(
                                  cache[j].mu_star, LVs[j], Lsig, true, true);

            arma::uword ur = (arma::uword)r;
            beta_all.slice(ur)  = b.rows(0, p-1);
            sigma_all.slice(ur) = sigma;
        }
    }); // end spBPS_par_for

    return List::create(Named("beta") = beta_all,
                         Named("sigma") = sigma_all,
                         Named("model") = model_idx);
}

// ═══════════════════════════════════════════════════════════════════════════
// BPS_pred_streaming_v2
//
// Two-phase streaming prediction — NEVER allocates u×u matrices.
//
// Phase 1 (OMP over K subsets):
//   - Pre-sample (sigma, beta, w) for every draw
//   - Fit each needed model j once per subset k
//   - Thread-safe: each thread writes to non-overlapping draw indices
//
// Phase 2 (sequential u-batches B×B; OMP over K within each batch):
//   - For each batch b: allocate only Bb×q×R (≤1.6 MB for B=500, R=200)
//   - chol(V_u_b) computed ONCE per (k,j,batch), shared across all draws
//   - Use pre-sampled (sigma, beta, w) — no re-sampling
//   - Writes Wu_b, Yu_b for the Bb sites of this batch
//
// return_summary=FALSE: return Wu[u×q×R], Yu[u×q×R], MY[u×q×R]
// return_summary=TRUE:  return mu_Yu, sd_Yu, q_lo, q_hi, mu_Wu  [u×q each]
//                       using exact quantiles (sorted R values per site/batch)
// ═══════════════════════════════════════════════════════════════════════════
//' @keywords internal
// [[Rcpp::export]]
List BPS_pred_streaming_v2(
        const List& Y_list, const List& X_list, const List& crd_list,
        const List& W_local_list,
        const arma::uvec& subset_ind,   // R-length, 0-indexed
        const arma::mat&  crd_u,        // u×d
        const arma::mat&  X_u,          // u×p
        const List& priors,
        const List& hyperpar,
        int  pred_batch_size,           // B; 0 = auto (min(500,u))
        int  n_cores,
        bool return_summary,
        double prob_lo,
        double prob_hi) {

    // ── Unpack ────────────────────────────────────────────────────────────
    int K = Y_list.size();
    int R = (int)subset_ind.n_elem;
    int u = X_u.n_rows, p = X_u.n_cols;
    int B = (pred_batch_size > 0) ? std::min(pred_batch_size, u) : std::min(500, u);
    int n_batches = (int)std::ceil((double)u / B);

    std::vector<arma::mat> Y_v(K), X_v(K), crd_v(K);
    std::vector<arma::vec> W_v(K);
    for (int k = 0; k < K; ++k) {
        Y_v[k]   = as<arma::mat>(Y_list[k]);
        X_v[k]   = as<arma::mat>(X_list[k]);
        crd_v[k] = as<arma::mat>(crd_list[k]);
        arma::vec wk = as<arma::vec>(W_local_list[k]);
        double sw = arma::accu(wk); if (sw > 0.0) wk /= sw;
        W_v[k] = wk;
    }

    arma::mat mu_B = as<arma::mat>(priors["mu_B"]);
    arma::mat V_r  = as<arma::mat>(priors["V_r"]);
    arma::mat Psi  = as<arma::mat>(priors["Psi"]);
    double nu      = as<double>(priors["nu"]);
    arma::mat Grid = expand_grid_cpp(as<arma::vec>(hyperpar["alpha"]),
                                      as<arma::vec>(hyperpar["phi"]));
    int J  = Grid.n_rows;
    int q  = Y_v[0].n_cols;

    arma::mat iV_r;
    if (!arma::inv_sympd(iV_r, V_r)) { V_r.diag() += 1e-10; arma::inv_sympd(iV_r, V_r); }
    spBPS::random::seed_from_R();

    // ── Build subset→draws mapping ──────────────────────────────────────────
    std::vector<std::vector<int>> subset_draws(K);
    for (int r = 0; r < R; ++r)
        subset_draws[(int)subset_ind(r)].push_back(r);

    // ── Per-draw pre-sampled storage ─────────────────────────────────────────
    std::vector<arma::mat> all_sigma(R), all_Lsig(R), all_beta(R), all_w(R);
    std::vector<int>       draw_j_arr(R, -1);

    // ── Per (k,j) fit cache (computed on demand in Phase 1) ────────────────
    // K×J vectors; each k processed by one thread → no race condition
    std::vector<std::vector<bool>>         fit_done(K, std::vector<bool>(J, false));
    std::vector<std::vector<FitResult_v2>> fit_cache(K, std::vector<FitResult_v2>(J));
    std::vector<std::vector<arma::mat>>    LVstar_cache(K, std::vector<arma::mat>(J));

    // ══════════════════════════════════════════════════════════════════════
    // PHASE 1: sample model indices + (sigma, beta, w) for all R draws
    // OMP over K: each thread handles one subset; writes to disjoint r indices
    // ══════════════════════════════════════════════════════════════════════
    spBPS_par_for(K, n_cores, [&](int k) {
        if (subset_draws[k].empty()) return;

        int n_k = Y_v[k].n_rows;
        int R_k = (int)subset_draws[k].size();
        arma::uvec kmod = spBPS::random::sample_indices(J, R_k, W_v[k], true);

        for (int ri = 0; ri < R_k; ++ri) {
            int r = subset_draws[k][ri];
            int j = (int)kmod(ri);
            draw_j_arr[r] = j;

            // Fit model j on subset k (once per k,j pair — thread-local)
            if (!fit_done[k][j]) {
                fit_cache[k][j]   = fit_pure_v2(Y_v[k], X_v[k], crd_v[k],
                                                  mu_B, iV_r, Psi,
                                                  nu, Grid(j,0), Grid(j,1));
                LVstar_cache[k][j] = safe_chol(fit_cache[k][j].V_star, "lower");
                fit_done[k][j]    = true;
            }

            // Pre-sample (sigma, beta, w) for this draw
            arma::mat sigma = spBPS::random::inverse_wishart(
                                  fit_cache[k][j].nu_star,
                                  fit_cache[k][j].Psi_star);
            arma::mat Lsig  = safe_chol(sigma, "lower");
            arma::mat b     = spBPS::random::matrix_normal(
                                  fit_cache[k][j].mu_star,
                                  LVstar_cache[k][j], Lsig, true, true);

            all_sigma[r] = sigma;
            all_Lsig[r]  = Lsig;
            all_beta[r]  = b.rows(0, p-1);
            all_w[r]     = b.rows(p, p + n_k - 1);  // latent field at training sites
        }
    }); // end Phase 1 spBPS_par_for

    // ══════════════════════════════════════════════════════════════════════
    // PHASE 2: u-batch prediction
    // Outer loop: sequential over u-batches (B sites at a time)
    // Inner loop: OMP over K active subsets
    // chol(V_u_b) computed ONCE per (k,j,batch) — shared across all draws in (k,j)
    // ══════════════════════════════════════════════════════════════════════
    arma::cube Wu_all, Yu_all, MY_all;
    arma::mat  mu_Yu, sd_Yu, q_lo_out, q_hi_out, mu_Wu;

    if (!return_summary) {
        Wu_all.set_size(u, q, R); Wu_all.zeros();
        Yu_all.set_size(u, q, R); Yu_all.zeros();
        MY_all.set_size(u, q, R); MY_all.zeros();
    } else {
        mu_Yu.set_size(u, q);  mu_Yu.zeros();
        sd_Yu.set_size(u, q);  sd_Yu.zeros();
        q_lo_out.set_size(u, q);
        q_hi_out.set_size(u, q);
        mu_Wu.set_size(u, q);  mu_Wu.zeros();
    }

    for (int b_u = 0; b_u < n_batches; ++b_u) {
        int b_start = b_u * B;
        int b_end   = std::min(b_start + B, u);
        int Bb      = b_end - b_start;
        arma::span rows(b_start, b_end - 1);

        arma::mat crd_u_b = crd_u.rows(rows);       // Bb×d
        arma::mat X_u_b   = X_u.rows(rows);          // Bb×p
        arma::mat d_u_bb  = arma_dist(crd_u_b);      // Bb×Bb  (NO u×u!)

        // Per-batch result cubes — small: Bb×q×R ≤ 1.6 MB
        arma::cube Wu_b_all(Bb, q, R, arma::fill::zeros);
        arma::cube Yu_b_all(Bb, q, R, arma::fill::zeros);
        arma::cube MY_b_all(Bb, q, R, arma::fill::zeros);

        // OMP over K active subsets
        spBPS_par_for(K, n_cores, [&](int k) {
            if (subset_draws[k].empty()) return;

            // Cross-distances for this batch (Bb×n_k) — recomputed each batch
            arma::mat d_cross_b = cross_dist_mat(crd_u_b, crd_v[k]);

            for (int j = 0; j < J; ++j) {
                if (!fit_done[k][j]) continue;

                // Collect draws in this (k,j) pair
                std::vector<int> dkj;
                for (int r : subset_draws[k])
                    if (draw_j_arr[r] == j) dkj.push_back(r);
                if (dkj.empty()) continue;

                double phi_j   = Grid(j, 1);
                double alpha_j = Grid(j, 0);
                double a_j     = 1.0 / alpha_j - 1.0;
                const FitResult_v2& fit = fit_cache[k][j];

                // Spatial terms — ONCE per (k,j,batch) ─────────────────────
                arma::mat Rphi_us_b = arma::exp(-phi_j * d_cross_b);     // Bb×n_k
                arma::mat M_u_b     = Rphi_us_b * fit.iRphi_s;           // Bb×n_k
                arma::mat Rphi_u_bb = arma::exp(-phi_j * d_u_bb);        // Bb×Bb
                arma::mat V_u_b     = Rphi_u_bb - M_u_b * Rphi_us_b.t();
                V_u_b = 0.5 * (V_u_b + V_u_b.t());                       // symmetrize
                arma::mat L_Vu_b    = safe_chol(V_u_b, "lower");         // ONCE per (k,j,batch)!

                // Posterior mean for this batch (same for all draws in (k,j))
                arma::mat M_gam_b = arma::join_horiz(X_u_b, M_u_b);     // Bb×(p+n_k)
                arma::mat my_b    = M_gam_b * fit.mu_star;               // Bb×q

                for (int r : dkj) {
                    arma::uword ur = (arma::uword)r;
                    // Use PRE-SAMPLED (sigma,beta,w) — no re-sampling per batch
                    arma::mat Wu_b = spBPS::random::matrix_normal(
                                         M_u_b * all_w[r], L_Vu_b,
                                         all_Lsig[r], true, true);
                    arma::mat Yu_b = X_u_b * all_beta[r] + Wu_b
                        + std::sqrt(a_j)
                          * spBPS::random::normal_matrix(Bb, q)
                          * all_Lsig[r].t();

                    Wu_b_all.slice(ur) = Wu_b;
                    Yu_b_all.slice(ur) = Yu_b;
                    MY_b_all.slice(ur) = my_b;   // same for all r in (k,j)
                }
            }
        }); // end K-loop

        // ── Accumulate batch results ────────────────────────────────────────
        if (!return_summary) {
            // Copy into global arrays (sequential, safe)
            for (int r = 0; r < R; ++r) {
                arma::uword ur = (arma::uword)r;
                Wu_all.slice(ur).rows(rows) = Wu_b_all.slice(ur);
                Yu_all.slice(ur).rows(rows) = Yu_b_all.slice(ur);
                MY_all.slice(ur).rows(rows) = MY_b_all.slice(ur);
            }
        } else {
            // Compute exact statistics per site from all R draws
            for (int i = 0; i < Bb; ++i) {
                for (int qq = 0; qq < q; ++qq) {
                    arma::vec y_vals(R), w_vals(R);
                    for (int r = 0; r < R; ++r) {
                        y_vals(r) = Yu_b_all(i, qq, r);
                        w_vals(r) = Wu_b_all(i, qq, r);
                    }
                    int si = b_start + i;
                    mu_Yu(si, qq) = arma::mean(y_vals);
                    sd_Yu(si, qq) = arma::stddev(y_vals);
                    mu_Wu(si, qq) = arma::mean(w_vals);
                    // Exact quantiles via sort
                    y_vals = arma::sort(y_vals);
                    int idx_lo = std::max(0, (int)std::round(prob_lo * (R-1)));
                    int idx_hi = std::min(R-1, (int)std::round(prob_hi * (R-1)));
                    q_lo_out(si, qq) = y_vals(idx_lo);
                    q_hi_out(si, qq) = y_vals(idx_hi);
                }
            }
        }
    } // end u-batch loop

    // ── Return ────────────────────────────────────────────────────────────
    if (!return_summary) {
        return List::create(Named("Wu") = Wu_all,
                             Named("Yu") = Yu_all,
                             Named("MY") = MY_all);
    } else {
        return List::create(Named("mu_Yu") = mu_Yu,
                             Named("sd_Yu") = sd_Yu,
                             Named("q_lo")  = q_lo_out,
                             Named("q_hi")  = q_hi_out,
                             Named("mu_Wu") = mu_Wu);
    }
}

//' @keywords internal
// [[Rcpp::export]]
List fit_cpp_MvT_v2(const List& data, const List& priors,
                     const arma::mat& coords, const List& hyperpar,
                     const arma::mat& iVr_pre) {
    FitResult_v2 r = fit_pure_v2(
        as<arma::mat>(data["Y"]), as<arma::mat>(data["X"]), coords,
        as<arma::mat>(priors["mu_B"]), iVr_pre,
        as<arma::mat>(priors["Psi"]), as<double>(priors["nu"]),
        as<double>(hyperpar["alpha"]), as<double>(hyperpar["phi"]));
    return List::create(Named("V_star")=r.V_star, Named("mu_star")=r.mu_star,
                         Named("Psi_star")=r.Psi_star, Named("nu_star")=r.nu_star,
                         Named("iRphi_s")=r.iRphi_s);
}

//' @keywords internal
// [[Rcpp::export]]
List post_draws_MvT_v2(const List& poster, const int R, const bool par, const int p) {
    arma::mat mu_star  = as<arma::mat>(poster["mu_star"]);
    arma::mat V_star   = as<arma::mat>(poster["V_star"]);
    arma::mat Psi_star = as<arma::mat>(poster["Psi_star"]);
    double    nu_star  = as<double>(poster["nu_star"]);
    arma::mat mu_full  = mu_star;
    if (par) {
        arma::uvec ip = arma::linspace<arma::uvec>(0, p-1, p);
        mu_star = mu_star.rows(0, p-1);
        V_star  = V_star.submat(ip, ip);
    }
    arma::mat L_V = safe_chol(V_star, "lower");
    List out(R);
    for (int r = 0; r < R; ++r) {
        arma::mat sigma = spBPS::random::inverse_wishart(nu_star, Psi_star);
        arma::mat beta  = spBPS::random::matrix_normal(mu_star, L_V, sigma, true, false);
        out(r) = List::create(Named("beta")=beta, Named("sigma")=sigma,
                               Named("mu_star")=mu_full);
    }
    return out;
}

//' @keywords internal
// [[Rcpp::export]]
double d_pred_cpp_MvT_v2(const List& data, const arma::mat& X_u,
                          const arma::mat& Y_u, const arma::mat& cross_d,
                          const List& hyperpar, const List& poster,
                          const arma::mat& Lc_Psi) {
    FitResult_v2 fit;
    fit.iRphi_s = as<arma::mat>(poster["iRphi_s"]);
    fit.V_star  = as<arma::mat>(poster["V_star"]);
    fit.mu_star = as<arma::mat>(poster["mu_star"]);
    fit.Psi_star= as<arma::mat>(poster["Psi_star"]);
    fit.nu_star = as<double>(poster["nu_star"]);
    return d_pred_pure_v2(Y_u.row(0), X_u.row(0), cross_d.row(0), fit,
                           as<double>(hyperpar["alpha"]), as<double>(hyperpar["phi"]));
}

//' @keywords internal
// [[Rcpp::export]]
arma::vec dens_kcv_MvT_v2(const List& data, const List& priors,
                            const arma::mat& coords, const List& hyperpar, int K) {
    arma::mat iV_r; arma::inv_sympd(iV_r, as<arma::mat>(priors["V_r"]));
    spBPS::random::seed_from_R();
    return dens_kcv_pure_v2(
        as<arma::mat>(data["Y"]), as<arma::mat>(data["X"]), coords,
        as<arma::mat>(priors["mu_B"]), iV_r,
        as<arma::mat>(priors["Psi"]), as<double>(priors["nu"]),
        as<double>(hyperpar["alpha"]), as<double>(hyperpar["phi"]), K);
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat models_dens_MvT_v2(const List& data, const List& priors,
                               const arma::mat& coords, const List& hyperpar, int K) {
    arma::mat iV_r; arma::inv_sympd(iV_r, as<arma::mat>(priors["V_r"]));
    arma::mat Grid = expand_grid_cpp(as<arma::vec>(hyperpar["alpha"]),
                                      as<arma::vec>(hyperpar["phi"]));
    spBPS::random::seed_from_R();
    arma::mat out;
    for (arma::uword j = 0; j < Grid.n_rows; ++j)
        out = arma::join_horiz(out,
            dens_kcv_pure_v2(
                as<arma::mat>(data["Y"]), as<arma::mat>(data["X"]), coords,
                as<arma::mat>(priors["mu_B"]), iV_r,
                as<arma::mat>(priors["Psi"]), as<double>(priors["nu"]),
                Grid(j,0), Grid(j,1), K));
    return out;
}

//' @keywords internal
// [[Rcpp::export]]
List BPS_weights_MvT_v2(const List& data, const List& priors,
                         const arma::mat& coords, const List& hyperpar, int K) {
    arma::mat epd = models_dens_MvT_v2(data, priors, coords, hyperpar, K);
    spBPS::stacking::SolverConfig cfg;
    cfg.method = spBPS::stacking::SolverMethod::MirrorDescent;
    cfg.tolerance = 1e-10; cfg.max_iterations = 2000;
    auto res = spBPS::stacking::solve_weights(epd, cfg);
    arma::mat W = arma::conv_to<arma::mat>::from(res.weights);
    W.elem(arma::find(W < 0.0)).zeros();
    double s = arma::accu(W); if (s > 0.0) W /= s;
    arma::mat Grid = expand_grid_cpp(as<arma::vec>(hyperpar["alpha"]),
                                      as<arma::vec>(hyperpar["phi"]));
    return List::create(Named("Grid")=arma::join_horiz(W,Grid),
                         Named("W")=W, Named("epd")=epd);
}

//' @keywords internal
// [[Rcpp::export]]
List BPS_combine_v2(const List& fit_list, const int K, const double rp) {
    List W_list(K);
    arma::mat out3 = as<arma::mat>(as<List>(fit_list[0])[0]);
    W_list[0] = as<arma::mat>(as<List>(fit_list[0])[1]);
    for (int i = 1; i < K; ++i) {
        out3 = arma::join_vert(out3, as<arma::mat>(as<List>(fit_list[i])[0]));
        W_list[i] = as<arma::mat>(as<List>(fit_list[i])[1]);
    }
    int n = out3.n_rows, N = (int)std::floor(rp * n);
    arma::uvec idx = sample_index(n, N, arma::ones<arma::vec>(n)/n);
    arma::mat out4(idx.n_elem, K, arma::fill::zeros);
    for (int k = 0; k < K; ++k)
        out4.col(k) = out3.rows(idx) * as<arma::mat>(W_list[k]);
    spBPS::stacking::SolverConfig cfg;
    cfg.method = spBPS::stacking::SolverMethod::MirrorDescent;
    cfg.tolerance = 1e-10; cfg.max_iterations = 2000;
    auto res = spBPS::stacking::solve_weights(out4, cfg);
    arma::mat Wbps = arma::conv_to<arma::mat>::from(res.weights);
    double thr = 1.0 / (2.0 * Wbps.n_elem);
    Wbps.elem(arma::find(Wbps < thr)).zeros();
    double s = arma::accu(Wbps); if (s > 0.0) Wbps /= s;
    return List::create(Named("W")=Wbps, Named("W_list")=W_list);
}

//' @keywords internal
// [[Rcpp::export]]
List compute_all_local_weights_v2(const List& Y_list, const List& X_list,
                                   const List& crd_list, const List& priors,
                                   const List& hyperpar, int cv_folds, int n_cores) {
    int K = Y_list.size();
    std::vector<arma::mat> Y_v(K), X_v(K), crd_v(K);
    for (int k = 0; k < K; ++k) {
        Y_v[k]=as<arma::mat>(Y_list[k]);
        X_v[k]=as<arma::mat>(X_list[k]);
        crd_v[k]=as<arma::mat>(crd_list[k]);
    }
    arma::mat mu_B=as<arma::mat>(priors["mu_B"]);
    arma::mat Psi =as<arma::mat>(priors["Psi"]);
    double nu=as<double>(priors["nu"]);
    arma::mat V_r_loc = as<arma::mat>(priors["V_r"]);
    arma::mat iV_r;
    if (!arma::inv_sympd(iV_r,V_r_loc)) { V_r_loc.diag()+=1e-10; arma::inv_sympd(iV_r,V_r_loc); }
    arma::mat Grid = expand_grid_cpp(as<arma::vec>(hyperpar["alpha"]),
                                      as<arma::vec>(hyperpar["phi"]));
    int J = Grid.n_rows;
    spBPS::stacking::SolverConfig cfg;
    cfg.method=spBPS::stacking::SolverMethod::MirrorDescent;
    cfg.tolerance=1e-10; cfg.max_iterations=2000;
    std::vector<arma::mat> epd_v(K), W_v(K);

    // if(n_thr > 1) prevents spinlock overhead when n_cores == 1.
    // omp_set_max_active_levels(1) is ONLY called when actually using OMP;
    // calling it unconditionally restricts BLAS threads even in serial mode.
#ifdef _OPENMP
    int n_thr = std::max(1, n_cores);
    if (n_thr > 1) omp_set_max_active_levels(1);
#pragma omp parallel for schedule(dynamic) num_threads(n_thr) if(n_thr > 1)
#else
    (void)n_cores;
#endif
    for (int k = 0; k < K; ++k) {
        arma::mat epd_k(Y_v[k].n_rows, J);
        for (int j = 0; j < J; ++j)
            epd_k.col(j) = dens_kcv_pure_v2(Y_v[k], X_v[k], crd_v[k],
                                              mu_B, iV_r, Psi, nu,
                                              Grid(j,0), Grid(j,1), cv_folds);
        // Force single-threaded OMP inside solve_weights.
        // compute_density/compute_gradient have #pragma omp parallel
        // without num_threads — on Windows each team creation costs
        // ~0.025ms, totalling ~22s for 2000 iters × K × J × folds.
#ifdef _OPENMP
        int _prev_threads = omp_get_max_threads();
        omp_set_num_threads(1);
#endif
        auto res = spBPS::stacking::solve_weights(epd_k, cfg);
#ifdef _OPENMP
        omp_set_num_threads(_prev_threads);
#endif
        arma::mat Wk = arma::conv_to<arma::mat>::from(res.weights);
        Wk.elem(arma::find(Wk < 0.0)).zeros();
        double s = arma::accu(Wk); if (s > 0.0) Wk /= s;
        epd_v[k]=epd_k; W_v[k]=Wk;
    }
    List fit_list(K);
    for (int k = 0; k < K; ++k)
        fit_list[k] = List::create(epd_v[k], W_v[k]);
    return fit_list;
}

//' @keywords internal
// [[Rcpp::export]]
List BPS_post_MvT_v2(const List& data, const arma::mat& X_u,
                      const List& priors, const arma::mat& coords,
                      const arma::mat& crd_u, const List& hyperpar,
                      const arma::vec& W, const int R,
                      const arma::mat& d_u,
                      const List& Rphi_u_list,   // J pre-computed exp(-phi_j*d_u)
                      int n_cores) {
    arma::mat Y   = as<arma::mat>(data["Y"]);
    arma::mat X   = as<arma::mat>(data["X"]);
    arma::mat mu_B= as<arma::mat>(priors["mu_B"]);
    arma::mat V_r = as<arma::mat>(priors["V_r"]);
    arma::mat Psi = as<arma::mat>(priors["Psi"]);
    double nu     = as<double>(priors["nu"]);
    arma::mat Grid= expand_grid_cpp(as<arma::vec>(hyperpar["alpha"]),
                                     as<arma::vec>(hyperpar["phi"]));
    int J  = Grid.n_rows;
    int u  = X_u.n_rows, p = X_u.n_cols;
    int n  = Y.n_rows,   q = Y.n_cols;

    arma::mat iV_r;
    if (!arma::inv_sympd(iV_r,V_r)){V_r.diag()+=1e-10;arma::inv_sympd(iV_r,V_r);}

    spBPS::random::seed_from_R();
    arma::vec W_norm = W / arma::accu(W);
    arma::uvec kmod  = spBPS::random::sample_indices(J, R, W_norm, true);

    // Cross-distance u×n (computed once per call)
    arma::mat d_cross = cross_dist_mat(crd_u, coords);

    // Output arrays — contiguous memory, no list overhead
    arma::cube Wu_all(u,q,R), Yu_all(u,q,R), MY_all(u,q,R);
    arma::cube beta_all(p,q,R), sigma_all(q,q,R);
    arma::ivec model_idx(R);

    for (int j = 0; j < J; ++j) {
        // Collect draw indices for model j
        arma::uvec draws_j;
        for (int r = 0; r < R; ++r) {
            if ((int)kmod(r) == j) {
                arma::uvec tmp(1); tmp(0) = (arma::uword)r;
                draws_j = arma::join_vert(draws_j, tmp);
            }
        }
        if (draws_j.is_empty()) continue;

        model_idx.elem(draws_j) = arma::ivec(draws_j.n_elem,
                                              arma::fill::value(j));

        double alpha_j = Grid(j,0), phi_j = Grid(j,1);
        double a_j     = 1.0/alpha_j - 1.0;

        // Fit model j on this subset's data
        FitResult_v2 fit = fit_pure_v2(Y, X, coords, mu_B, iV_r, Psi,
                                        nu, alpha_j, phi_j);
        arma::mat L_Vstar_j = safe_chol(fit.V_star, "lower");

        // Spatial terms — Rphi_u from pre-computed list
        arma::mat Rphi_us = arma::exp(-phi_j * d_cross);          // u×n
        arma::mat M_u     = Rphi_us * fit.iRphi_s;                // u×n
        arma::mat Rphi_u  = as<arma::mat>(Rphi_u_list[j]);        // u×u — pre-computed!
        arma::mat V_u     = 0.5*(Rphi_u - M_u*Rphi_us.t());
        V_u               = V_u + V_u.t();                        // ensure symmetry
        arma::mat L_Vu    = safe_chol(V_u, "lower");

        arma::mat M_gamma = arma::join_horiz(X_u, M_u);           // u×(p+n)

        // Draw loop — OMP-parallelised over R_j draws for this model j
        draw_group_j(fit, alpha_j, a_j, M_u, L_Vu, L_Vstar_j, M_gamma,
                     X_u, draws_j, p, n, n_cores,
                     Wu_all, Yu_all, MY_all, beta_all, sigma_all);
    }

    // Return compact arrays — no list-of-R-lists overhead
    return List::create(
        Named("predictive") = List::create(
            Named("Wu") = Wu_all,
            Named("Yu") = Yu_all,
            Named("MY") = MY_all),
        Named("posterior") = List::create(
            Named("beta")  = beta_all,
            Named("sigma") = sigma_all,
            Named("model") = model_idx));
}

// ═══════════════════════════════════════════════════════════════════════════
// PRE-COMPUTED BATCH STRUCTURE
// Holds chol(V_u_b) for each batch — computed once per model j, shared across
// all R_j draws. This is the key fix: Cholesky NOT inside the draw loop.
// ═══════════════════════════════════════════════════════════════════════════
struct BatchCache {
    std::vector<arma::mat> L_Vu_b;   // chol(V_u_b) per batch
    std::vector<arma::mat> M_u_b;    // kriging weights per batch (B×n)
    std::vector<arma::mat> X_u_b;    // covariates per batch (B×p)
    std::vector<arma::span> spans;   // row spans for each batch
    int n_batches;
};

static BatchCache build_batch_cache(
        const arma::mat& M_u,       // u×n
        const arma::mat& Rphi_us,   // u×n
        const arma::mat& Rphi_u,    // u×u (pre-computed)
        const arma::mat& X_u,       // u×p
        int u, int B) {

    BatchCache bc;
    bc.n_batches = (int)std::ceil((double)u / B);
    bc.L_Vu_b.resize(bc.n_batches);
    bc.M_u_b.resize(bc.n_batches);
    bc.X_u_b.resize(bc.n_batches);
    bc.spans.reserve(bc.n_batches);

    for (int b = 0; b < bc.n_batches; ++b) {
        int b_start = b * B;
        int b_end   = std::min(b_start + B, u) - 1;
        arma::span rows(b_start, b_end);

        bc.spans.push_back(rows);
        arma::mat M_b   = M_u.rows(rows);          // B×n
        arma::mat Rph_b = Rphi_us.rows(rows);       // B×n
        arma::mat Ru_b  = Rphi_u.submat(rows,rows); // B×B
        arma::mat V_b   = Ru_b - M_b * Rph_b.t();  // B×B
        V_b = 0.5 * (V_b + V_b.t());               // symmetrize

        bc.M_u_b[b] = M_b;
        bc.X_u_b[b] = X_u.rows(rows);
        bc.L_Vu_b[b] = safe_chol(V_b, "lower");    // ← paid ONCE per j per batch
    }
    return bc;
}

//' @keywords internal
// [[Rcpp::export]]
List BPS_post_MvT_v2_batch(const List& data, const arma::mat& X_u,
                             const List& priors, const arma::mat& coords,
                             const arma::mat& crd_u, const List& hyperpar,
                             const arma::vec& W, const int R,
                             const arma::mat& d_u,
                             const List& Rphi_u_list,
                             int n_cores,
                             int pred_batch_size) {
    arma::mat Y   = as<arma::mat>(data["Y"]);
    arma::mat X   = as<arma::mat>(data["X"]);
    arma::mat mu_B= as<arma::mat>(priors["mu_B"]);
    arma::mat V_r = as<arma::mat>(priors["V_r"]);
    arma::mat Psi = as<arma::mat>(priors["Psi"]);
    double nu     = as<double>(priors["nu"]);
    arma::mat Grid= expand_grid_cpp(as<arma::vec>(hyperpar["alpha"]),
                                     as<arma::vec>(hyperpar["phi"]));
    int J = Grid.n_rows;
    int u = X_u.n_rows, p = X_u.n_cols, n = Y.n_rows, q = Y.n_cols;
    int B = std::min(pred_batch_size, u);

    arma::mat iV_r;
    if (!arma::inv_sympd(iV_r,V_r)){V_r.diag()+=1e-10;arma::inv_sympd(iV_r,V_r);}

    spBPS::random::seed_from_R();
    arma::vec W_norm = W / arma::accu(W);
    arma::uvec kmod  = spBPS::random::sample_indices(J, R, W_norm, true);

    arma::mat d_cross = cross_dist_mat(crd_u, coords);   // u×n

    arma::cube Wu_all(u,q,R,arma::fill::zeros),
               Yu_all(u,q,R,arma::fill::zeros),
               MY_all(u,q,R,arma::fill::zeros);
    arma::cube beta_all(p,q,R), sigma_all(q,q,R);
    arma::ivec model_idx(R);

    for (int j = 0; j < J; ++j) {
        arma::uvec draws_j;
        for (int r = 0; r < R; ++r) {
            if ((int)kmod(r) == j) {
                arma::uvec tmp(1); tmp(0)=(arma::uword)r;
                draws_j = arma::join_vert(draws_j, tmp);
            }
        }
        if (draws_j.is_empty()) continue;
        model_idx.elem(draws_j) = arma::ivec(draws_j.n_elem,
                                              arma::fill::value(j));

        double alpha_j=Grid(j,0), phi_j=Grid(j,1), a_j=1.0/alpha_j-1.0;
        FitResult_v2 fit = fit_pure_v2(Y,X,coords,mu_B,iV_r,Psi,nu,alpha_j,phi_j);
        arma::mat L_Vstar_j = safe_chol(fit.V_star,"lower");

        arma::mat Rphi_us = arma::exp(-phi_j * d_cross);          // u×n
        arma::mat M_u     = Rphi_us * fit.iRphi_s;                // u×n
        arma::mat Rphi_u  = as<arma::mat>(Rphi_u_list[j]);        // u×u pre-computed

        // ── Precompute ALL batch Cholesky factors ONCE per model j ──────────
        // This is the key fix: O(B³ × n_batches) paid once, NOT per draw.
        BatchCache bc = build_batch_cache(M_u, Rphi_us, Rphi_u, X_u, u, B);

        // Posterior mean (full u×q — same for all draws)
        arma::mat M_gamma = arma::join_horiz(X_u, M_u);
        arma::mat mu_bar  = M_gamma * fit.mu_star;

        arma::uword n_dj = draws_j.n_elem;

#ifdef _OPENMP
        if (n_cores > 1) omp_set_max_active_levels(1);
#pragma omp parallel for schedule(static) num_threads(n_cores) if(n_cores > 1)
#endif
        for (std::int64_t ri = 0; ri < (std::int64_t)n_dj; ++ri) {
            int r = (int)draws_j((arma::uword)ri);
            arma::uword ur = (arma::uword)r;

            // Sample sigma and (beta, w) — O(q³) + O((p+n)² q)
            arma::mat sigma = spBPS::random::inverse_wishart(fit.nu_star, fit.Psi_star);
            arma::mat L_sig = safe_chol(sigma,"lower");
            arma::mat b     = spBPS::random::matrix_normal(
                                  fit.mu_star, L_Vstar_j, L_sig, true, true);
            arma::mat beta  = b.rows(0, p-1);
            arma::mat w     = b.rows(p, p+n-1);

            beta_all.slice(ur)  = beta;
            sigma_all.slice(ur) = sigma;
            MY_all.slice(ur)    = mu_bar;

            // ── Sample W_u and Y_u batch by batch ─────────────────────────
            // L_Vu_b is pre-computed outside this loop → no Cholesky here.
            // Each batch is O(B² q) matrix multiply — much cheaper than O(u² q).
            arma::mat Wu_full(u,q), Yu_full(u,q);
            for (int b_idx = 0; b_idx < bc.n_batches; ++b_idx) {
                const arma::span& rows = bc.spans[b_idx];
                const arma::mat&  M_b  = bc.M_u_b[b_idx];   // B×n
                const arma::mat&  Xb   = bc.X_u_b[b_idx];   // B×p
                const arma::mat&  LVb  = bc.L_Vu_b[b_idx];  // B×B chol

                // W_u batch: MN(M_b * w, V_u_b, sigma) — L_sig + LVb reused
                arma::mat Wu_b = spBPS::random::matrix_normal(
                                     M_b * w, LVb, L_sig, true, true);
                // Y_u batch: diagonal noise
                arma::mat Yu_b = Xb * beta + Wu_b
                    + std::sqrt(a_j)
                      * spBPS::random::normal_matrix(Wu_b.n_rows, q)
                      * L_sig.t();

                Wu_full.rows(rows) = Wu_b;
                Yu_full.rows(rows) = Yu_b;
            }
            Wu_all.slice(ur) = Wu_full;
            Yu_all.slice(ur) = Yu_full;
        }
    }

    return List::create(
        Named("predictive") = List::create(
            Named("Wu")=Wu_all, Named("Yu")=Yu_all, Named("MY")=MY_all),
        Named("posterior") = List::create(
            Named("beta")=beta_all, Named("sigma")=sigma_all,
            Named("model")=model_idx));
}


// ═══════════════════════════════════════════════════════════════════════════
// FULLY PARALLEL PREDICTION — OMP over K subsets
//
// Architecture:
//   OMP parallel for → K subsets  (each subset independent)
//     for j in unique models of this subset:
//       precompute L_Vu once per j (or per j × batch)
//       for r in draws of model j:
//         sample sigma, beta, w, Wu, Yu
//
// Why this beats OMP-on-draws:
//   With K=4000, R=200: ~200 active subsets, each with 1-2 draws.
//   OMP-on-draws gives 1-2 elements per thread → useless.
//   OMP-on-subsets gives ~200/n_cores subsets per thread → efficient.
//
// All inputs (d_u, Rphi_u_list, X_u) are read-only → thread-safe.
// Output arrays written to different draw-slices → no race conditions.
// ═══════════════════════════════════════════════════════════════════════════

//' @keywords internal
// [[Rcpp::export]]
List BPS_pred_all_subsets_v2(
        const List& Y_list, const List& X_list, const List& crd_list,
        const List& W_local_list,
        const arma::uvec& subset_ind,    // R-length: draw r → subset k
        const arma::mat&  d_u,
        const List&       Rphi_u_list,   // J pre-computed u×u matrices
        const arma::mat&  X_u,
        const arma::mat&  crd_u,
        const List& priors,
        const List& hyperpar,
        int pred_batch_size,
        int n_cores) {

    int K = Y_list.size();
    int R = (int)subset_ind.n_elem;
    int u = X_u.n_rows, p = X_u.n_cols;

    // ── Extract K subsets into C++ before parallel ──────────────────────
    std::vector<arma::mat> Y_v(K), X_v(K), crd_v(K);
    std::vector<arma::vec> W_v(K);
    for (int k = 0; k < K; ++k) {
        Y_v[k]   = as<arma::mat>(Y_list[k]);
        X_v[k]   = as<arma::mat>(X_list[k]);
        crd_v[k] = as<arma::mat>(crd_list[k]);
        W_v[k]   = as<arma::vec>(W_local_list[k]);
    }

    // ── Extract shared data ──────────────────────────────────────────────
    arma::mat mu_B = as<arma::mat>(priors["mu_B"]);
    arma::mat V_r  = as<arma::mat>(priors["V_r"]);
    arma::mat Psi  = as<arma::mat>(priors["Psi"]);
    double    nu   = as<double>(priors["nu"]);
    arma::mat Grid = expand_grid_cpp(as<arma::vec>(hyperpar["alpha"]),
                                      as<arma::vec>(hyperpar["phi"]));
    int J = Grid.n_rows;
    int q = Y_v[0].n_cols;
    bool do_batch = (pred_batch_size > 0 && pred_batch_size < u);
    int  B        = do_batch ? pred_batch_size : u;

    arma::mat iV_r;
    if (!arma::inv_sympd(iV_r,V_r)){V_r.diag()+=1e-10;arma::inv_sympd(iV_r,V_r);}

    // Pre-extract Rphi_u matrices from R list (read-only, shared across threads)
    std::vector<arma::mat> Rphi_u_v(J);
    for (int j = 0; j < J; ++j)
        Rphi_u_v[j] = as<arma::mat>(Rphi_u_list[j]);

    // ── Build subset→draw mapping (which draws go to which subset) ───────
    // subset_ind[r] = k means draw r uses subset k (0-indexed)
    std::vector<std::vector<int>> subset_draws(K);
    for (int r = 0; r < R; ++r)
        subset_draws[(int)subset_ind(r)].push_back(r);

    // ── Output arrays (pre-allocated, written by different threads) ──────
    arma::cube Wu_all(u,q,R,arma::fill::zeros),
               Yu_all(u,q,R,arma::fill::zeros),
               MY_all(u,q,R,arma::fill::zeros);
    arma::cube beta_all(p,q,R), sigma_all(q,q,R);
    arma::ivec model_idx(R, arma::fill::zeros);

    spBPS::random::seed_from_R();

    spBPS_par_for(K, n_cores, [&](int k) {
        const std::vector<int>& draws_k = subset_draws[k];
        if (draws_k.empty()) return;

        int n_k = Y_v[k].n_rows;
        const arma::mat& Yk   = Y_v[k];
        const arma::mat& Xk   = X_v[k];
        const arma::mat& crdk = crd_v[k];

        // Normalise local weights
        arma::vec Wk = W_v[k];
        double sw = arma::accu(Wk); if (sw > 0.0) Wk /= sw;

        // Sample model index for each draw assigned to this subset
        arma::vec p_j(J);
        for (int j = 0; j < J; ++j) p_j(j) = Wk(j);
        arma::uvec kmod = spBPS::random::sample_indices(
            J, draws_k.size(), p_j, true);

        // Cross-distance crd_u × crdk  (u×n_k)
        arma::mat d_cross = cross_dist_mat(crd_u, crdk);

        // Process draws grouped by model j
        for (int j = 0; j < J; ++j) {
            // Which draws (within draws_k) use model j?
            std::vector<int> dj_local;
            for (int ri = 0; ri < (int)draws_k.size(); ++ri)
                if ((int)kmod(ri) == j) dj_local.push_back(ri);
            if (dj_local.empty()) continue;

            double alpha_j = Grid(j,0), phi_j = Grid(j,1);
            double a_j     = 1.0/alpha_j - 1.0;

            // Fit model j on subset k
            FitResult_v2 fit = fit_pure_v2(Yk,Xk,crdk,mu_B,iV_r,Psi,
                                            nu,alpha_j,phi_j);
            arma::mat L_Vstar_j = safe_chol(fit.V_star,"lower");

            arma::mat Rphi_us = arma::exp(-phi_j * d_cross);     // u×n_k
            arma::mat M_u     = Rphi_us * fit.iRphi_s;           // u×n_k

            // Posterior mean (u×q, same for all draws in j)
            arma::mat M_gamma = arma::join_horiz(X_u, M_u);
            arma::mat mu_bar  = M_gamma * fit.mu_star;

            if (!do_batch) {
                // ── Full u×u Cholesky — exact correlations ───────────────
                arma::mat V_u = arma::mat(Rphi_u_v[j])
                                - M_u * Rphi_us.t();
                V_u = 0.5*(V_u + V_u.t());
                arma::mat L_Vu = safe_chol(V_u,"lower");

                for (int ri : dj_local) {
                    int r = draws_k[ri];
                    arma::uword ur = (arma::uword)r;
                    model_idx(ur) = j;

                    arma::mat sigma = spBPS::random::inverse_wishart(
                                          fit.nu_star, fit.Psi_star);
                    arma::mat L_sig = safe_chol(sigma,"lower");
                    arma::mat b = spBPS::random::matrix_normal(
                                      fit.mu_star,L_Vstar_j,L_sig,true,true);
                    arma::mat beta = b.rows(0,p-1);
                    arma::mat w    = b.rows(p,p+n_k-1);

                    arma::mat Wu = spBPS::random::matrix_normal(
                                       M_u*w, L_Vu, L_sig, true, true);
                    arma::mat Yu = X_u*beta + Wu
                        + std::sqrt(a_j)
                          *spBPS::random::normal_matrix(u,q)*L_sig.t();

                    Wu_all.slice(ur)    = Wu;
                    Yu_all.slice(ur)    = Yu;
                    MY_all.slice(ur)    = mu_bar;
                    beta_all.slice(ur)  = beta;
                    sigma_all.slice(ur) = sigma;
                }
            } else {
                // ── Batched Cholesky — precompute all L_Vu_b once per j ──
                BatchCache bc = build_batch_cache(
                    M_u, Rphi_us, Rphi_u_v[j], X_u, u, B);

                for (int ri : dj_local) {
                    int r = draws_k[ri];
                    arma::uword ur = (arma::uword)r;
                    model_idx(ur) = j;

                    arma::mat sigma = spBPS::random::inverse_wishart(
                                          fit.nu_star, fit.Psi_star);
                    arma::mat L_sig = safe_chol(sigma,"lower");
                    arma::mat b = spBPS::random::matrix_normal(
                                      fit.mu_star,L_Vstar_j,L_sig,true,true);
                    arma::mat beta = b.rows(0,p-1);
                    arma::mat w    = b.rows(p,p+n_k-1);

                    arma::mat Wu_full(u,q), Yu_full(u,q);
                    for (int bi = 0; bi < bc.n_batches; ++bi) {
                        const arma::span& rows = bc.spans[bi];
                        arma::mat Wu_b = spBPS::random::matrix_normal(
                            bc.M_u_b[bi]*w, bc.L_Vu_b[bi], L_sig, true, true);
                        arma::mat Yu_b = bc.X_u_b[bi]*beta + Wu_b
                            + std::sqrt(a_j)
                              *spBPS::random::normal_matrix(Wu_b.n_rows,q)
                              *L_sig.t();
                        Wu_full.rows(rows) = Wu_b;
                        Yu_full.rows(rows) = Yu_b;
                    }
                    Wu_all.slice(ur)    = Wu_full;
                    Yu_all.slice(ur)    = Yu_full;
                    MY_all.slice(ur)    = mu_bar;
                    beta_all.slice(ur)  = beta;
                    sigma_all.slice(ur) = sigma;
                }
            }
        } // j loop
    }); // end K-loop (spBPS_par_for)

    return List::create(
        Named("predictive") = List::create(
            Named("Wu")=Wu_all, Named("Yu")=Yu_all, Named("MY")=MY_all),
        Named("posterior") = List::create(
            Named("beta")=beta_all, Named("sigma")=sigma_all,
            Named("model")=model_idx));
}

//' @keywords internal
// [[Rcpp::export]]
List BPS_postdraws_MvT_v2(const List& data, const List& priors,
                           const arma::mat& coords, const List& hyperpar,
                           const arma::vec& W, const int R, bool par) {
    arma::mat Y=as<arma::mat>(data["Y"]);
    arma::mat X=as<arma::mat>(data["X"]);
    arma::mat mu_B=as<arma::mat>(priors["mu_B"]);
    arma::mat V_r_=as<arma::mat>(priors["V_r"]);
    arma::mat Psi=as<arma::mat>(priors["Psi"]);
    double nu=as<double>(priors["nu"]);
    int p=V_r_.n_cols, q=Y.n_cols;
    arma::mat iV_r; arma::inv_sympd(iV_r,V_r_);
    arma::mat Grid=expand_grid_cpp(as<arma::vec>(hyperpar["alpha"]),
                                    as<arma::vec>(hyperpar["phi"]));
    int J=Grid.n_rows;
    arma::vec W_norm=W/arma::accu(W);
    spBPS::random::seed_from_R();

    // Return arrays for efficiency
    arma::cube beta_arr(p,q,R), sigma_arr(q,q,R);
    arma::ivec model_arr(R);

    for (int r = 0; r < R; ++r) {
        arma::uvec km = spBPS::random::sample_indices(J,1,W_norm,true);
        int j=(int)km(0);
        model_arr(r) = j;
        FitResult_v2 fit=fit_pure_v2(Y,X,coords,mu_B,iV_r,Psi,nu,Grid(j,0),Grid(j,1));
        arma::mat mu_use, V_use;
        if (par) {
            arma::uvec ip=arma::regspace<arma::uvec>(0,p-1);
            mu_use=fit.mu_star.rows(0,p-1); V_use=fit.V_star.submat(ip,ip);
        } else { mu_use=fit.mu_star; V_use=fit.V_star; }
        arma::mat sigma=spBPS::random::inverse_wishart(fit.nu_star,fit.Psi_star);
        arma::mat beta=spBPS::random::matrix_normal(mu_use,V_use,sigma);
        beta_arr.slice(r)=beta; sigma_arr.slice(r)=sigma;
    }
    return List::create(Named("beta")=beta_arr, Named("sigma")=sigma_arr,
                         Named("model")=model_arr);
}


// ═══════════════════════════════════════════════════════════════════════════
// BPS_PseudoBMA_v2 — OMP-parallel version of BPS_PseudoBMA
//
// Statistically IDENTICAL to BPS_PseudoBMA. The outer j-loop is
// parallelised: each j writes exclusively to out4.row(j), sp is local
// per iteration → no race conditions.
//
// R API (as<>, List[]) is NOT thread-safe: all extractions happen
// BEFORE the OMP region in plain sequential C++.
//
// Expected speedup: ~n_cores × (sequential cost of BPS_PseudoBMA)
// minus pre-extraction overhead O(K) — negligible vs O(K²) loop.
// ═══════════════════════════════════════════════════════════════════════════
//' @keywords internal
// [[Rcpp::export]]
List BPS_PseudoBMA_v2(const List& fit_list, int n_cores = 1) {

    int K = fit_list.size();

    // ── Pre-extract ALL matrices before OMP (R API not thread-safe) ──────
    std::vector<arma::mat> epd_v(K);   // n_k × J  log-density matrix per subset
    std::vector<arma::mat> W_v(K);     // J × 1    local weight vector per subset
    List W_list_out(K);                // for the return value (compatible with BPS_PseudoBMA output)

    for (int i = 0; i < K; ++i) {
        List ls  = as<List>(fit_list[i]);
        epd_v[i] = as<arma::mat>(ls[0]);
        W_v[i]   = as<arma::mat>(ls[1]);
        W_list_out[i] = W_v[i];
    }

    int nn  = (int)epd_v[0].n_rows;
    arma::mat out4(K, K, arma::fill::zeros);

    // ── OMP parallel for over j ───────────────────────────────────────────
    // Each j: reads epd_v[k] (shared, read-only) and W_v[j] (per-j read-only),
    // writes out4.row(j) (exclusively owned by this j) → no race condition.
    // sp declared inside the loop → private to each thread automatically.
#ifdef _OPENMP
    int n_thr = std::max(1, n_cores);
    if (n_thr > 1) omp_set_max_active_levels(1);
#pragma omp parallel for schedule(static) num_threads(n_thr) if(n_thr > 1)
#else
    (void)n_cores;
#endif
    for (int j = 0; j < K; ++j) {
        arma::mat sp(nn, K, arma::fill::zeros);
        const arma::mat& W_j = W_v[j];       // read-only reference, no copy
        for (int k = 0; k < K; ++k) {
            sp.col(k) = arma::log(epd_v[k] * W_j);
        }
        sp.elem(arma::find_nonfinite(sp)).zeros();
        out4.row(j) = arma::sum(sp, 0);
    }

    // ── Compute elpd and pseudo-BMA weights (sequential, O(K)) ──────────
    out4 = out4 / nn;
    arma::mat out8 = arma::trans(arma::mean(out4, 0));
    arma::mat Waic = arma::exp(out8);
    Waic = Waic / arma::sum(Waic, 0).eval()(0, 0);
    Waic.elem(arma::find(Waic <= 1.0 / (2.0 * (double)Waic.size()))).zeros();
    Waic /= arma::sum(Waic, 0).eval()(0, 0);

    return List::create(Named("W")      = Waic,
                         Named("W_list") = W_list_out);
}

//' @keywords internal
// [[Rcpp::export]]
double dmatrix_t_test(const arma::mat& X, const arma::mat& Lambda,
                       const arma::mat& SigmaR, const arma::mat& SigmaC,
                       double nu, bool log_d = false) {
    return spBPS::random::dmatrix_t(X, Lambda, SigmaR, SigmaC, nu, log_d);
}
