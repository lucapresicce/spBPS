# packages
# library(spBPS)
# devtools::install(upgrade = "never", build_vignettes = FALSE, dependencies = FALSE)
library(mniw)
library(tictoc)
library(parallel)
library(doParallel)
library(foreach)

# dimensions
set.seed(2025)
n <- 1000
u <- 100
p <- 2
q <- 2

# parameters
B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
sigma2 <- matrix(c(1, 0.5, 0.5, 1), q, q)
# B <- matrix(c(-0.75, -1.1), p, q)
# sigma2 <- 1
alfa <- 0.8
phi <- 4

set.seed(97)
# generate sintethic data
crd <- matrix(runif((n+u) * 2), ncol = 2)
X_or <- cbind(rep(1, n+u), matrix(runif((p-1)*(n+u)), ncol = (p-1)))
D <- spBPS:::arma_dist(crd)
gc()
Rphi <- exp(-phi * D)
rm("D"); gc()
W_or <- matrix(0, n+u, q) + mniw::rMNorm(1, Lambda = matrix(0, n+u, q), SigmaR = Rphi, SigmaC = sigma2)
rm("Rphi"); gc()
Y_or <- X_or %*% B + W_or + mniw::rMNorm(1, Lambda = matrix(0, n+u, q), SigmaR = diag((1/alfa)-1, n+u), SigmaC = sigma2)
gc()

# sample data
crd_s <- crd[1:n, ]
X <- X_or[1:n, ]
W <- W_or[1:n, ]
Y <- matrix(Y_or[1:n, ], n, q)

# prediction data
crd_u <- crd[-(1:n), ]
X_u <- X_or[-(1:n), ]
W_u <- W_or[-(1:n), ]
Y_u <- matrix(Y_or[-(1:n), ], u, q)

# # hyperparameters values
# alfa_seq <- c(0.7, 0.8, 0.9)
# phi_seq <- c(3, 4, 5)
#
# # function for the fit loop
# fit_loop <- function(i) {
#
#   Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
#   p <- ncol(Xi); q <- ncol(Yi)
#   bps <- spBPS:::BPS_weights_MvT(data = list(Y = Yi, X = Xi),
#                                priors = list(mu_B = matrix(0, nrow = p, ncol = q),
#                                              V_r = diag(10, p),
#                                              Psi = diag(1, q),
#                                              nu = 3), coords = crd_i,
#                                hyperpar = list(alpha = alfa_seq, phi = phi_seq), K = 5)
#   w_hat <- bps$W
#   epd <- bps$epd
#
#   result <- list(epd, w_hat)
#   return(result)
#
# }
#
# # function for the pred loop
# pred_loop <- function(r) {
#
#   ind_s <- subset_ind[1]
#   Ys <- data_part$Y_list[[ind_s]]; Xs <- data_part$X_list[[ind_s]]; crds <- data_part$crd_list[[ind_s]]; Ws <- W_list[[ind_s]]
#   result <- spBPS:::BPS_post_MvT(data = list(Y = Ys, X = Xs), coords = crds,
#                                X_u = X_u, crd_u = crd_u,
#                                priors = list(mu_B = matrix(0, nrow = p, ncol = q),
#                                              V_r = diag(10, p),
#                                              Psi = diag(1, q),
#                                              nu = 3),
#                                hyperpar = list(alpha = alfa_seq, phi = phi_seq),
#                                W = Ws, R = 1)
#
#   return(result)
# }
#
# # subsetting data
# subset_size <- 200
# K <- n/subset_size
# data_part <- spBPS::subset_data(data = list(Y = Y, X = X, crd = crd_s), K = K)
#
# # timing
# tic("total")
#
# # number of clusters for parallel implementation
# n.core <- parallel::detectCores(logical=F)-1
#
# # list of function
# funs_fit <- lsf.str()[which(lsf.str() != "fit_loop")]
#
# # list of function
# funs_pred <- lsf.str()[which(lsf.str() != "pred_loop")]
#
# # starting cluster
# cl <- makeCluster(n.core)
# registerDoParallel(cl)
#
# # parallelized subset computation of GP in different cores
# tic("fit")
# obj_fit <- foreach(i = 1:K, .noexport = funs_fit) %dopar% { fit_loop(i) }
# fit_time <- toc()
#
# gc()
# # Combination using double BPS
# tic("comb")
# comb_bps <- spBPS:::BPS_combine(obj_fit, K, 1)
# comb_time <- toc()
# Wbps <- comb_bps$W
# W_list <- comb_bps$W_list
#
# gc()
# # parallelized subset computation of GP in different cores
# R <- 200
# subset_ind <- sample(1:K, R, T, Wbps)
# tic("prediction")
# predictions <- foreach(r = 1:R, .noexport = funs_pred) %dopar% { pred_loop(r) }
# prd_time <- toc()
#
# # closing cluster
# stopCluster(cl)
# gc()
#
# # timing
# tot_time <- toc()


# New unified function: spBPS()

# priors
priors <- list(mu_B = matrix(0, nrow = p, ncol = q),
               V_r = diag(10, p),
               Psi = diag(1, q),
               nu = 3)

# hyperparameters values
alfa_seq <- c(0.7, 0.8, 0.9)
phi_seq <- c(3, 4, 5)
hyperpar <- list(alpha = alfa_seq, phi = phi_seq)

tic("spBPS")
out <- spBPS::spBPS(data = list(Y = Y, X = X),
# out <- spBPS(data = list(Y = Y, X = X),
      priors = priors,
      coords = crd_s,
      hyperpar = hyperpar,
      subset_size = 200,
      cv_folds = 5,
      rp = 1,
      combine_method = "bps",
      draws = 200,
      newdata = list(X = X_u, coords = crd_u),
      cores = parallel::detectCores(logical=F)-1)
new_time <- toc()

# Compare results

# Compare times
tot_time$callback_msg
new_time$callback_msg

# Compare predictions
pred_matrix_legacy <- do.call(abind::abind, c(lapply(predictions, function(x) x$Pred[[1]]$Yu), along = 3))
pred_matrix_orch <- do.call(abind::abind, c(lapply(out$predictive, function(x) x$Yu), along = 3))

# Adjusting for disagreement
disagree_array <- do.call(abind::abind, c(lapply(predictions, function(x) x$Pred[[1]]$MY), along = 3))
pred_matrix_legacy <- pred_matrix_legacy - disagree_array
agreement <- apply(disagree_array, c(1,2), mean)
agree_array <- replicate(R, agreement, simplify = "array")
pred_matrix_legacy <- pred_matrix_legacy + agree_array

# Check equality
all.equal(pred_matrix_legacy, pred_matrix_orch, tolerance = 1e-5)

# Check root mean squared prediciton error (RMSPE) and coverage
post_mean_legacy <- apply(pred_matrix_legacy, c(1,2), mean)
post_qnt_legacy <- apply(pred_matrix_legacy, c(1,2), quantile, c(0.025, 0.975))
c(mean(Y_u[,1] >= post_qnt_legacy[1,,1] & Y_u[,1] <= post_qnt_legacy[2,,1]),
                mean(Y_u[,2] >= post_qnt_legacy[1,,2] & Y_u[,2] <= post_qnt_legacy[2,,2]))
mean(post_qnt_legacy[2,,]-post_qnt_legacy[1,,])
rmspe_legacy <- sqrt( colMeans( (Y_u - post_mean_legacy)^2 ) )
mean(rmspe_legacy)

post_mean_orch <- apply(pred_matrix_orch, c(1,2), mean)
post_qnt_orch <- apply(pred_matrix_orch, c(1,2), quantile, c(0.025, 0.975))
c(mean(Y_u[,1] >= post_qnt_orch[1,,1] & Y_u[,1] <= post_qnt_orch[2,,1]),
                mean(Y_u[,2] >= post_qnt_orch[1,,2] & Y_u[,2] <= post_qnt_orch[2,,2]))
mean(post_qnt_orch[2,,]-post_qnt_orch[1,,])
rmspe_orch <- sqrt( colMeans( (Y_u - post_mean_orch)^2 ) )
mean(rmspe_orch)

# Compare weights
Wbps
out$weights_global

# Compare posterior samples
beta_smp <- sapply(1:R, function(r){predictions[[r]]$Post[[1]]$beta[1:p,]}, simplify = "array")
(post_mean_beta <-  apply(beta_smp, c(1,2), mean))
post_low_beta <- apply(beta_smp, c(1,2), quantile, c(0.025))
post_upp_beta <- apply(beta_smp, c(1,2), quantile, c(0.975))

sigma_smp <- sapply(1:R, function(r){predictions[[r]]$Post[[1]]$sigma}, simplify = "array")
(post_mean_sigma <- apply(sigma_smp, c(1,2), mean))
post_low_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.025))
post_upp_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.975))

(post_mean_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(hyperpar$alpha, hyperpar$phi)) %*% W_list[[k]]) %*% Wbps)
(post_var_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(hyperpar$alpha, hyperpar$phi)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2))

# collecting
posterior_bps <- cbind(t(sapply(1:R, function(r)matrix(beta_smp[,,r]))),
                        t(sapply(1:R, function(r)matrix(sigma_smp[,,r])))[,-3],
                        cbind(rep(post_mean_hyp[2], R), rep(post_mean_hyp[2], R)))
colnames(posterior_bps) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))

# plotting bps
library(ggplot2)
library(bayesplot)
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_bps <- mcmc_recover_intervals(posterior_bps, true_par,
                                        prob = 0.95,
                                        prob_outer = 0.95,
                                        point_est = "mean",
                                        size = 4,
                                        alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[list(1,1)]), expression(beta[list(2,1)]), expression(beta[list(3,1)]), expression(beta[list(1,2)]), expression(beta[list(2,2)]), expression(beta[list(3,2)]),
                              expression(phi[1]), expression(phi[2]),
                              expression(Sigma[list(1,1)]), expression(Sigma[list(2,1)]), expression(Sigma[list(2,2)]))) +
  ggtitle("Credible intervals BPS",
          "with posterior means, true values, and 95% credible intervals")
post_int_bps


beta_smp <- sapply(1:R, function(r){out$posterior[[r]]$beta[1:p,]}, simplify = "array")
(post_mean_beta <-  apply(beta_smp, c(1,2), mean))
post_low_beta <- apply(beta_smp, c(1,2), quantile, c(0.025))
post_upp_beta <- apply(beta_smp, c(1,2), quantile, c(0.975))

sigma_smp <- sapply(1:R, function(r){out$posterior[[r]]$sigma}, simplify = "array")
(post_mean_sigma <- apply(sigma_smp, c(1,2), mean))
post_low_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.025))
post_upp_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.975))

(post_mean_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(hyperpar$alpha, hyperpar$phi)) %*% matrix(out$weights_local[[k]])) %*% matrix(out$weights_global))
(post_var_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(hyperpar$alpha, hyperpar$phi)) %*% matrix(out$weights_local[[k]]))^2 %*% matrix(out$weights_global) - (post_mean_hyp^2))

# collecting
posterior_bps <- cbind(t(sapply(1:R, function(r)matrix(beta_smp[,,r]))),
                        t(sapply(1:R, function(r)matrix(sigma_smp[,,r])))[,-3],
                        cbind(rep(post_mean_hyp[2], R), rep(post_mean_hyp[2], R)))
colnames(posterior_bps) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))

# plotting bps
library(ggplot2)
library(bayesplot)
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_bps2 <- mcmc_recover_intervals(posterior_bps, true_par,
                                        prob = 0.95,
                                        prob_outer = 0.95,
                                        point_est = "mean",
                                        size = 4,
                                        alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[list(1,1)]), expression(beta[list(2,1)]), expression(beta[list(3,1)]), expression(beta[list(1,2)]), expression(beta[list(2,2)]), expression(beta[list(3,2)]),
                              expression(phi[1]), expression(phi[2]),
                              expression(Sigma[list(1,1)]), expression(Sigma[list(2,1)]), expression(Sigma[list(2,2)]))) +
  ggtitle("Credible intervals BPS",
          "with posterior means, true values, and 95% credible intervals")
post_int_bps2
