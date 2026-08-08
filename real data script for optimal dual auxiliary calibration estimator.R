# ==============================================================================
# Comparative Application of Calibration Estimators to Real-World Datasets
# Author: Based on Obulezi (2026) Theoretical Framework
# Description: Fetches and constructs 2 real-world datasets (`trees` & `swiss`),
#              implements all 6 estimators (Horvitz-Thompson, Deville & Särndal,
#              Gupta et al., Singh et al., Abubakar et al., and Proposed Optimal),
#              and computes Point Estimates, Percent Bias, Empirical Variance,
#              MSE, and Relative Efficiency (RE %).
# ==============================================================================

# 0. Load Required Packages ----------------------------------------------------
required_pkgs <- c("MASS", "stats")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

set.seed(2026)

# ==============================================================================
# Safe Inversion Helper Function
# Prevents "system is computationally singular" errors using Moore-Penrose ginv
# ==============================================================================
safe_solve <- function(M, rhs) {
  tryCatch({
    solve(M, rhs)
  }, error = function(e) {
    # Fallback to Generalized Inverse (MASS::ginv) if matrix is singular/ill-conditioned
    MASS::ginv(M) %*% rhs
  })
}

# ==============================================================================
# Core Function: Compute All 6 Calibration Estimators for a Sample
# ==============================================================================
run_calibration_estimators <- function(data_df, y_col, x1_col, x2_col, strat_col = NULL, n_sample = 15) {
  N <- nrow(data_df)
  Y_pop  <- data_df[[y_col]]
  X1_pop <- data_df[[x1_col]]
  X2_pop <- data_df[[x2_col]]

  # Population Totals & Means
  Y_mean  <- mean(Y_pop)
  X1_mean <- mean(X1_pop)
  X2_mean <- mean(X2_pop)
  T_z     <- c(N, N * X1_mean, N * X2_mean)

  # Heteroscedastic weight structure: q_i = 1 / (X1_i / mean(X1))
  q_pop <- 1 / (X1_pop / mean(X1_pop))

  # Draw SRSWOR Sample
  s_idx <- sample(1:N, size = n_sample, replace = FALSE)

  y_s  <- Y_pop[s_idx]
  x1_s <- X1_pop[s_idx]
  x2_s <- X2_pop[s_idx]
  q_s  <- q_pop[s_idx]

  d_i <- N / n_sample  # Basic design weights

  # ----------------------------------------------------------------------------
  # 1. Unadjusted Baseline Estimator (Horvitz & Thompson, 1952)
  # ----------------------------------------------------------------------------
  y_HT <- mean(y_s)

  # ----------------------------------------------------------------------------
  # 2. Classical Univariate Calibration Estimator (Deville & Särndal, 1992)
  # ----------------------------------------------------------------------------
  z1_s  <- cbind(1, x1_s)
  T_z1  <- c(N, N * X1_mean)
  Z1_HT <- colSums(d_i * z1_s)

  M1      <- t(z1_s) %*% diag(d_i * q_s) %*% z1_s
  lambda1 <- safe_solve(M1, T_z1 - Z1_HT)
  w1      <- d_i * (1 + q_s * as.vector(z1_s %*% lambda1))
  y_cal1  <- sum(w1 * y_s) / N

  # ----------------------------------------------------------------------------
  # 3. Bivariate Generalized Calibration Estimator (Gupta et al., 2022)
  # ----------------------------------------------------------------------------
  z2_s  <- cbind(1, x1_s, x2_s)
  Z2_HT <- colSums(d_i * z2_s)
  M2    <- t(z2_s) %*% diag(d_i * q_s) %*% z2_s

  lambda_g <- safe_solve(M2, T_z - Z2_HT)
  w_gcal   <- d_i * (1 + q_s * as.vector(z2_s %*% lambda_g))
  y_Gcal   <- sum(w_gcal * y_s) / N

  # ----------------------------------------------------------------------------
  # 4. Robust Dual Auxiliary Calibration Estimator (Singh et al., 2023)
  # ----------------------------------------------------------------------------
  residuals_s <- y_s - lm(y_s ~ x1_s + x2_s)$fitted.values
  mad_res     <- median(abs(residuals_s - median(residuals_s))) + 1e-6
  u           <- residuals_s / mad_res
  c_huber     <- 1.345
  w_huber     <- ifelse(abs(u) <= c_huber, 1, c_huber / abs(u))

  q_s_robust <- q_s * w_huber
  M_robust   <- t(z2_s) %*% diag(d_i * q_s_robust) %*% z2_s
  lambda_R   <- safe_solve(M_robust, T_z - Z2_HT)
  w_Rcal     <- d_i * (1 + q_s_robust * as.vector(z2_s %*% lambda_R))
  y_Rcal     <- sum(w_Rcal * y_s) / N

  # ----------------------------------------------------------------------------
  # 5. Stratified Dual Auxiliary Calibration Estimator (Abubakar et al., 2024)
  # ----------------------------------------------------------------------------
  if (!is.null(strat_col)) {
    strat_var <- data_df[[strat_col]]
    h_levels  <- unique(strat_var)

    y_Scal_h <- c()
    W_h_vec  <- c()

    for (h in h_levels) {
      idx_h   <- which(strat_var == h)
      N_h     <- length(idx_h)
      W_h     <- N_h / N

      idx_s_h <- intersect(s_idx, idx_h)
      n_h     <- length(idx_s_h)

      if (n_h >= 3) {
        d_ih <- N_h / n_h
        y_sh <- Y_pop[idx_s_h]
        z_sh <- cbind(1, X1_pop[idx_s_h], X2_pop[idx_s_h])
        q_sh <- q_pop[idx_s_h]

        T_zh   <- c(N_h, N_h * mean(X1_pop[idx_h]), N_h * mean(X2_pop[idx_h]))
        Z_h_HT <- colSums(d_ih * z_sh)
        M_h    <- t(z_sh) %*% diag(d_ih * q_sh) %*% z_sh

        lam_h <- safe_solve(M_h, T_zh - Z_h_HT)
        w_h   <- d_ih * (1 + q_sh * as.vector(z_sh %*% lam_h))

        y_Scal_h <- c(y_Scal_h, sum(w_h * y_sh) / N_h)
        W_h_vec  <- c(W_h_vec, W_h)
      }
    }

    if (length(y_Scal_h) > 0) {
      y_Scal <- sum(W_h_vec * y_Scal_h) / sum(W_h_vec)
    } else {
      y_Scal <- y_HT
    }
  } else {
    # If no explicit strata column exists, stratify by X1 median split
    strat_var <- ifelse(X1_pop >= median(X1_pop), 1, 2)
    data_df$strat_tmp <- strat_var
    y_Scal <- run_calibration_estimators(data_df, y_col, x1_col, x2_col, "strat_tmp", n_sample)$y_Scal
  }

  # ----------------------------------------------------------------------------
  # 6. Proposed Optimal Dual Auxiliary Calibration Estimator (Obulezi, 2026)
  # ----------------------------------------------------------------------------
  B_q <- safe_solve(M2, t(z2_s) %*% diag(d_i * q_s) %*% y_s)
  y_prop_cal <- y_HT + (1 / N) * t(T_z - Z2_HT) %*% B_q

  return(list(
    y_HT       = y_HT,
    y_cal1     = y_cal1,
    y_Gcal     = y_Gcal,
    y_Rcal     = y_Rcal,
    y_Scal     = as.numeric(y_Scal),
    y_prop_cal = as.numeric(y_prop_cal),
    Y_true     = Y_mean
  ))
}

# ==============================================================================
# Monte Carlo Evaluation Wrapper
# ==============================================================================
evaluate_real_dataset <- function(dataset_df, y_col, x1_col, x2_col, dataset_name, sample_size = 15, R = 2000) {
  cat("\n======================================================================\n")
  cat(sprintf("Evaluating Real Dataset: %s (N = %d, Sample n = %d, R = %d)\n",
              dataset_name, nrow(dataset_df), sample_size, R))
  cat("======================================================================\n")

  Y_true <- mean(dataset_df[[y_col]])

  # Storage Matrix
  sim_results <- matrix(NA, nrow = R, ncol = 6)
  colnames(sim_results) <- c("Horvitz-Thompson", "Classical Univariate (Deville 1992)",
                             "Bivariate Generalized (Gupta 2022)", "Robust Dual (Singh 2023)",
                             "Stratified Dual (Abubakar 2024)", "Proposed Optimal Dual (Obulezi 2026)")

  for (r in 1:R) {
    res <- run_calibration_estimators(dataset_df, y_col, x1_col, x2_col, n_sample = sample_size)
    sim_results[r, ] <- c(res$y_HT, res$y_cal1, res$y_Gcal, res$y_Rcal, res$y_Scal, res$y_prop_cal)
  }

  # Compute Performance Metrics
  mean_estimates <- colMeans(sim_results)
  bias           <- mean_estimates - Y_true
  percent_rb     <- (bias / Y_true) * 100
  emp_var        <- apply(sim_results, 2, var)
  mse            <- apply(sim_results, 2, function(x) mean((x - Y_true)^2))
  re_percent     <- (mse[1] / mse) * 100  # Relative Efficiency wrt HT Baseline

  # Summary Data Frame
  summary_df <- data.frame(
    Estimator           = colnames(sim_results),
    Mean_Estimate       = round(mean_estimates, 4),
    Percent_RB          = round(percent_rb, 4),
    Empirical_Variance  = round(emp_var, 4),
    MSE                 = round(mse, 4),
    Relative_Efficiency = round(re_percent, 2)
  )

  cat(sprintf("\nTrue Population Mean (Y_bar) = %.4f\n\n", Y_true))
  print(summary_df, row.names = FALSE)
  return(summary_df)
}

# ==============================================================================
# Application to 2 Real Datasets
# ==============================================================================

# ------------------------------------------------------------------------------
# DATASET 1: Timber Volume Data (`trees`)
# Y  = Volume (Timber Volume in cubic feet)
# X1 = Girth (Tree diameter in inches)
# X2 = Height (Tree height in feet)
# ------------------------------------------------------------------------------
data("trees")
results_trees <- evaluate_real_dataset(
  dataset_df   = trees,
  y_col        = "Volume",
  x1_col       = "Girth",
  x2_col       = "Height",
  dataset_name = "Timber Volume Dataset (trees)",
  sample_size  = 12,
  R            = 2000
)

# ------------------------------------------------------------------------------
# DATASET 2: Swiss Fertility Data (`swiss`)
# Y  = Fertility (Ig standardized fertility measure)
# X1 = Education (% education beyond primary school)
# X2 = Agriculture (% males involved in agriculture)
# ------------------------------------------------------------------------------
data("swiss")
results_swiss <- evaluate_real_dataset(
  dataset_df   = swiss,
  y_col        = "Fertility",
  x1_col       = "Education",
  x2_col       = "Agriculture",
  dataset_name = "Swiss Fertility and Socioeconomic Dataset (swiss)",
  sample_size  = 20,
  R            = 2000
)
