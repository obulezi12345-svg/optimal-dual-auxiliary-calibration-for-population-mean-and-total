# ==============================================================================
# Monte Carlo Simulation Study: Dual-Auxiliary Calibration Estimation
# Author: Based on Obulezi (2026) Manuscript
# Description: Evaluates Point Estimation, Relative Efficiency, Variance
#              Estimation, Coverage Probabilities, and Power Curves.
# ==============================================================================
setwd("C:\\Users\\Dr. O. J. Obulezi\\Documents\\R projects\\calibrated weighted estimator")
# 0. Load Required Packages ----------------------------------------------------
required_pkgs <- c("ggplot2", "gridExtra", "MASS")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Set Seed for Reproducibility
set.seed(2026)

# 1. Simulation Parameters -----------------------------------------------------
N <- 10000                  # Population size
sample_sizes <- c(100, 300, 500, 1000) # Sample sizes (n)
R <- 5000                   # Monte Carlo replications
alpha_level <- 0.05         # Nominal significance level (5%)

# Population Auxiliary Characteristics
mu_x <- c(50, 60)
sigma_x1 <- 10
sigma_x2 <- 12
rho_x1x2 <- 0.5
cov_x1x2 <- rho_x1x2 * sigma_x1 * sigma_x2

Sigma_x <- matrix(c(sigma_x1^2, cov_x1x2,
                    cov_x1x2, sigma_x2^2), nrow = 2)

# Generate Finite Population
X_mat <- mvrnorm(n = N, mu = mu_x, Sigma = Sigma_x)
X1 <- X_mat[, 1]
X2 <- X_mat[, 2]

# True Population Model Coefficients
beta0 <- 15
beta1 <- 1.8
beta2 <- 2.2
sigma_eps <- 12

# Heteroscedastic error structure: Var(eps_i) = sigma_eps^2 * (X1_i / mean(X1))
eps <- rnorm(N, mean = 0, sd = sigma_eps * sqrt(X1 / mean(X1)))
Y <- beta0 + beta1 * X1 + beta2 * X2 + eps

# True Population Means and Totals
Y_mean <- mean(Y)
X1_mean <- mean(X1)
X2_mean <- mean(X2)
T_z <- c(N, N * X1_mean, N * X2_mean)

# Optimal tuning weight structure: q_i = 1 / Var(eps_i)
q_pop <- 1 / (X1 / mean(X1))

# 2. Main Simulation Loop ------------------------------------------------------
cat("Starting Monte Carlo Simulation (R =", R, "replications)...\n")

results_list <- list()

for (n in sample_sizes) {
  cat(sprintf("Running simulations for n = %d...\n", n))

  f <- n / N
  d_i <- N / n # Basic design weights under SRSWOR

  # Storage matrices for replications
  estimates <- matrix(NA, nrow = R, ncol = 3) # HT, Cal(X1), Cal(X1, X2)
  var_est   <- matrix(NA, nrow = R, ncol = 3) # Analytical variance estimates

  for (r in 1:R) {
    # Draw Sample under SRSWOR
    s_idx <- sample(1:N, size = n, replace = FALSE)

    y_s  <- Y[s_idx]
    x1_s <- X1[s_idx]
    x2_s <- X2[s_idx]
    q_s  <- q_pop[s_idx]

    # --- 1. Horvitz-Thompson Estimator ---
    y_HT <- mean(y_s)
    v_HT <- ((1 - f) / n) * var(y_s)

    # --- 2. Single Auxiliary Calibration Estimator (X1 only) ---
    z1_s <- cbind(1, x1_s)
    T_z1 <- c(N, N * X1_mean)
    Z1_HT <- colSums(d_i * z1_s)

    M1 <- t(z1_s) %*% diag(d_i * q_s) %*% z1_s
    lambda1 <- solve(M1, T_z1 - Z1_HT)
    w1 <- d_i * (1 + q_s * as.vector(z1_s %*% lambda1))

    y_cal1 <- sum(w1 * y_s) / N
    e1 <- y_s - z1_s %*% solve(M1, t(z1_s) %*% diag(d_i * q_s) %*% y_s)
    g1 <- w1 / d_i
    v_cal1 <- ((1 - f) / n) * (1 / (n - 1)) * sum((g1 * e1)^2)

    # --- 3. Dual Auxiliary Calibration Estimator (X1, X2) ---
    z2_s <- cbind(1, x1_s, x2_s)
    Z2_HT <- colSums(d_i * z2_s)

    M2 <- t(z2_s) %*% diag(d_i * q_s) %*% z2_s
    lambda2 <- solve(M2, T_z - Z2_HT)
    w2 <- d_i * (1 + q_s * as.vector(z2_s %*% lambda2))

    y_cal2 <- sum(w2 * y_s) / N
    e2 <- y_s - z2_s %*% solve(M2, t(z2_s) %*% diag(d_i * q_s) %*% y_s)
    g2 <- w2 / d_i
    v_cal2 <- ((1 - f) / n) * (1 / (n - 1)) * sum((g2 * e2)^2)

    # Store results
    estimates[r, ] <- c(y_HT, y_cal1, y_cal2)
    var_est[r, ]   <- c(v_HT, v_cal1, v_cal2)
  }

  # Calculate Metrics for Sample Size n
  mean_est <- colMeans(estimates)
  bias     <- mean_est - Y_mean
  percent_rb <- (bias / Y_mean) * 100

  emp_var  <- apply(estimates, 2, var)
  mean_vhat <- colMeans(var_est)
  mse      <- apply(estimates, 2, function(x) mean((x - Y_mean)^2))

  # Relative Efficiency wrt Horvitz-Thompson (%)
  re <- (mse[1] / mse) * 100

  # Coverage Probability (95% CI)
  z_crit <- qnorm(1 - alpha_level / 2)
  ci_lower <- estimates - z_crit * sqrt(var_est)
  ci_upper <- estimates + z_crit * sqrt(var_est)
  cp <- colMeans(ci_lower <= Y_mean & ci_upper >= Y_mean) * 100

  # Summary Data Frame
  df_summary <- data.frame(
    Sample_Size = n,
    Estimator = c("Horvitz-Thompson", "Calibrated (X1)", "Proposed Calibrated (X1, X2)"),
    Mean_Estimate = mean_est,
    Percent_RB = percent_rb,
    Empirical_Variance = emp_var,
    Analytical_Variance = mean_vhat,
    MSE = mse,
    Relative_Efficiency = re,
    Coverage_Probability = cp
  )

  results_list[[as.character(n)]] <- df_summary
}

# Combine all performance results
final_results_df <- do.call(rbind, results_list)
rownames(final_results_df) <- NULL

# Save Simulation Table to Working Directory
write.csv(final_results_df, "simulation_results_performance.csv", row.names = FALSE)
cat("Saved 'simulation_results_performance.csv' to working directory.\n")


# 3. Power Function Evaluation under Pitman Alternatives ---------------------
cat("Evaluating Power Curves under Pitman Local Alternatives...\n")

n_power <- 300 # Fixed sample size for power curves
delta_seq <- seq(-15, 15, length.out = 31)
power_results <- data.frame()

z_crit <- qnorm(1 - alpha_level / 2)
f_p <- n_power / N
d_i_p <- N / n_power

for (delta in delta_seq) {
  # True mean under local alternative
  Y_alt <- Y_mean + delta / sqrt(n_power)

  rejections <- c(HT = 0, Cal1 = 0, Cal2 = 0)

  for (r in 1:1000) { # 1000 reps per delta step
    s_idx <- sample(1:N, size = n_power, replace = FALSE)

    y_s  <- Y[s_idx]
    x1_s <- X1[s_idx]
    x2_s <- X2[s_idx]
    q_s  <- q_pop[s_idx]

    # 1. HT
    y_HT <- mean(y_s)
    v_HT <- ((1 - f_p) / n_power) * var(y_s)
    z_stat_HT <- (y_HT - Y_alt) / sqrt(v_HT)

    # 2. Cal1
    z1_s <- cbind(1, x1_s)
    T_z1 <- c(N, N * X1_mean)
    Z1_HT <- colSums(d_i_p * z1_s)
    M1 <- t(z1_s) %*% diag(d_i_p * q_s) %*% z1_s
    lambda1 <- solve(M1, T_z1 - Z1_HT)
    w1 <- d_i_p * (1 + q_s * as.vector(z1_s %*% lambda1))
    y_cal1 <- sum(w1 * y_s) / N
    e1 <- y_s - z1_s %*% solve(M1, t(z1_s) %*% diag(d_i_p * q_s) %*% y_s)
    v_cal1 <- ((1 - f_p) / n_power) * (1 / (n_power - 1)) * sum(((w1/d_i_p) * e1)^2)
    z_stat_cal1 <- (y_cal1 - Y_alt) / sqrt(v_cal1)

    # 3. Cal2
    z2_s <- cbind(1, x1_s, x2_s)
    Z2_HT <- colSums(d_i_p * z2_s)
    M2 <- t(z2_s) %*% diag(d_i_p * q_s) %*% z2_s
    lambda2 <- solve(M2, T_z - Z2_HT)
    w2 <- d_i_p * (1 + q_s * as.vector(z2_s %*% lambda2))
    y_cal2 <- sum(w2 * y_s) / N
    e2 <- y_s - z2_s %*% solve(M2, t(z2_s) %*% diag(d_i_p * q_s) %*% y_s)
    v_cal2 <- ((1 - f_p) / n_power) * (1 / (n_power - 1)) * sum(((w2/d_i_p) * e2)^2)
    z_stat_cal2 <- (y_cal2 - Y_alt) / sqrt(v_cal2)

    # Test decision
    if (abs(z_stat_HT) > z_crit)   rejections["HT"]   <- rejections["HT"] + 1
    if (abs(z_stat_cal1) > z_crit) rejections["Cal1"] <- rejections["Cal1"] + 1
    if (abs(z_stat_cal2) > z_crit) rejections["Cal2"] <- rejections["Cal2"] + 1
  }

  power_vec <- rejections / 1000

  power_results <- rbind(power_results, data.frame(
    Delta = delta,
    HT = power_vec["HT"],
    Cal1 = power_vec["Cal1"],
    Cal2 = power_vec["Cal2"]
  ))
}

# Save Power Results Table to CSV
write.csv(power_results, "simulation_results_power.csv", row.names = FALSE)
cat("Saved 'simulation_results_power.csv' to working directory.\n")


# 4. Generate Publication-Ready Plots ----------------------------------------
cat("Generating and exporting publication-ready plots...\n")

# Plot 1: Relative Efficiency vs Sample Size
p1 <- ggplot(final_results_df, aes(x = factor(Sample_Size), y = Relative_Efficiency,
                                   fill = Estimator, group = Estimator)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_manual(values = c("Horvitz-Thompson" = "#999999",
                               "Calibrated (X1)" = "#E69F00",
                               "Proposed Calibrated (X1, X2)" = "#0072B2")) +
  labs(title = "Relative Efficiency across Sample Sizes",
       subtitle = "Reference Baseline: Horvitz-Thompson = 100%",
       x = "Sample Size (n)", y = "Relative Efficiency (%)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))

# Save Plot 1
ggsave("plot_relative_efficiency.png", plot = p1, width = 8, height = 5, dpi = 300)
ggsave("plot_relative_efficiency.pdf", plot = p1, width = 8, height = 5)

# Plot 2: Empirical Power Curves under Local Alternatives
power_long <- reshape2::melt(power_results, id.vars = "Delta",
                             variable.name = "Estimator", value.name = "Power")
levels(power_long$Estimator) <- c("Horvitz-Thompson", "Calibrated (X1)", "Proposed Calibrated (X1, X2)")

p2 <- ggplot(power_long, aes(x = Delta, y = Power, color = Estimator, linetype = Estimator)) +
  geom_line(size = 1.1) +
  geom_point(size = 2) +
  geom_hline(yintercept = alpha_level, linetype = "dashed", color = "red", alpha = 0.7) +
  scale_color_manual(values = c("Horvitz-Thompson" = "#999999",
                                "Calibrated (X1)" = "#E69F00",
                                "Proposed Calibrated (X1, X2)" = "#0072B2")) +
  labs(title = "Empirical Power Curves under Pitman Local Alternatives",
       subtitle = expression(paste("Evaluating Power ", Pi(delta), " at Significance Level ", alpha, " = 0.05")),
       x = expression(paste("Local Alternative Shift Parameter (", delta, ")")),
       y = "Empirical Power Function") +
  annotate("text", x = 10, y = 0.08, label = expression(paste("Type I Error Rate (", alpha, " = 0.05)")), color = "red") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# Save Plot 2
ggsave("plot_power_curves.png", plot = p2, width = 8, height = 5, dpi = 300)
ggsave("plot_power_curves.pdf", plot = p2, width = 8, height = 5)

cat("\n======================================================================\n")
cat("Simulation Complete! The following files are saved in your working directory:\n")
cat(" 1. simulation_results_performance.csv (Estimation metrics & CPs)\n")
cat(" 2. simulation_results_power.csv       (Power grid data)\n")
cat(" 3. plot_relative_efficiency.png/pdf    (Relative efficiency bar charts)\n")
cat(" 4. plot_power_curves.png/pdf           (Asymptotic power curves)\n")
cat("======================================================================\n")
