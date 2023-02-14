# Simulation study in the multivariate framework
# Load packages
library(MASS)
library(parallel)

# Dimension
d <- 1
# d <- 2
# d <- 3
# d <- 4
# d <- 5

# Convergence orders
alpha <- min(c(d, 4)) / (2 * d + 8)

# Sample data parameters N(0, sigma)
mu <- rep(0, d)
sigma <- diag(1, nrow = d, ncol = d)

# Simulation parameters
M <- 100
k <- 5:14
J.d2 <- matrix(1, nrow = d^2, ncol = d^2)

# Compute the empirical rates
r.n.list <- mclapply(X = k, FUN = function(j) {
  # define sample sizes
  n <- 2^j

  # initialize convergence rates
  r.n <- vector(mode = "numeric", length = M)

  # compute the H_{AMISE} for N(0, sigma) distribution
  H.amise <- (4 / (d + 2))^(2 / (d + 4)) * n^(-2 / (d + 4)) * sigma


  for (i in 1:M) {
    if (d == 1) {
      # generate sample data for each Montecarlo's replicate
      xs <- drop(mvrnorm(n = n, mu = mu, Sigma = sigma))

      # compute H.ucv
      H.ucv <- tryCatch(as.matrix(ks::hucv(xs)^2), error = function(e) NA)
    } else {
      # generate sample data for each Montecarlo's replicate
      xs <- mvrnorm(n = n, mu = mu, Sigma = sigma)

      # compute H.ucv
      H.ucv <- tryCatch(ks::Hucv.diag(xs), error = function(e) NA)
      # Hucv.diag() to speed up optimization since sigma is diagonal
    }
    # compute de convergence rate for each H.ucv and for each M
    r.n[i] <- norm(H.ucv - H.amise, type = "F") /
      norm(J.d2 %*% ks::vec(H.amise), type = "F")
  }
  r.n
}, mc.cores = 2)

# Mean rates
r.n <- do.call(rbind, r.n.list)
mean_H.ucv <- rowMeans(r.n, na.rm = TRUE)

# Congruence of powers
(mod <- lm(log(mean_H.ucv) ~ log(2^k)))
-alpha


# Compare empirical rate with theoretical rate
matplot(k, log(mean_H.ucv), type = "l", col = 1, lwd = 2, ylim = c(-3, 2))
matlines(k, log(2^(-alpha * k)), type = "l", col = "red", lwd = 2, lty = 1)
legend("topleft", legend = c("Empirical", "Theoretical"), col = 1:2, lwd = 2)

# Check whether the slopes of the above lines are equal using confidence intervals
ci_beta1 <- confint(mod, level = 0.95)
beta1 <- coef(mod)[2]
lower_bound <- ci_beta1[2, 1]
upper_bound <- ci_beta1[2, 2]
beta1
lower_bound
upper_bound
-alpha

plot(c(0, 1),
  c(0, 0),
  type = "n",
  xlab = "",
  ylab = "",
  xaxt = "n",
  ylim = c(ci_beta1[2, 1] - 0.15, ci_beta1[2, 2] + 0.15)
)
abline(h = ci_beta1[2, 1], col = "black", lty = 2)
abline(h = ci_beta1[2, 2], col = "black", lty = 2)
abline(h = -alpha, col = "red", lty = 1, lwd = 2)
legend(
  "topright",
  legend = c("CI of empirical slope", "Theoretical slope"),
  col = c("black", "red"),
  lty = c(2, 1),
  lwd = c(1, 2),
  bty = "n"
)


# Plots
n <- 2^k
col <- colorRampPalette(c("yellow", "blue", "yellow"))
matplot(k, log(r.n),
  type = "l",
  xlab = "k",
  ylim = c(-8, 2),
  col = col(12)
)
matlines(k, log(mean_H.ucv), type = "l", col = "red", lwd = 2, lty = 1)
legend(
  x = "topright",
  legend = "log(mean(H.ucv))",
  inset = 0.05,
  lty = 1,
  col = "red",
  lwd = 2,
  bty = "n",
  cex = 0.80
)
