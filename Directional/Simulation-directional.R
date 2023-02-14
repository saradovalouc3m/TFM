# Simulation study in the directional framework
# Load packages
library(MASS)
library(parallel)
library(rotasym)

# Functions
bw_dir_rot <- function(data, kappa) {
  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Set parameters
  n <- nrow(data)
  q <- ncol(data) - 1
  # kappa <- kappa_ml(data = data)

  # Cases for different dimensions
  if (q == 1) {
    num <- 4 * sqrt(pi) * besselI(nu = 0, x = kappa)^2
    den <- kappa * (2 * besselI(nu = 1, x = 2 * kappa) +
      3 * kappa * besselI(nu = 2, x = 2 * kappa)) * n
  } else if (q == 2) {
    num <- 8 * sinh(kappa)^2
    den <- (-2 * kappa * cosh(2 * kappa) +
      (1 + 4 * kappa^2) * sinh(2 * kappa)) * n
    # Caution! Typo in G-P (2013) in equation (6) for q = 2. It displays an
    # incorrect extra kappa:
    # * WRONG:    kappa * (1 + 4 * kappa^2) * sinh(2 * kappa)
    # * CORRECT:          (1 + 4 * kappa^2) * sinh(2 * kappa)
  } else {
    num <- 4 * sqrt(pi) * besselI(nu = (q - 1) / 2, x = kappa)^2
    den <- kappa^((q + 1) / 2) *
      (2 * q * besselI(nu = (q + 1) / 2, x = 2 * kappa) +
        (2 + q) * kappa * besselI(nu = (q + 3) / 2, x = 2 * kappa)) * n
  }

  return((num / den)^(1 / (q + 4)))
}


direc_as.matrix <- function(h, d) {
  H <- h^
    {
      2
    } * diag(1, nrow = d + 1, ncol = d + 1)
  return(H)
}


# Dimension
d <- 1
# d <- 2

# Conjetured convergence orders
alpha <- min(c(d, 4)) / (2 * d + 8)

# Sample data parameters vMF(0, sigma)
mu <- c(replicate(d, 0), 1)
kappa <- 2

# Simulation parameters
k <- 5:10
M <- 5000
J.d2 <- matrix(1, nrow = (d + 1)^2, ncol = (d + 1)^2)

# Compute the empirical rates
r.n.list <- mclapply(X = k, FUN = function(j) {
  # define sample sizes
  n <- 2^j

  # initialize convergence rates
  r.n <- vector(mode = "numeric", length = M)

  # initialising time progress bar
  pb <- txtProgressBar(
    min = 0,
    max = M,
    style = 3,
    width = 50,
    char = "="
  )
  for (i in 1:M) {
    # generate sample data for each Montecarlo's replicate
    set.seed(i)
    xs <- rotasym::r_vMF(n = n, mu = mu, kappa = kappa)

    # compute H_AMISE for vMF(0, sigma)
    bw.amise.dir <- bw_dir_rot(data = xs, kappa = kappa)
    H.amise <- direc_as.matrix(h = bw.amise.dir, d = d)

    # compute H.ucv (H=h^{2}*I_{d+1})
    bw.ucv.dir <- tryCatch(
      DirStats::bw_dir_lscv(
        data = xs, optim_par = bw.amise.dir, optim_lower = 0.1 * bw.amise.dir,
        optim_upper = 10 * bw.amise.dir, plot_it = TRUE
      )$h_opt,
      error = function(e) NA
    )
    H.ucv <- direc_as.matrix(h = bw.ucv.dir, d = d)

    # compute de convergence rate for each H.ucv and for each M
    r.n[i] <- norm(H.ucv - H.amise, type = "F") / norm(J.d2 %*% ks::vec(H.amise), type = "F")

    # ending time progress bar
    setTxtProgressBar(pb, i)
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
matplot(k, log(mean_H.ucv), type = "l", col = 1, lwd = 2, ylim = c(-5, 2))
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
  col = col(12)
)
matlines(k, log(mean_H.ucv), type = "l", col = "red", lwd = 2, lty = 1) # Empirical rate
legend(
  x = "bottomleft",
  legend = "log(mean(h.lscv))",
  inset = 0.05,
  lty = 1,
  col = "red",
  lwd = 2,
  bty = "n",
  cex = 0.80
)
