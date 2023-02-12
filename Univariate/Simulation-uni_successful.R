# Proof of the asymtotic normality of the h.ucv bandwidth
# Function to calculate h.ucv (Scott, 2015)
ucv <- function(h, data) {
  n <- length(data)
  X <- replicate(n, data)
  delta <- ((X - t(X)) / h)
  delta <- delta * delta
  term1 <- 1 / (2 * n * h * sqrt(pi))
  total <- 0
  sum <- 1 / (n^2 * h * sqrt(pi))
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        term2 <- exp(-delta[i, j] / 4) - sqrt(8.0) * exp(-delta[i, j] / 2)
        total <- total + sum * term2
      }
    }
  }
  return(term1 + total)
}


# Optimization of the above function
ucv2 <- function(h, data) {
  n <- length(data)
  exp_mdelta <- exp(-dist(data)^2 / (2 * h^2))
  s2 <- sqrt(8.0) * sum(exp_mdelta)
  exp_mdelta <- sqrt(exp_mdelta)
  s1 <- sum(exp_mdelta)
  return(1 / (2 * n * h * sqrt(pi)) * (1 + 2 * (s1 - s2) / n))
}


# Improvements over x10
a <- rnorm(1e4)
microbenchmark::microbenchmark(ucv(0.15, a), ucv2(0.15, a), times = 10)
# Unit: milliseconds
# expr       min        lq      mean    median        uq      max neval cld
# ucv(0.15, a) 8511.0627 8545.2277 8649.8559 8606.8295 8768.3198 8861.712    10   b
# ucv2(0.15, a)  532.0203  547.8815  682.3954  581.2923  714.4529 1340.112    10  a


# Optimize
optimize_h_opt <- function(x) {
  n <- length(x)
  h_opt <- optim(
    par = 0.5 * n^{
      -1 / 5
    }, fn = function(h) ucv2(h, data = x), lower = 1e-5,
    upper = 10 * n^{
      -1 / 5
    },
    method = "L-BFGS-B"
  )$par
  return(h_opt)
}


# Small simulation of the distribution of h.ucv
n <- 32
M <- 1000
h.ucv <- numeric(M)
for (i in 1:M) {
  set.seed(i)
  xs <- rnorm(n)
  h.ucv[i] <- optimize_h_opt(xs)
}
hist(h.ucv, freq = FALSE, main = "n=32", xlab = "h.ucv")
rug(h.ucv)


# Load the .mat files containing the results of the simulations for large sample sizes
library(R.matlab)
data <- list.files(path = "~/Desktop/UC3M, Statistics for Data Science/TFM/R scripts/Unidimensional case", pattern = "*.mat")
data <- paste0(sort(as.integer(substring(data, 1, nchar(data) - 4))), ".mat")
data
xlabs <- c("n=32", "n=128", "n=1024", "n=65536")

# Simulation with M=1e3 Montecarlo's replicates
M <- length(data)
total_data <- list()
mean_hucv <- numeric(M)
par(mfrow = c(2, 2))
for (i in seq_len(M)) {
  # each .mat data file
  data_i <- readMat(data[i])
  total_data <- append(total_data, list(data_i))

  # compute mean(h.ucv) for each Montecarlo's replicate
  mean_hucv[i] <- mean(total_data[[i]]$hucv)

  # histogram
  hist(data_i$hucv, freq = FALSE, main = xlabs[i], xlab = "h.ucv")
  rug(data_i$hucv)
}

# Construct the theorical bandwidth
K_gaussian <- function(t) 1 / (sqrt(2 * pi)) * exp(-t^2 / 2)
R_K_gaussian <- integrate(function(t) K_gaussian(t)^2, lower = -Inf, upper = Inf)$value
mu2_K_gaussian <- integrate(function(t) t^2 * K_gaussian(t), lower = -Inf, upper = Inf)$value
K_gaussian_resc <- function(t, sigma) K_gaussian(t / sigma) / sigma
library(Deriv)
p <- function(c) c * Deriv(K_gaussian_resc, "t")(c, sqrt(2)) - 2 * c * Deriv(K_gaussian_resc, "t")(c, 1)
R_p <- integrate(function(c) p(c)^2, lower = -2, upper = 2)$value
sigma <- 1
R_curvature <- 3 / (8 * pi^0.5 * sigma^5)
R_f <- 1 / (2 * pi^0.5 * sigma^2)
k <- 5:16
n <- 2^k
h_opt <- (R_K_gaussian / (n * mu2_K_gaussian^2 * R_curvature))^(1 / 5)

# Compare both bandwidths
plot(k, mean_hucv,
  type = "l",
  ylab = "",
  col = "red",
  lwd = 2
)
lines(k, h_opt,
  col = "black",
  lwd = 2,
  lty = 2
)
legend(
  x = "topright",
  # inset = 0.05,
  legend = c(latex2exp::TeX("$mean(\\hat{h}_{UCV})$"), latex2exp::TeX("$h_{AMISE}$")),
  lty = c(1, 2),
  inset = 0.05,
  bty = "n",
  col = c("red", "black"),
  lwd = c(2, 1)
)
