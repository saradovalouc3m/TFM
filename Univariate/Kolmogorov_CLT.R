# Simulation of exponential sample data
M <- 1e4
k <- 1:16
statistic_values_tcl <- numeric(length(k))
for (j in k) {
  x.mean <- numeric(M)
  n <- 2^j
  for (i in 1:M) {
    set.seed(i)

    # generate sample data
    x <- rexp(n)

    # compute the mean of M replicates for each sample size
    x.mean[i] <- mean(x)
  }

  # Kolmogorov-Smirnov test (H_0: N(0,1) - TCL)
  statistic_values_tcl[j] <- ks.test(
    x = x.mean, y = "pnorm",
    mean = 1, sd = 1 / sqrt(n)
  )$statistic
}

# Evolution of the statistic with increasing sample size
plot(k, log(statistic_values_tcl),
  ylim = c(-5.0, 0),
  type = "l",
  col = "lightblue",
  lwd = 5,
  ylab = latex2exp::TeX("$log(D_{n})$")
)
