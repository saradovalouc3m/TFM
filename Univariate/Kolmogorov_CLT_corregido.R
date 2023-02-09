# Simulation of sample data
# --------------------------------------------
M <- 1e4
k=1:16
statistic_values <- numeric(length(k))
for (j in k){
  x.mean <- numeric(M)
  n=2^j
  for (i in 1:M) {
    set.seed(i)
    x <- rexp(n)
    x.mean[i] <- mean(x)
  }
  #Kolmogorov-Smirnov test (H_0: N(0,1) - TCL)
  # --------------------------------------------
  statistic_values[j]<- ks.test(x = x.mean, y = "pnorm",
                                mean=1, sd=1/sqrt(n))$statistic
}

plot(k, log(statistic_values),
     type="l",
     col = "lightblue",
     lwd = 5,
     xlim=c(1,16),
     ylab = latex2exp::TeX("value of $log(D_{n})$"))

