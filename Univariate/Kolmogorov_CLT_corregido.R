# Simulation of sample data
# --------------------------------------------
k=5:16
x.mean <- numeric(M)
pvalues_H1 <- numeric(length(k))
statistic_values <- numeric(length(k))
for (j in k){
  for (i in 1:M) {
    set.seed(i)
    n=2^j
    x <- rexp(n)
    x.mean[i] <- mean(x)
  }
  #Kolmogorov-Smirnov test (H_0: N(0,1) - TCL)
  # --------------------------------------------
  pvalues_H1[j-4]<- ks.test(x = x.mean, y = "pnorm", mean=1, sd=1/2^j)$p.value
  statistic_values[j-4]<- ks.test(x = x.mean, y = "pnorm", mean=1, sd=1/2^j)$statistic
}

pvalues_H1
statistic_values
# [1] 3.447732e-02 1.542118e-03 8.845586e-07 2.687046e-03 1.263373e-02 2.076407e-04 2.098753e-04 8.838619e-06 2.075903e-04 2.075903e-04 2.075903e-04 2.075903e-04
# [9] 0.3938482 0.5195456 0.7176259 0.4999082 0.4392366 0.5833262 0.5830114 0.6666667 0.5833333 0.5833333 0.5833333 0.5833333


plot(k, log(statistic_values), 
     ylim = c(-5.0,0),
     type="l",
     col = "lightblue", 
     lwd = 5,
     ylab = latex2exp::TeX("$log(D_{n})$"))

