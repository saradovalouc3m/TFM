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
  pvalues_H1[j-4]<- ks.test(x = x.mean, y = "pnorm", mean=1, sd=1/sqrt(2^j))$p.value
  statistic_values[j-4]<- ks.test(x = x.mean, y = "pnorm", mean=1, sd=1/sqrt(2^j))$statistic
}

pvalues_H1
statistic_values
# [1] 0.02205273 0.01914100 0.02984804 0.03438274 0.03375189 0.04589654 0.03634512 0.02641757
# [9] 0.03835536 0.03661112 0.04670277 0.01936466


plot(k, log(statistic_values), 
     type="l",
     col = "lightblue", 
     lwd = 5,
     ylab = latex2exp::TeX("value of $ log(D_{n})$"))
