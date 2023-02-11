# Simulation of sample data
# --------------------------------------------
M <- 1e4
<<<<<<< HEAD
k = 1:16
statistic_values_tcl <- numeric(length(k))
=======
k=1:16
statistic_values <- numeric(length(k))
>>>>>>> 0c9ff93b8476f4637f93faed6181691a5eae09aa
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
<<<<<<< HEAD
  statistic_values_tcl[j]<- ks.test(x = x.mean, y = "pnorm",
                                    mean=1, sd=1/sqrt(n))$statistic
}

plot(k, log(statistic_values_tcl),
     ylim = c(-5.0,0),
=======
  statistic_values[j]<- ks.test(x = x.mean, y = "pnorm",
                                mean=1, sd=1/sqrt(n))$statistic
}

plot(k, log(statistic_values),
>>>>>>> 0c9ff93b8476f4637f93faed6181691a5eae09aa
     type="l",
     col = "lightblue",
     lwd = 5,
     xlim=c(1,16),
     ylab = latex2exp::TeX("value of $log(D_{n})$"))

