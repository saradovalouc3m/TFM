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

# Sample data
#------------------------------------------------------
n <- c(2^(5:16))
h_opt <- (R_K_gaussian / (n * mu2_K_gaussian^2 * R_curvature))^(1 / 5)
mean.41 <- h_opt
sd.41 <- sqrt(2 * R_p * R_f / (25 * n^2 * h_opt^7 * mu2_K_gaussian^4 * R_curvature^2))


# Data reading
# -------------------------------------------------------
library(R.matlab)
data <- list.files(path = "~/Desktop/UC3M, Statistics for Data Science/TFM/R scripts/Unidimensional case/Results_simulation", pattern = "*.mat")
data <- paste0(sort(as.integer(substring(data, 1, nchar(data)-4))), ".mat")
data


# Kolmogorov-Smirnov test por each file
# -------------------------------------------------------
M <- length(data)
total_data <- list()
mean_hucv <- numeric(M)
pvalues_H1 <- numeric(M)
statistic_values <- numeric(M)
for (i in seq_len(M)) {
  # Each .mat data file
  data_i <- readMat(data[i])
  total_data <- append(total_data, list(data_i))
  mean_hucv[i] <- mean(data_i$hucv)
  
  
  # Kolmogorov-Smirnov test for H_0: F = N(mean.41, sd.41)
  pvalues_H1[i] <- ks.test(x = as.vector(data_i$hucv), y = "pnorm", mean = mean.41[i], sd = sd.41[i])$p.value
  statistic_values[i] <- ks.test(x = as.vector(data_i$hucv), y = "pnorm", mean = mean.41[i], sd = sd.41[i])$statistic
}

total_data
pvalues_H1
statistic_values
mean_hucv


# Histogram of the p-values and the value of the statistic 
# depending on the sample size
# -------------------------------------------------------
hist(pvalues_H1,
       breaks = seq(0, 1, l = 20), probability = TRUE,
       main = latex2exp::TeX("$H_1$")
  )
abline(h = 1, col = 2)


hist(statistic_values,
     probability = TRUE, main = latex2exp::TeX("$H_1$")
)


plot(k, log(statistic_values),
     ylim = c(-5.0,-1.2),
     type = "l",
     col = "lightblue", 
     lwd = 5,
     ylab = latex2exp::TeX("value of $ log(D_{n})$"))




