# Internal C source code used in stats::bw.ucv()
# https://github.com/aviralg/r/blob/9ca112c1b9b6a9751e3ba2944cca78cc3657e4c5/src/library/stats/src/bandwidths.c
# SEXP bw_ucv(SEXP sn, SEXP sd, SEXP cnt, SEXP sh)
# {
#   double h = asReal(sh), d = asReal(sd), sum = 0.0, term, u;
#   int n = asInteger(sn), nbin = LENGTH(cnt);
#   double *x = REAL(cnt);
#   for (int i = 0; i < nbin; i++) {
#     double delta = i * d / h;
#     delta *= delta;
#     if (delta >= DELMAX) break;
#     term = exp(-delta / 4.) - sqrt(8.0) * exp(-delta / 2.);
#     sum += term * x[i];
#   }
#   u = (0.5 + sum/n) / (n * h * M_SQRT_PI);
#   // = 1 / (2 * n * h * sqrt(PI)) + sum / (n * n * h * sqrt(PI));
#   return ScalarReal(u);
# }

# Translation to R
DELMAX <- 1000
ucv <- function(sn, sd, cnt, sh, mod = FALSE) {
  h <- sh
  d <- sd
  sum <- 0.0
  term <- 0
  u <- 0
  n <- sn
  nbin <- length(cnt)
  x <- cnt
  for (i in 1:nbin) {
    delta <- (i - 1) * d / h
    delta <- delta * delta
    if (!mod & delta >= DELMAX) break
    term <- exp(-delta / 4) - sqrt(8.0) * exp(-delta / 2)
    sum <- sum + term * x[i]
  }
  u <- (0.5 + sum / n) / (n * h * sqrt(pi))
  return(u)
}

# Test of the objective function
n <- 2000
x <- rnorm(n)
nb <- 1e3
Z <- stats:::bw_pair_cnts(x, nb, n > nb / 2)
d <- Z[[1L]]
cnt <- Z[[2L]]
ucv_C <- .Call(stats:::C_bw_ucv, n, d, cnt, h)
ucv_R <- function(h, mod = FALSE) ucv(n, d, cnt, h, mod)
h <- seq(0.01, 1, l = 100)
plot(h, sapply(h, ucv_C))
points(h, sapply(h, ucv_R), col = 2)
points(h, sapply(h, ucv_R, mod = TRUE), col = 3)

# Selector based on (6.67) on Scott
bw_ucv_Scott <- function(x, opt = TRUE,
                         nb = 1e4, # CAREFUL: the number of bins has to increase with n!!!
                         h_grid = 10^seq(log10(0.01), log10(1), l = 100),
                         show_plot = FALSE) {
  n <- length(x)
  Z <- stats:::bw_pair_cnts(x, nb, n > nb / 2)
  d <- Z[[1L]]
  cnt <- Z[[2L]]
  ucv_C <- function(h) .Call(stats:::C_bw_ucv, n, d, cnt, h)
  if (opt) {
    h_opt <- optim(
      par = 0.5 * n^{
        -1 / 5
      }, fn = ucv_C, lower = 1e-5,
      upper = 10 * n^{
        -1 / 5
      },
      method = "L-BFGS-B"
    )$par
  } else {
    ucv_grid <- sapply(h_grid, ucv_C)
    h_opt <- h_grid[which.min(ucv_grid)]
    if (show_plot) {
      plot(h_grid, ucv_grid)
      abline(v = h_opt, col = 2)
    }
  }
  return(h_opt)
}

# Test
bw_ucv_Scott(rnorm(1e3))

# Small simulation -- it still looks non-normalish
set.seed(1234)
n <- 1e3
M <- 1e4
h.ucv <- numeric(M)


pb <- txtProgressBar(style = 3)
for (i in 1:M) {
  set.seed(i)
  xs <- rnorm(n)
  h.ucv[i] <- bw_ucv_Scott(xs, nb = 1e3)
  setTxtProgressBar(pb, value = i / M)
}
hist(h.ucv, freq = FALSE)
rug(h.ucv)

# Important takeaways:
# - Optimizing over a grid of bandwidths gives discrete values which is not
#   adequate to estimate the distribution of h_UCV by Monte Carlo.
# - The number of bins (defaults to nb = 1000L in bw.ucv()!) needs to be
#   adapted to n. It seems that is affecting the result! We will compare what happens
#   with nb = 1e3 to nb = 1e4.
# - The distribution is highly skewed, even for n = 1e4.
