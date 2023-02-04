# 0) Calcular H_{AMISE} para una vMF(0, Sigma).
# 1) Para diferentes valores de n = 2^k, k = 5, ..., 14, generar M = 500 muestras de una vMF(0, Sigma) (usar rotasym::r_VMF())
#   (si ves que tarda mucho, usa M = 200).
# 2) Para cada una de las muestras, calcular \hat{H}_{UCV} (identificación para el caso direccional: H=h^{2}*I_{d+1})
# 3) Para cada \hat{H}_{UCV}, calcular r(n). Es decir, tendríamos r_1(n), ..., r_M(n) (tasas correspondientes a la dimensión d de S^{d} \subset R^{d+1}).
# 4) Calcular la ventana del AMISE óptima en el caso vMF: se puede modificar DirStats::bw_dir_rot() para pasarle el kappa real (en vez de que la estime),
#   mientras que para calcular el LSCV se puede usar DirStats::bw_dir_lscv()
# 5) Calcular C_aprox = mean(r_1(N), ..., r_M(N)) * N^{\alpha} para N = 2^14 (el tamaño muestral más grande que explores).
# 6) Representamos para cada i = 1, ..., M, las trayectorias n -> r_i(n).
#   Usa transparencias en los colores de las curvas para no saturar el gráfico, pues vas a graficar M curvas.
# 7) Representar las curvas n -> - C_aprox * n^{-\alpha} y n -> C_aprox * n^{-\alpha} en rojo. Deberían contener al grueso de las curvas en 5).


# Packages
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

  return((num / den) ^ (1 / (q + 4)))

}

direc_as.matrix <- function(h,d) {
  H=h^{2}*diag(1,nrow=d+1,ncol=d+1)
  return(H)
}



# Simulation
d <- 1
# d <- 2

k <- 5:10
M <- 5000

# Simulation
mu <- c(replicate(d, 0), 1)
kappa <- 2
J.d2 <- matrix(1, nrow = (d+1)^2, ncol = (d+1)^2)


r.n.list <- mclapply(X = k, FUN = function(j) {
  n=2^j
  r.n = vector(mode="numeric", length=M)
  pb <- txtProgressBar(min = 0,      
                       max = M, 
                       style = 3,    
                       width = 50,   
                       char = "=")   
  for (i in 1:M) {
    set.seed(i)
    xs <- rotasym::r_vMF(n = n, mu = mu, kappa = kappa)
    bw.amise.dir <- bw_dir_rot(data = xs, kappa = kappa)
    bw.ucv.dir <- tryCatch(DirStats::bw_dir_lscv(
      data = xs, optim_par = bw.amise.dir, optim_lower = 0.1 * bw.amise.dir,
      optim_upper = 10 * bw.amise.dir, plot_it = TRUE)$h_opt,
      error = function(e) NA)
    H.amise <- direc_as.matrix(h = bw.amise.dir, d = d)
    H.ucv <- direc_as.matrix(h = bw.ucv.dir, d = d)
    r.n[i] <- norm(H.ucv-H.amise,type="F") / norm(J.d2%*%ks::vec(H.amise),type="F")
    
    setTxtProgressBar(pb, i)
  }
  r.n
}, mc.cores = 2)

# Mean rates
r.n <- do.call(rbind, r.n.list)
mean_H.ucv <- rowMeans(r.n, na.rm = TRUE)

# Congruence of powers
(mod <- lm(log(mean_H.ucv) ~ k))


# Compare empirical rate with theoretical rate
alpha <- min(c(d, 4)) / (2 * d + 8)
matplot(k, log(mean_H.ucv), type = "l", col = 1, lwd = 2, ylim = c(-5, 2))
matlines(k, log(2^(-alpha * k)), type = "l", col = "red", lwd = 2, lty = 1)
legend("topleft", legend = c("Empirical", "Theoretical"), col = 1:2, lwd = 2)


# Plots
n=2^k
col <- colorRampPalette(c("yellow", "blue", "yellow"))


matplot(k,log(r.n),
        type = "l", 
        xlab="k",
        col=col(12))
matlines(k,log(mean_H.ucv), type = "l", col = "red", lwd = 2, lty = 1)  # Empirical rate
legend(x = "bottomleft",     
       legend = "log(mean(h.lscv))", 
       inset = 0.05,
       lty = 1,          
       col = "red",          
       lwd = 2,
       bty = "n",
       cex = 0.80) 
