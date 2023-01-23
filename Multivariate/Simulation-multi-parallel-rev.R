library(MASS)
library(parallel)
# 0) Calcular H_{AMISE} para una N(0, Sigma).
# 1) Para diferentes valores de n = 2^k, k = 5, ..., 14, generar M = 500 muestras de una N(0, Sigma) (si ves que tarda mucho, usa M = 200).
# 2) Para cada una de las muestras, calcular \hat{H}_{UCV}.
# 3) Para cada \hat{H}_{UCV}, calcular r(n). Es decir, tendríamos r_1(n), ..., r_M(n).
# 4) Calcular C_aprox = mean(r_1(N), ..., r_M(N)) * N^{\alpha} para N = 2^14 (el tamaño muestral más grande que explores).
# 5) Representamos para cada i = 1, ..., M, las trayectorias n -> r_i(n). Usa transparencias en los colores de las curvas para no saturar el gráfico, pues vas a graficar M curvas.
# 6) Representar las curvas n -> - C_aprox * n^{-\alpha} y n -> C_aprox * n^{-\alpha} en rojo. Deberían contener al grueso de las curvas en 5).

# 7) CORRECCIONES: La ventana H_AMISE no tienes que calcularla para cada muestra que simulas.
# Puedes calcular su expresión exacta, simplemente considerando Sigma, en lugar de "Cov(..)".
# Puedes sacarla del bucle, y calcularla en el paso 0, tal y como te habíamos indicado.
# Además, puedes tomar matrices diagonales (agilizará los cálculos)
# El caso d=1 está contenido como caso particular.
# El operador vec() se puede aplicar a un vector (matriz d×1) o escalar (matriz 1×1).
# La norma de Frobenious igualmente se puede aplicar a un vector y es simplememente la norma euclídea.
# En R tendrás que tener cuidado con los atributos de vector/matriz yadaptar ligeramente el código, pero matemáticamente es un caso particular.

d <- 1
# d <- 2
# d <- 3
# d <- 4
# d <- 5
alpha <- min(c(d, 4)) / (2 * d + 8)
M <- 100
k <- 5:14

mu <- rep(0, d)
sigma <- diag(1, nrow = d, ncol = d)
J.d2 <- matrix(1, nrow = d^2, ncol = d^2)

# Empirical rates
r.n.list <- mclapply(X = k, FUN = function(j) {
  n <- 2^j
  r.n <- vector(mode = "numeric", length = M)
  H.amise <- (4 / (d + 2))^(2 / (d + 4)) * n^(-2 / (d + 4)) * sigma
  for (i in 1:M) {
    if (d == 1) {
      xs <- drop(mvrnorm(n = n, mu = mu, Sigma = sigma))
      H.ucv <- tryCatch(as.matrix(ks::hucv(xs)), error = function(e) NA)
    } else {
      xs <- mvrnorm(n = n, mu = mu, Sigma = sigma)
      H.ucv <- tryCatch(ks::Hucv.diag(xs), error = function(e) NA)
      # Hucv.diag() to speed up optimization since sigma is diagonal
    }
    r.n[i] <- norm(H.ucv - H.amise, type = "F") /
      norm(J.d2 %*% ks::vec(H.amise), type = "F")
  }
  r.n
# }, mc.cores = 4)
}, mc.cores = 2)
  
  
# Mean rates
r.n <- do.call(rbind, r.n.list)
mean_H.ucv <- rowMeans(r.n, na.rm = TRUE)

# Congruence of powers
(mod <- lm(log(mean_H.ucv) ~ k))
-alpha

# Compare empirical rate with theoretical rate
matplot(k, log(mean_H.ucv), type = "l", col = 1, lwd = 2, ylim = c(-3, 2))
matlines(k, log(2^(-alpha * k)), type = "l", col = "red", lwd = 2, lty = 1)
legend("topleft", legend = c("Empirical", "Theoretical"), col = 1:2, lwd = 2)


# Plots
n=2^k
col <- colorRampPalette(c("yellow", "blue", "yellow"))


matplot(k,log(r.n),
        type = "l", 
        xlab="k", 
        ylim=c(-6,2),
        col=col(12))
matlines(k,log(mean_H.ucv), type = "l", col = "red", lwd = 2, lty = 1) 
legend(x = "bottomright",     
       legend = "log(mean(H.ucv))", 
       inset = 0.05,
       lty = 1,          
       col = "red",          
       lwd = 2,
       bty = "n",
       cex = 0.80)  

