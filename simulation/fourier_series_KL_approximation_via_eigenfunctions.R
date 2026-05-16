library(glue)
f <- glue::glue

set.seed(42)

# setup
n <- 200          # number of subjects
M <- 200
s <- seq(0, 1, length.out = M)
ds <- s[2] - s[1]

sigma1_sq <- 1.00
sigma2_sq <- 0.25

# simulate curves
u1 <- rnorm(n, 0, sqrt(sigma1_sq))
u2 <- rnorm(n, 0, sqrt(sigma1_sq))
u3 <- rnorm(n, 0, sqrt(sigma2_sq))
u4 <- rnorm(n, 0, sqrt(sigma2_sq))

# each row is one subject's curve X_i(s)
X <- outer(u1, sin(2*pi*s)) +
  outer(u2, cos(2*pi*s)) +
  outer(u3, sin(4*pi*s)) +
  outer(u4, cos(4*pi*s))   # n x M matrix

# theoretical covariance kernel
K_theory <- sigma1_sq * outer(sin(2*pi*s), sin(2*pi*s)) +
  sigma1_sq * outer(cos(2*pi*s), cos(2*pi*s)) +
  sigma2_sq * outer(sin(4*pi*s), sin(4*pi*s)) +
  sigma2_sq * outer(cos(4*pi*s), cos(4*pi*s))

# estimate from data
K_hat <- cov(X)    # M x M

# solve eigenvalue problem
eig_theory <- eigen(K_theory * ds, symmetric = TRUE)
eig_data <- eigen(K_hat, symmetric = TRUE)

phi_theory <- eig_theory$vectors
phi_data <- eig_data$vectors
l2_norms <- sqrt(colSums(phi_data^2) * ds)    # ∫ φ(s)^2 ds ≈ 1
phi_data <- sweep(phi_data, 2, l2_norms, "/")

lam_theory <- eig_theory$values    # should be: 1, 1, 0.25, 0.25, 0, 0, ...
lam_data <- eig_data$values * ds

# compare eigenvalues
f("Theoretical eigenvalues (top 4): {paste(round(lam_theory[1:4], 4), collapse=' ')}")
f("Estimated eigenvalues  (top 4): {paste(round(lam_data[1:4], 4), collapse=' ')}")

# plot eigenfunctions: estimated vs theoretical
par(mfrow = c(2,2))
theoretical_phi <- list(
  sqrt(2)*sin(2*pi*s),
  sqrt(2)*cos(2*pi*s),
  sqrt(2)*sin(4*pi*s),
  sqrt(2)*cos(4*pi*s)
)

for(k in 1:4){
  # fix sign flip (eigenvectors defined up to ±1)
  sign_fix <- sign(sum(phi_data[,k] * theoretical_phi[[k]]))

  plot(s, theoretical_phi[[k]], type='l', col='red', lwd=2,
       main=paste("Eigenfunction", k),
       ylab=expression(phi(s)), xlab="s")
  lines(s, sign_fix * phi_data[,k], col='blue', lty=2, lwd=2)
  legend("topright", c("Theoretical","Estimated"),
         col=c("red","blue"), lty=1:2, cex=0.7)
}


##### CORRECT FOR SIGNS ########################################################
#### Uncomment after correction
# 
# # compute FPC scores
# X_centered <- sweep(X, 2, colMeans(X), "-")
# xi <- X_centered %*% phi_data[, 1:4] * ds   # n x 4 scores
# 
# # verify: scores should recover original u's
# par(mfrow=c(2,2))
# plot(u1, xi[,1], main="Score 1 vs u1 (should be diagonal)",
#      xlab="True u1", ylab="ξ_i1"); abline(0,1,col='red')
# plot(u2, xi[,2], main="Score 2 vs u2",
#      xlab="True u2", ylab="ξ_i2"); abline(0,1,col='red')
# plot(u3, xi[,3], main="Score 3 vs u3",
#      xlab="True u3", ylab="ξ_i3"); abline(0,1,col='red')
# plot(u4, xi[,4], main="Score 4 vs u4",
#      xlab="True u4", ylab="ξ_i4"); abline(0,1,col='red')
# 
# # variance explained
# total_var <- sum(lam_data[lam_data > 0])
# pct_var <- cumsum(lam_data) / total_var * 100
# f("Cumulative variance explained by top 4 FPCs: {round(pct_var[4], 1)}%")   # should be ~100%
