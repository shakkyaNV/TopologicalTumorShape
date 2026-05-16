# ## Permutation invariance of the persistence image pipeline
#
# This script demonstrates algebraically that permuting the *column order* of
# the persistence image matrix X has no effect on Cox model estimates or PC
# scores, and that the spatial beta back-projection is fully recoverable.
#
# ---
#
# ### What is a permutation matrix P?
#
# P is an orthogonal matrix (P'P = PP' = I) whose columns are a reordering of
# the identity columns. Multiplying X on the right by P reorders the *columns*
# of X without changing any individual value:
#
#     X_perm = X · P
#
# In R: `X_perm = X[, perm_idx]` where `perm_idx = sample(ncol(X))`.
#
# ---
#
# ### Why are PC scores invariant under column permutation?
#
# Let `X = U · D · V'` be the thin SVD of the centred X (shape n×p).
#
# The SVD of `X_perm = X · P` is:
#
#     X_perm = U · D · (P'V)'
#
# so the right singular vectors of X_perm are `V_perm = P'V`
# (P merely permutes the *rows* of V, which correspond to pixel positions).
#
# PC scores (projection of each subject onto the first k eigenfunctions):
#
#     scores_perm = X_perm · V_perm[:,1:k]
#                = X · P · (P'V)[:,1:k]
#                = X · (P · P') · V[:,1:k]
#                = X · V[:,1:k]          (since P'P = I)
#                = scores_orig
#
# The scores are *identical* to those from the original X.
# This holds for any orthogonal P, regardless of dimension.
#
# ---
#
# ### Why is the raw beta back-projection spatially scrambled?
#
# The Cox coefficient `cox_coef` is a scalar that weights the first PC score.
# The coefficient function (beta on the birth-death plane) is:
#
#     beta_orig     = V[:,1:k]     · cox_coef       (p-vector)
#     beta_perm_raw = V_perm[:,1:k] · cox_coef_perm  (same p-vector, rows in
#                   = (P'V)[:,1:k] · cox_coef         permuted order)
#
# Because V_perm = P'V permutes the *rows* of V (= pixel positions), the
# elements of `beta_perm_raw` are in the wrong spatial order — the coefficient
# at pixel position j in the original layout is now at position `perm_idx[j]`.
#
# ---
#
# ### How does the inverse permutation recover the original layout?
#
#     inv_perm_idx  = order(perm_idx)   # R's `order()` inverts a permutation
#     beta_recovered = beta_perm_raw[ inv_perm_idx ]
#
# This maps each element back to its original pixel position, restoring the
# spatial pattern. `beta_recovered` should be numerically identical to
# `beta_orig` (up to machine epsilon and SVD sign choices).
#
# ---
#
# ### Note on SVD sign ambiguity
#
# Singular vectors are defined only up to a sign flip (±1 per column).
# After permuting X, the solver may flip the sign of one or more eigenvectors.
# We correct for this by aligning the sign of each V_perm[:,j] with V[:,j]
# before computing scores and beta vectors, so the comparison is meaningful.

library(spatstat)
library(tidyverse)
library(tools)
library(survival)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(viridis)

# File load Functions (Except the loop, else verbatim)
getCurrentFileLocation <- function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0) this_file <- rstudioapi::getSourceEditorContext()$path
  return(dirname(this_file))
}

# Parameters (verbatim) - (From simulation_1_results.R) and (Simulation_2_result.R)
tt      <- 1 # There should be data for 30 iterations, but repo only has 2 (0 and 1)
sigma0  <- 2;  sigma1 <- 2
minmax0 <- c(-50, 50);  minmax1 <- c(-10, 90)
k0 <- 1;  k1 <- 1

#  Load PD files (verbatim)
filename_control <- list.files(
  path    = paste0(getCurrentFileLocation(), "/simulation_1_pd"),
  pattern = paste0("control_iter_", tt - 1, "_"),
  full.names = TRUE
)
filename_test <- list.files(
  path    = paste0(getCurrentFileLocation(), "/simulation_1_pd"),
  pattern = paste0("feature_iter_", tt - 1, "_"),
  full.names = TRUE
)
filename <- c(filename_control, filename_test)

pd.total <- NULL
for (ii in 1:length(filename)) {
  pd       <- read.table(filename[ii], col.names = c("dim", "birth", "death"))
  pd$group <- substr(basename(filename[ii]), 1, 1)
  pd$id    <- strtoi(strsplit(file_path_sans_ext(basename(filename[ii])), "_")[[1]][5])
  pd       <- pd %>% filter(!(death == Inf & dim == 1))
  pd[pd$death == Inf, 3] <- pd[pd$death == Inf, 2]
  pd.total <- rbind(pd.total, pd)
}

pd0.total <- pd.total %>% filter(dim == 0)
pd1.total <- pd.total %>% filter(dim == 1)

pixnum0 <- ceiling(minmax0[2]) - floor(minmax0[1])  # = 100, For Dimension 0, it can be seen that many features (b, d) occures in the middle so we choose -50, 50
pixnum1 <- ceiling(minmax1[2]) - floor(minmax1[1])  # = 100 For Dimension 1, many features (b, d) occur relateively later to dim 0, so we choose -10, 90

gid <- unique(pd.total[, 4:5]) %>% arrange(group, id)
N   <- nrow(gid)

# MDW weight function (verbatim) 
mdw <- function(pd) {
  apply(cbind(pd$death - pd$birth, abs(pd$birth), abs(pd$death)), 1, max)
}

# persistence image matrices (verbatim from simulation_1_result.R) 
pfmdw0 <- as.data.frame(matrix(NA, N, pixnum0 * (pixnum0 + 1) / 2))
pfmdw1 <- as.data.frame(matrix(NA, N, pixnum1 * (pixnum1 + 1) / 2))

# Grid from -49.5 t0 49.5 and -9.5 to 89.5
cunit0 <- (minmax0[2] - minmax0[1]) / (2 * pixnum0)
gridx0 <- seq(minmax0[1] + cunit0, minmax0[2] - cunit0, length.out = pixnum0)
gridy0 <- seq(minmax0[1] + cunit0, minmax0[2] - cunit0, length.out = pixnum0)

cunit1 <- (minmax1[2] - minmax1[1]) / (2 * pixnum1)
gridx1 <- seq(minmax1[1] + cunit1, minmax1[2] - cunit1, length.out = pixnum1)
gridy1 <- seq(minmax1[1] + cunit1, minmax1[2] - cunit1, length.out = pixnum1)



for (ii in 1:nrow(gid)) {
  pd0 <- pd0.total %>% filter(group == gid$group[ii], id == gid$id[ii])
  pd1 <- pd1.total %>% filter(group == gid$group[ii], id == gid$id[ii])

  # H0 persistence image 
  dim0.x  <- pd0$birth;  dim0.y <- pd0$death
  pd0.ppp <- ppp(dim0.x, dim0.y, minmax0, minmax0) ### spatstat.ppp object
  mdw0    <- mdw(pd0)
  pf0     <- density(pd0.ppp, sigma = sigma0, dimyx = c(pixnum0, pixnum0), weights = mdw0)
  df0     <- data.frame(x = rep(gridx0, each = pixnum0),
                        y = rep(gridy0, pixnum0),
                        z = as.vector(pf0$v))
  dfmdw0  <- df0 %>% filter(y >= x) %>% mutate(newz = pmax(z, 0))
  pfmdw0[ii, ] <- dfmdw0$newz

  # H1 Image
  dim1.x  <- pd1$birth;  dim1.y <- pd1$death
  pd1.ppp <- ppp(dim1.x, dim1.y, minmax1, minmax1)
  mdw1    <- mdw(pd1)
  pf1     <- density(pd1.ppp, sigma = sigma1, dimyx = c(pixnum1, pixnum1), weights = mdw1)
  df1     <- data.frame(x = rep(gridx1, each = pixnum1),
                        y = rep(gridy1, pixnum1),
                        z = as.vector(pf1$v))
  dfmdw1  <- df1 %>% filter(y >= x) %>% mutate(newz = pmax(z, 0))
  pfmdw1[ii, ] <- dfmdw1$newz
}

X0 <- pfmdw0 # pixel x pixel: im: describing value for each pixel for control and feature [1, ] control, [2, ] feature for 0th homology
X1 <- pfmdw1 # same same

################################################################################
# The two grids can be visualized as a raster by following: 
# df_grid <- data.frame(
#   x = rep(gridx0, each = pixnum0),
#   y = rep(gridy0, pixnum0)
# ) %>% filter(y >= x)
# df_both <- bind_rows(
#   df_grid %>% mutate(z = as.numeric(pfmdw0[1, ]), subject = "control"),
#   df_grid %>% mutate(z = as.numeric(pfmdw0[2, ]), subject = "feature")
# )
# 
# ggplot(df_both, aes(x, y, fill = z)) +
#   geom_raster() +
#   geom_abline(slope = 1, intercept = 0, colour = "grey80", linetype = "dashed") +
#   scale_fill_viridis_c() +
#   facet_wrap(~ subject) +
#   theme_tufte() +
#   labs(x = "birth", y = "death", fill = "density")
################################################################################

# Survival data simulation (verbatim)
simulWeib <- function(N, lambda, rho, beta, rateC) {
  x <- c(rep(0, N / 2), rep(1, N / 2))
  v   <- runif(n = N)
  Tlat <- (-log(v) / (lambda * exp(x * beta)))^(1 / rho)
  C <- rexp(n = N, rate = rateC)
  time  <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)
  data.frame(id = 1:N, time = time, status = status, x = x)
}
N = 2
set.seed(20)
dat <- simulWeib(N = N, lambda = 0.01, rho = 1, beta = 0.8, rateC = 0.001)
dat$age <- rpois(N, 40)
dat$sex <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))


############## Original pipeline (not verbatim but basiaally the idea is hte same)
scaleX0  <- scale(X0, center = TRUE, scale = FALSE)
svd0  <- svd(scaleX0)
Eigvec0  <- svd0$v              # p × min(N,p) right singular vectors
est_Eigval0 <- svd0$d^2           # eigenvalues of (X0-mu)'(X0-mu)
scores_orig <- as.matrix(scaleX0) %*% Eigvec0[, 1:k0]   # N × k0

comdf_orig     <- dat
comdf_orig$dim0.1 <- scores_orig[, 1]
cox_orig <- coxph(Surv(time, status) ~ dim0.1, data = comdf_orig)


#######################Permuted pipeline
# Fixed column permutation of X0 
set.seed(42)
perm_idx     <- sample(ncol(X0))
inv_perm_idx <- order(perm_idx)
X0_perm     <- X0[, perm_idx]
scaleX0_perm  <- scale(X0_perm, center = TRUE, scale = FALSE)
svd0_perm <- svd(scaleX0_perm)
Eigvec0_perm   <- svd0_perm$v
est_Eigval0_perm <- svd0_perm$d^2
scores_perm      <- as.matrix(scaleX0_perm) %*% Eigvec0_perm[, 1:k0]

# Sign correction: SVD eigenvectors are defined up to ±1; align perm with orig
for (j in seq_len(k0)) {
  if (sum(Eigvec0_perm[, j] * Eigvec0[, j]) < 0) {
    Eigvec0_perm[, j] <- -Eigvec0_perm[, j]
    scores_perm[, j]  <- -scores_perm[, j]
  }
}

comdf_perm   <- dat
comdf_perm$dim0.1 <- scores_perm[, 1]
cox_perm <- coxph(Surv(time, status) ~ dim0.1, data = comdf_perm)

#### Numerical Comparisons and Checks

tribble(
~Type, ~Original, ~Permuted,
"Eigen Value 1", est_Eigval0[1], est_Eigval0_perm[1],
"Eigen value 2", est_Eigval0[2], est_Eigval0_perm[2],
"PC", scores_orig[1, 1], scores_perm[1, 1],
"Cox Coef", coef(cox_orig)["dim0.1"], coef(cox_perm)["dim0.1"]
) %>% 
  mutate(Diff = round(Original - Permuted))


# Beta back-projection
cox_coef_orig <- coef(cox_orig)["dim0.1"]
cox_coef_perm <- coef(cox_perm)["dim0.1"]

beta_orig <- as.vector(Eigvec0[, 1:k0, drop = FALSE] %*% cox_coef_orig)
beta_perm_raw  <- as.vector(Eigvec0_perm[, 1:k0, drop = FALSE] %*% cox_coef_perm)  # scrambled
beta_recovered <- beta_perm_raw[inv_perm_idx]                                      # restored by reverse indexing

glue::glue("Floating point miscalculation between original and recorvered Beta:", max(abs(beta_orig - beta_recovered)), "\n")

# ── H0 upper-triangle grid for plotting (same ordering as columns of X0) ──────
df_grid <- data.frame(
  x = rep(gridx0, each = pixnum0),
  y = rep(gridy0, pixnum0)
) %>% filter(y >= x)

zlim <- max(abs(c(beta_orig, beta_perm_raw, beta_recovered)))

df_A <- df_grid; df_A$z <- beta_orig
df_B <- df_grid; df_B$z <- beta_perm_raw
df_C <- df_grid; df_C$z <- beta_recovered


# beta back-projection 
make_beta_panel <- function(df, title) {
  ggplot(df, aes(x, y, fill = z)) +
    geom_raster() +
    geom_abline(slope = 1, intercept = 0, colour = "grey80", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    scale_fill_gradient2(
      low    = "#2166ac",
      mid    = "white",
      high   = "#d6604d",
      limits = c(-zlim, zlim)
    ) +
    theme_tufte() +
    labs(title = title, x = "birth", y = "death", fill = "β") +
    theme(legend.position = "bottom")
}

pA <- make_beta_panel(df_A, "A: H0 Original PI weighetd with Beta Coef from Cox")
pB <- make_beta_panel(df_B, "B: H0 Sample (permuted) with β weighted")
pC <- make_beta_panel(df_C, "C: H0 Recovered (inverse permutation)")

grid.arrange(
  pA, pB, pC,
  ncol = 3,
  top  = "Beta back-projection on H0 birth-death plane (upper triangle)"
)

