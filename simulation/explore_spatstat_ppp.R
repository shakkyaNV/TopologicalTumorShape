# ═══════════════════════════════════════════════════════════════════════════════
# EXPLORATION NOTES — spatstat ppp and density.ppp in the persistence image pipeline
# ═══════════════════════════════════════════════════════════════════════════════
# ### 1. PPP object
# 
# ppp(x, y, window) creates a planar point pattern:
#   - x, y    : coordinates of each point (here: birth, death values from PD)
# - window  : the observation domain - a rectangle defined by [xmin,xmax] × [ymin,ymax]
# (here: minmax × minmax, e.g. [-10,90] × [-10,90] for H1)
# - weights : optional per-point weights passed to density() separately
# 
# The window is important - it defines the domain over which the KDE is evaluated.
# Points outside the window are rejected. The window is also what sets the
# axis ranges of the output image.
# 
# ###  2. KDE: density.ppp (not stats:::density)
# 
# density(ppp_obj, sigma, dimyx, weights) computes a weighted KDE:
#   
#   f̂(u,v) = Σᵢ  wᵢ · K_σ(u − xᵢ, v − yᵢ)
# 
# where K_σ is an isotropic 2D Gaussian:
#   
#   K_σ(Δx,Δy) = 1/(2πσ²) · exp( −(Δx² + Δy²) / 2σ² )
# 
# Each persistence diagram point places a Gaussian bump on the birth-death plane,
# centred at (birth, death), height scaled by its weight wᵢ (MDW here).
# The total density at any pixel = sum of all bump contributions.
# 
# Parameters:
#   sigma  - standard deviation of the Gaussian kernel. This IS the smoothing
# parameter - they are the same quantity, not two separate knobs.
# Larger σ → wider bumps → smoother image (less spatial detail).
# Smaller σ → narrower bumps → spikier image (more detail, noisier).
# There is no way to set kernel width and smoothing separately in KDE -
#   they are definitionally identical.
# 
# dimyx  - c(nrow, ncol): resolution of the output pixel grid.
# Independent of sigma. Finer dimyx = more pixels, not more smoothing.
# 
# Output: an im (image) object. Key slot: $v - a pixnum×pixnum matrix of KDE values.
# 
# ### 3. Effect of Smoothing (Note how sigma here acts as both as smoothing and sigma in Gaussian kernel)
# 
# The 2×2 grid (σ = 0.5, 2, 5, 10) visualizes the smoothing tradeoff:
#   σ = 0.5 : very narrow bumps - individual PD points visible as isolated peaks
# σ = 2   : used in the paper - moderate smoothing, structure preserved
# σ = 5   : heavy smoothing - spatial detail lost, broad region highlighted
# σ = 10  : extreme smoothing - nearly uniform blob
# 
# The σ used in Moon et al. (1.8 for lung, 2 for simulations) is chosen by
# cross-validation (see brain_cv_*.R, lung_cv_*.R) to maximise predictive
# performance of the Cox model - not by any geometric criterion.
# 
# ### 4. The Diagonal Bleeding
# 
# Persistence diagram points near the diagonal (birth ≈ death, low persistence)
# produce Gaussian bumps that spill slightly BELOW the diagonal (death < birth),
# the physically impossible region. The width of this bleed ∝ σ.
# 
# Two steps fix this:
#   (a) filter(y >= x)   - removes all pixels where death < birth entirely.
# These rows are ABSENT from df_upper (not set to 0).
# In ggplot, absent pixels show as blank background,
# not as viridis minimum colour.
# (b) pmax(value, 0)   - clips the small negative values that spatstat's
#                            FFT-based kernel evaluator can produce near the
#                            boundary. These are numerical artifacts, not real density.
# 
# ### 5. $v matrix layout and coordinate correspondence 
# 
#   $v[i, j] = KDE value at  x = xcol[j]  (birth),  y = yrow[i]  (death)
# 
# Both xcol and yrow are in INCREASING order:
#   xcol[1]   = min birth  →  col 1  = LEFT  edge
#   xcol[100] = max birth  →  col 100 = RIGHT edge
#   yrow[1]   = min death  →  row 1  = BOTTOM edge
#   yrow[100] = max death  →  row 100 = TOP  edge
# 
#   as.data.frame(dens_obj) unpacks $v into (x, y, value) triplets:
#     $v[i, j]  ←→  df row where  x == xcol[j]  AND  y == yrow[i]
# 
#   So the entire first ROW of $v ($v[1, ]):
#     = all df rows where y == yrow[1]  (y fixed at min death, x varies)
#     = a horizontal slice along the bottom of the birth-death plane
# 
# ### 6. Flattening order: arrange(x, y) 
# 
#   arrange(x, y) → x is outer/slow, y is inner/fast  →  column-major traversal
#   (x = birth = column index;  for each column, walk all rows bottom → top)
# 
#   Traversal sequence in $v terms:
#     $v[1,1], $v[2,1], ..., $v[100,1],   ← col 1 (min birth), rows 1→100
#     $v[1,2], $v[2,2], ..., $v[100,2],   ← col 2, rows 1→100
#     ...
#     $v[1,100], ..., $v[100,100]          ← col 100 (max birth), rows 1→100
# 
#   After filter(y >= x), lower-triangle pixels are absent, so each column's
# run is truncated - col j only contributes rows where yrow[i] >= xcol[j].
# 
# The ordering must be IDENTICAL when: (see [[Permutation Invariance in SVD]])
# (a) building each subject's row in matrix X  (flattening the persistence image)
#     (b) reshaping beta.est back to the grid      (reconstructing β̂(b,d) for plotting)
#   Any inconsistency would spatially scramble the coefficient function heatmap
#   (but would NOT affect Cox model estimates, which only see scalar PC scores).
# ═══════════════════════════════════════════════════════════════════════════════

library(spatstat)
library(tidyverse)
library(tools)
library(ggthemes)
library(gridExtra)
library(viridis)
library(glue)
f <- glue::glue


# $$ 1: Load persistent Diagram (verbatim from Moon)
getCurrentFileLocation <- function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0) this_file <- rstudioapi::getSourceEditorContext()$path
  return(dirname(this_file))
}

pd_path <- file.path(getCurrentFileLocation(), "simulation_1_pd",
                     "feature_iter_0_num_0_pd.txt")

pd <- read.table(pd_path, col.names = c("dim", "birth", "death"))
pd <- pd %>%
  filter(!(death == Inf & dim == 1)) %>%
  mutate(death = ifelse(death == Inf, birth, death))

pd0 <- pd %>% filter(dim == 0)
pd1 <- pd %>% filter(dim == 1)

cat("H0 points:", nrow(pd0), "\n")
cat("H1 points:", nrow(pd1), "\n")


# SS 2. Raw persistence diagram
#    (birth on x, death on y; points above the diagonal are "alive")


plot_pd <- function(pd_df, dim_label, colour) {
  diag_lim <- range(c(pd_df$birth, pd_df$death))
  ggplot(pd_df, aes(x = birth, y = death)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(colour = colour, alpha = 0.7, size = 2) +
    coord_equal() +
    labs(title = paste("Persistence diagram —", dim_label),
         x = "Birth", y = "Death") +
    theme_tufte(base_size = 12)
}

if (nrow(pd0) > 0) p_pd0 <- plot_pd(pd0, "H0 (components)", "#2166ac")
if (nrow(pd1) > 0) p_pd1 <- plot_pd(pd1, "H1 (loops)",      "#d6604d")


# SS 3. Build ppp objects

minmax0 <- c(-50, 50)
minmax1 <- c(-10, 90)
sigma <- 2
pixnum  <- 100

mdw <- function(pd_df) {
  apply(cbind(pd_df$death - pd_df$birth, abs(pd_df$birth), abs(pd_df$death)), 1, max)
}


build_ppp_and_density <- function(pd_df, minmax, pixnum, sigma) {
  ppp_obj <- ppp(pd_df$birth, pd_df$death, minmax, minmax)
  w <- mdw(pd_df)
  dens   <- density(ppp_obj, sigma = sigma, dimyx = c(pixnum, pixnum), weights = w)
  list(ppp = ppp_obj, dens = dens, weights = w)
}

out1 <- build_ppp_and_density(pd1, minmax1, pixnum, sigma) # Makse sure to Only build ppp for dimensions that have data


# SS 4: patstat ppp object
#    Visualizes raw point patterns (make sure it's not stats::density)


plot_ppp_as_ggplot <- function(ppp_obj, weights, minmax, title) {
  df <- data.frame(x = ppp_obj$x, y = ppp_obj$y, w = weights)
  ggplot(df, aes(x, y, size = w, colour = w)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_viridis_c(option = "C") +
    scale_size_continuous(range = c(1, 5)) +
    coord_equal(xlim = minmax, ylim = minmax) +
    labs(title = title, x = "Birth", y = "Death",
         colour = "MDW", size = "MDW") +
    theme_tufte(base_size = 12)
}

p_ppp1 <- plot_ppp_as_ggplot(out1$ppp, out1$weights, minmax1,
                              "ppp object — H1 (point size = MDW weight)")


# SS 5  KDE density matrix (pf$v)
#    This is the raw 2D output before any filtering


im_to_df <- function(dens_obj) {
  # spatstat im $v[i,j]: row i = yrow[i] (yrow increasing, so row 1 = min death = bottom)
  #                      col j = xcol[j] (xcol increasing, so col 1 = min birth = left)
  # $v[1,1]=bottom-left  $v[1,100]=bottom-right
  # $v[100,1]=top-left   $v[100,100]=top-right
  as.data.frame(dens_obj) # unpacks all 100x100 cells to (x, y, value) triplets
}

df_dens1 <- im_to_df(out1$dens)

p_kde_full <- ggplot(df_dens1, aes(x, y, fill = value)) +
  geom_raster() +
  geom_abline(slope = 1, intercept = 0, colour = "white", linetype = "dashed") +
  scale_fill_viridis_c(option = "D") +
  coord_equal() +
  labs(title = "Full KDE density — pf$v (2D matrix, 100×100)",
       x = "Birth", y = "Death", fill = "Density") +
  theme_tufte(base_size = 12)


# SS 6: CLip Negatives: Upper triangle only (y >= x): (Essentially feature matrix X0)


df_upper1 <- df_dens1 %>%
  filter(y >= x) %>%
  mutate(value = pmax(value, 0))

p_kde_upper <- ggplot(df_upper1, aes(x, y, fill = value)) +
  geom_raster() +
  geom_abline(slope = 1, intercept = 0, colour = "white", linetype = "dashed") +
  scale_fill_viridis_c(option = "D") +
  coord_equal() +
  labs(title = "Upper triangle only (y ≥ x) — input to feature matrix X",
       x = "Birth", y = "Death", fill = "Density") +
  theme_tufte(base_size = 12)


# SS 7 - Effect of sigma on the KDE - Bandwidth and smoothness


sigma_vals <- c(0.5, 2, 5, 10)

sigma_plots <- lapply(sigma_vals, function(s) {
  d   <- density(out1$ppp, sigma = s, dimyx = c(pixnum, pixnum), weights = out1$weights)
  df  <- as.data.frame(d) %>% filter(y >= x) %>% mutate(value = pmax(value, 0))
  ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    geom_abline(slope = 1, intercept = 0, colour = "white", linetype = "dashed") +
    scale_fill_viridis_c(option = "D", guide = "none") +
    coord_equal() +
    labs(title = paste0("σ = ", s), x = "Birth", y = "Death") +
    theme_tufte(base_size = 10)
})


# SS 8 - The flattened vector (1D representation of the upper triangle)


flat_vec <- df_upper1 %>% arrange(x, y) %>% pull(value)  # one consistent ordering (see explore_permutation_invariance.R for permuted matrices)
p_flat <- ggplot(data.frame(index = seq_along(flat_vec), z = flat_vec),
                 aes(x = index, y = z)) +
  geom_line(colour = "#2166ac", linewidth = 0.4) +
  labs(title = paste0("Flattened upper-triangle vector (length = ",
                      length(flat_vec), " = ", pixnum, "×", pixnum+1, "/2)"),
       x = "Pixel index", y = "KDE value") +
  theme_tufte(base_size = 12)



### Render all plots
if (nrow(pd0) > 0) print(p_pd0)
if (nrow(pd1) > 0) print(p_pd1)
print(p_ppp1)
print(p_kde_full)
print(p_kde_upper)

# Sigma comparison grid
grid.arrange(grobs = sigma_plots, ncol = 2,
             top = "Effect of KDE bandwidth ($\sigma$) on persistence image — H1, upper triangle")

# Flattened vector
print(p_flat)


# SS (TEST)

f("PPP Object:         {out1$ppp} \n")
f("PPP Density Object: {out1$dens} \n")
f("SVD - Dimensions:   {dim(out1$dens$v)} \n")
f("SVD - Class:        {class(out1$dens$v} \n")
f("Flattened VEctor: grid size: {pixnum*pixnum} \n")
f("Flattened Vector: Grid size after clipping: {pixnum * (pixnum+1)/2}, (includes diagonal pixels) \n")
