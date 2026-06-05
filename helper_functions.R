##################################################
## Project:         Fire in Aude department (France)
##
## Script purpose:  helper functions for processing fire mask and making model-assisted estimation
##
## Date:            2025-09-08
##
## Authors:         Alexander Massey (alexander.massey@ign.fr)
##                  Cedric Vega (Cedric.Vega@ign.fr)
##
## Notes:           Contains 4 functions and 1 pure helper for MODEL-ASSISTED ESTIMATION
##                    - one-phase estimation
##                    - two-phase estimation with (potentially true) external model
##                    - two-phase estimation accounting for internal fit using g-weight
##                        - based on a first order Taylor approximation 
##                        - for technical details with derivation see https://doi.org/10.3929/ethz-a-010579388
##                        - Mandallaz's cluster sampling estimators follow exact same 
##                          weighting scheme as non-uniform i.i.d. sampling
##
##                  Contains 1 function for PROCESSING FIRE MASK
##                    - vectorize_burn_areas
##
##################################################





#################################
##  MODEL-ASSISTED ESTIMATION  ##
#################################

# No auxilary information
one_phase_non_uniform <- function(y, w, alpha = 0.05, 
                                  ci_method = c("normal", "t"),
                                  type = c("mean", "total"),
                                  surface_area = NULL) {
  ci_method <- match.arg(ci_method)
  type <- match.arg(type)
  
  # --- Remove NAs consistently ---
  keep <- !is.na(y) & !is.na(w)
  y <- y[keep]
  w <- w[keep]
  
  n <- length(y)
  if (n < 2) stop("Need at least 2 non-missing observations")
  
  # --- Design-based mean ---
  mu_hat <- weighted.mean(y, w = w)
  
  # --- Average weight ---
  bar_W <- mean(w)
  
  # --- Variance estimator ---
  var_hat <- (1 / n) * (1 / (n - 1)) *
    sum(((w / bar_W)^2) * (y - mu_hat)^2)
  
  # --- Scale if total requested ---
  if (type == "total") {
    if (is.null(surface_area)) stop("surface_area must be provided for totals")
    mu_hat <- mu_hat * surface_area
    var_hat <- var_hat * (surface_area^2)
  }
  
  # --- Confidence intervals ---
  se_hat <- sqrt(var_hat)
  if (ci_method == "normal") {
    crit <- qnorm(1 - alpha / 2)
  } else {
    crit <- qt(1 - alpha / 2, df = n - 1)
  }
  ci_lower <- mu_hat - crit * se_hat
  ci_upper <- mu_hat + crit * se_hat
  
  # --- Return ---
  list(
    estimate = mu_hat,
    variance = var_hat,
    se = se_hat,
    ci = c(lower = ci_lower, upper = ci_upper),
    n = n,
    type = type,
    ci_method = ci_method,
    alpha = alpha
  )
}



# True external model
sae_external <- function(ext_model_obj, data, small_area, exhaustive_mean_df, 
                         surface_area, type = c("mean", "total"), 
                         alpha = 0.05, ci_method = c("normal", "t")) {
  type <- match.arg(type)
  ci_method <- match.arg(ci_method)
  
  # --- Extract fitted coefficients from model ---
  beta_hat <- coef(ext_model_obj)
  
  # --- Synthetic prediction part ---
  y_syn <- predict(ext_model_obj, newdata = exhaustive_mean_df)
  
  # --- Residual correction for small area ---
  mf <- model.frame(ext_model_obj, data = data)   # align to model
  y <- model.response(mf)
  X <- model.matrix(ext_model_obj, data = data)
  weights <- model.weights(mf)
  
  idx_sa <- which(small_area == 1)
  X_sa <- X[idx_sa, , drop = FALSE]
  y_sa <- y[idx_sa]
  w_sa <- weights[idx_sa]
  
  # Residuals in small area
  res_sa <- y_sa - as.numeric(X_sa %*% beta_hat)
  
  # Weighted mean of residuals
  res_bar <- weighted.mean(res_sa, w_sa)
  
  # --- Final estimate (per hectare) ---
  y_ext_hat <- y_syn + res_bar
  var_ext <- {
    n_sa <- length(idx_sa)
    W_bar <- mean(w_sa)
    (1 / n_sa) * (1 / (n_sa - 1)) *
      sum(((w_sa / W_bar)^2) * (res_sa - res_bar)^2)
  }
  
  # Scale both estimate and variance if total is requested
  # (variance of total = variance of mean * area^2)
  if (type == "total") {
    y_ext_hat <- y_ext_hat * surface_area
    var_ext <- var_ext * (surface_area^2)
  }
  
  # --- Confidence intervals ---
  se_ext <- sqrt(var_ext)
  if (ci_method == "normal") {
    crit <- qnorm(1 - alpha / 2)
  } else {
    n_sa <- length(idx_sa)
    crit <- qt(1 - alpha / 2, df = n_sa - 1)
  }
  ci_lower <- y_ext_hat - crit * se_ext
  ci_upper <- y_ext_hat + crit * se_ext
  
  list(
    coefficients = beta_hat,
    estimate = y_ext_hat,
    variance = var_ext,
    se = se_ext,
    ci = c(lower = ci_lower, upper = ci_upper),
    n_sa = length(idx_sa),
    n = length(ext_model_obj$residuals),
    type = type,
    ci_method = ci_method,
    alpha = alpha
  )
}


# REG internal
# Produce estimates and variances for exhaustive two-phase small area estimation
# NOTE: This is sometimes referred to as just a "regression estimator" because two-phase
# is sometimes only considered to be when the first phase is a sample but we
# use the nomenclature used in the r package "forestinventory" where a census is a 
# special type of sample (obviously). Anybody who actually programs these estimators in the forest inventory
# context will likely prefer the latter interpretation.

# model_obj should use formula with indicator variable on small area included
# data should contain plots in both extended area and small area (SA) of interest.
# exhaustive_mean_vector should be "true" mean vector over wall-to-wall pixels in small area (no intercept or SA indicator)
# weights are inverse inclusion probabilities conditional on uniform dense sample of photo-interpretation (i.e. no surface area of forest)
# surface_area is the known surface area of the target small area in hectares
g_weight_internal <- function(model_obj, data, small_area, 
                              exhaustive_mean_df,
                              ci_method = c("normal", "t"),
                              alpha = 0.05,
                              type = c("mean", "total"),
                              surface_area = NULL,
                              allow_bias = FALSE) {
  ci_method <- match.arg(ci_method)
  type <- match.arg(type)
  
  # --- Extract model pieces ---
  mf <- model.frame(model_obj, data = data)
  y <- model.response(mf)
  X <- model.matrix(model_obj, data = data)
  w <- model.weights(mf)
  
  # restrict to non-missing
  keep <- !is.na(y) & !is.na(w)
  y <- y[keep]; X <- X[keep, , drop = FALSE]; w <- w[keep]
  small_area <- small_area[keep]
  
  # --- Global regression fit ---
  n_total <- nrow(X)
  As2inv <- solve( (t(X) %*% (w * X)) / n_total )
  theta_hat <- drop( ((w * y) %*% X / n_total) %*% As2inv )  # regression coefficients
  
  # fitted and residuals
  y_hat <- as.vector(X %*% theta_hat)
  res <- y - y_hat
  
  # covariance matrix of coefficients
  WR_sq <- (w^2) * (res^2)  # this is M^2(x) * R_hat(x)^2 in Eq. 11 of Hill et al. 2018 (DOI: 10.3390/rs10071052)
  middle_term <- (t(X) %*% (WR_sq * X)) / (n_total^2)
  cov_theta <- As2inv %*% middle_term %*% As2inv
  
  # --- Small area subset ---
  idx_sa <- which(small_area == 1)
  X_sa <- X[idx_sa, , drop = FALSE]
  y_sa <- y[idx_sa]
  w_sa <- w[idx_sa]
  n_sa <- length(idx_sa)
  
  res_sa <- res[idx_sa]
  res_bar_sa <- weighted.mean(res_sa, w_sa) # should be zero by construction
  
  if(!allow_bias){
    if (round(res_bar_sa,10)!=0) stop("Mean residual not zero, did you forget the small area indicator variable?")
  }
  
  # --- Synthetic part (using exhaustive mean vector) ---
  # extract RHS of formula (predictors only)
  rhs_formula <- reformulate(attr(terms(model_obj), "term.labels"))
  # build design row from exhaustive means based on hs_formula (response is unknown...)
  Z_bar <- as.numeric(model.matrix(rhs_formula, data = exhaustive_mean_df))
  names(Z_bar) <- colnames(X)
  
  syn_part <- sum(Z_bar * theta_hat)
  
  
  # --- Final estimate ---
  est <- syn_part
  if (type == "total") {
    if (is.null(surface_area)) stop("surface_area must be provided for type = 'total'")
    est <- est * surface_area
  }
  
  # --- Variance parts ---
  g_var <- drop(t(Z_bar) %*% cov_theta %*% Z_bar)
  ext_var <- (1/n_sa) * (1/(n_sa - 1)) *
    sum(((w_sa / mean(w_sa))^2) * (res_sa - res_bar_sa)^2)
  
  
  if (type == "total") {
    g_var <- g_var * (surface_area^2)
    ext_var <- ext_var * (surface_area^2)
  }
  
  se <- sqrt(g_var)
  
  # --- CI ---
  if (ci_method == "normal") {
    z <- qnorm(1 - alpha/2)
    ci <- est + c(-1, 1) * z * se
  } else {
    df <- n_sa - 1
    tval <- qt(1 - alpha/2, df = df)
    ci <- est + c(-1, 1) * tval * se
  }
  
  # --- Output ---
  list(
    coefficients = theta_hat,
    estimate = est,
    variance = g_var,
    ext_variance = ext_var,
    se = se,
    ci = ci,
    n_sa = n_sa,
    n = n_total,
    type = type,
    ci_method = ci_method
  )
}




# Make output cleaner for viewing
print_estimates <- function(res, alpha = 0.05) {
  # Ensure necessary fields exist
  est <- res$estimate %||% res$mean   # works for both functions
  var <- res$variance
  se  <- res$se
  n_sa   <- res$n_sa %||% NA
  n   <- res$n %||% NA
  type <- res$type %||% "mean"
  ci   <- res$ci
  
  # Build data.frame
  out <- data.frame(
    Type = type,
    Estimate = est,
    Variance = var,
    SE = se,
    n_sa = n_sa,
    n = n,
    `Lower CI` = ci[1],
    `Upper CI` = ci[2],
    row.names = NULL
  )
  
  # Nice printing
  print(out, digits = 4)
  invisible(out)
}



# helper: like purrr::%||% (use fallback if NULL)
`%||%` <- function(x, y) if (!is.null(x)) x else y






############################
##  PROCESSING FIRE MASK  ##
############################

# ==============================================================================
# Function for detecting fire area from dNBR+
# ==============================================================================
#' Function for extracting a fire mask from a normalize burn ratio difference
#' Input  tif file of dNBR+ (see main article for definition)
#' Output the fire mask
vectorize_burn_areas <- function(NormBurnRatioPlusDiff,
                                 min_nbr_thr = -0.05,
                                 min_patch_size_ha = 10,
                                 buffer_dist = 100,
                                 min_group_area_ha = 50,
                                 min_nbr_value = -0.4,
                                 output_dir = "./output",
                                 crs_target = 2154,
                                 plot_results = FALSE) {
  
  # Verify inputs
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' required")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' required")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required")
  }
  
  # Create output directory if necessary
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  
  # ===========================================================================
  # 1. Thresholding and filtering
  # ===========================================================================
  message("Initial mask : threshold on the Normalize Burn Ration Plus difference...")
  NormBurnRatioPlusDiff_th <- NormBurnRatioPlusDiff < min_nbr_thr
  if (plot_results) plot(NormBurnRatioPlusDiff_th, main = "Normalize Burn Ration Plus difference thresholded")
  
  # Filter isolated pixels (morphological opening)
  struct <- matrix(1, nrow = 5, ncol = 5)
  NormBurnRatioPlusDiff_fil <- terra::focal(NormBurnRatioPlusDiff_th, w = struct, fun = min, na.policy = "omit")
  NormBurnRatioPlusDiff_fil <- terra::focal(NormBurnRatioPlusDiff_fil, w = struct, fun = max, na.policy = "omit")
  if (plot_results) plot(NormBurnRatioPlusDiff_fil, main = "Normalize Burn Ration Plus difference filtered")
  
  # ===========================================================================
  # 2. Detect and filter patched
  # ===========================================================================
  message("Patch detection...")
  NormBurnRatioPlusDiff_fil <- terra::patches(NormBurnRatioPlusDiff_fil, zeroAsNA = TRUE, directions = 8)
  if (plot_results) plot(NormBurnRatioPlusDiff_fil, main = "Patches detected")
  
  # Filter patches according to patch size
  message(paste("Filter patches <", min_patch_size_ha, "ha..."))
  NormBurnRatioPlusDiff_fil <- terra::zonal(terra::cellSize(NormBurnRatioPlusDiff_fil, unit = "ha"),
                                            NormBurnRatioPlusDiff_fil, sum, as.raster = TRUE)
  NormBurnRatioPlusDiff_fil <- terra::ifel(NormBurnRatioPlusDiff_fil < min_patch_size_ha, NA, NormBurnRatioPlusDiff_fil)
  if (plot_results) plot(NormBurnRatioPlusDiff_fil, main = "Patches filtered")
  
  # ===========================================================================
  # 3. Convert to polygons
  # ===========================================================================
  message("Convert to polygons...")
  NormBurnRatioPlusDiff_poly <- terra::as.polygons(NormBurnRatioPlusDiff_fil,
                                                   dissolve = TRUE,
                                                   values = TRUE,
                                                   na.rm = TRUE)
  
  # ===========================================================================
  # 4. Aggregate polygons in buffer distance
  # ===========================================================================
  message("Buffer and aggregate...")
  buffers <- terra::buffer(NormBurnRatioPlusDiff_poly, width = buffer_dist)
  r <- terra::rasterize(buffers, NormBurnRatioPlusDiff_fil, field = 1)
  
  # Identify connected patches
  clumps <- terra::patches(r, directions = 8)
  groups <- terra::as.polygons(clumps, dissolve = TRUE, values = TRUE)
  names(groups) <- "poly_Id"
  
  # ===========================================================================
  # 5. Surface filter
  # ===========================================================================
  message(paste("Filter polygons <", min_group_area_ha, "ha..."))
  groups$area <- terra::expanse(groups, unit = "ha")
  groups <- groups[groups$area > min_group_area_ha, ]
  groups$area <- NULL
  if (plot_results) plot(groups, main = "Polygons filtered by surface")
  
  # ===========================================================================
  # 6. Intensity filter
  # ===========================================================================
  message(paste("Filter polygons with NBR min <", min_nbr_value, "..."))
  min_vals <- terra::extract(NormBurnRatioPlusDiff, groups, fun = min, na.rm = TRUE)
  groups$min <- min_vals[, 2]
  groups <- groups[groups$min < min_nbr_value, ]
  if (plot_results) plot(groups, main = "Polygons filtered by burn ratio intensity")
  
  # ===========================================================================
  # 7. Intersect and aggregate
  # ===========================================================================
  message("Intersect and aggregate...")
  burn_mask <- terra::intersect(NormBurnRatioPlusDiff_poly, groups)
  burn_mask <- terra::aggregate(burn_mask, by = "poly_Id")
  burn_mask$area <- terra::expanse(burn_mask, unit = "ha")
  
  # ===========================================================================
  # 8. Products : Fire mask
  # ===========================================================================
  output_file <- file.path(output_dir, "fire_mask.gpkg")
  message(paste("Export:", output_file))
  terra::writeVector(burn_mask, output_file, filetype = "GPKG", overwrite = TRUE)
  
  
  # ===========================================================================
  # Return outputs
  # ===========================================================================
  message("Processing ended!")
}

