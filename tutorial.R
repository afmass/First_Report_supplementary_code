##################################################
## Project:         Fire in Aude department (France)
##
## Script purpose:  Design-based estimation of the productive forest timber volume in
##                  areas affected by the fire using multi-source forest inventory
##
## Date:            2025-10-22
##
## Authors:         Alexander Massey (alexander.massey@ign.fr)
##                  Cedric Vega (Cedric.Vega@ign.fr)
##
## Notes:           -The fires began on June 26, 2025 and lasted until Aug 10, 2025.
##                  -Canopy height auxiliary information provided by stereo-photogrammetry (2024).
##                  -The forest mask was produced by IGN (2015)
##                  -Plot level data is from the French National Forest Inventory (2018-2024)
##################################################


# List of required packages
required_packages <- c("dplyr", "terra", "arrow", "sf", "here")

# Install any missing packages
installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) {
    message(paste("Installing missing package:", pkg))
    install.packages(pkg)
  }
}


library(dplyr)
library(terra)
library(arrow)
library(sf)
library(here)  # for relative paths to proper work directories


# Source helper functions
source(here("helper_functions.R"))




###############################
# DETECT FIRE AREA FROM dNBR+ #
###############################

# Read the raster file
NormBurnRatioPlusDiff <- rast(here("data", "NormBurnRatioPlusDiff.tif"))

# Inspect it
plot(NormBurnRatioPlusDiff)

# Create fire mask and save it in ./data folder
vectorize_burn_areas(NormBurnRatioPlusDiff, output_dir = "./data")

# Read fire mask output back into memory
fire_mask <- sf::st_read(here("data", "fire_mask.gpkg"))  # read fire mask into memory

# Inspect it
plot(fire_mask[1], col = "darkgray")   

############################################
# CALCULATE EXHAUSTIVE MEAN IN BURNED AREA #
############################################

canopy_metrics <- terra::rast(here("data", "canopy_metrics_bdf_30m.tif"))


# Rasterize fire_mask
fire_mask_vect <- vect(fire_mask) %>%
  project(canopy_metrics)  # use same coordinate system

# Intersect with canopy height model
canopy_in_fire <- mask(crop(canopy_metrics, fire_mask_vect), fire_mask_vect)


# Exhaustive means from wall-to-wall auxiliary info in fire_mask for model predictor variables
# Computed with 165645 pixels (means are negligible from theoretical true exhaustive vector)
exhaustive_mean <- as.data.frame(canopy_in_fire,
                                 xy = FALSE,
                                 na.rm = FALSE) %>%
  filter(!is.na(gap_area_percent)) %>%
  mutate(logmean = ifelse(is.na(mean), 0, log(mean)),
         sd = ifelse(is.na(sd), 0, sd),
         closed = ifelse(CODE_TFV_NUM %in% 2:26, 1,0),
         open = ifelse(CODE_TFV_NUM %in% 27:30, 1,0),
         other = ifelse(!(CODE_TFV_NUM %in% 2:30) | is.na(CODE_TFV_NUM), 1,0),
         logmean_in_closed_forest = logmean*closed,
         sd_in_closed_forest = sd*closed,
         gap_area_percent_in_closed_forest = gap_area_percent*closed,
         fire = 1
  ) %>%
  dplyr::select(closed, open, other, sd_in_closed_forest, 
                logmean_in_closed_forest, gap_area_percent_in_closed_forest, 
                fire) %>%
  colMeans %>% t %>%
  as.data.frame





############################################
# LOAD NFI AND PREPARE DATA FOR ESTIMATION #
############################################


# Get surface area of the burned area from fire mask
fire_surface_area_ha <- st_area(fire_mask) %>% 
  units::set_units(ha) %>% 
  sum %>% 
  as.numeric # 14111.6 ha in total



# French NFI data in burned area as well as expanded area to borrow strength for model fitting.
# NOTE: logmean has been set to 0 if gap area was 100% or mean canopy height was below 5 meters
# We only need the volume model in areas with volume (i.e. the closed forest) so we use the interaction variable below.
# The other areas (i.e. open forest, heathland and non-forest) are zeros or misclassifications from photo-interpretation.
nfi_plots <- read.csv(here("data", "NFI_data.csv")) %>%
  mutate(
    logmean_in_closed_forest = logmean * closed,
    gap_area_percent_in_closed_forest = gap_area_percent * closed
  )


# Subsetted datasets to within and outside fire area
nfi_plots_not_in_fire <- nfi_plots[nfi_plots$fire != 1,] # subset to plots outside fire
nfi_plots_in_fire <- nfi_plots[nfi_plots$fire == 1,] # subset to plots outside fire





##########################
## FIT ASSISTING MODELS ##
##########################

# Fitting only in the fire area problematic because very few plots in closed/open forests
table(nfi_plots_in_fire$forest_map)



# Model for REG External
# fire could be included but as fire=0 outside of fire it produces an NA and warning messages but doesn't change results)
mod_external <- lm(vol ~ open + closed + logmean_in_closed_forest + gap_area_percent_in_closed_forest, 
                   data = nfi_plots_not_in_fire, weights = w)
summary(mod_external)  



# Model for REG Internal
mod_internal <- lm(vol ~ open + closed + logmean_in_closed_forest + gap_area_percent_in_closed_forest + fire, 
                   data = nfi_plots, weights = w)
summary(mod_internal)


# The adjusted R2 for mod_internal is inflated from predicting lots of zeros so check R^2 on closed forest only
summary(lm(vol ~ open + closed + logmean_in_closed_forest + gap_area_percent_in_closed_forest + fire, 
           data = nfi_plots %>% filter(closed==1), weights = w))$adj.r.squared


# Demonstration for why indicator variable on small area works
# Update nfi_plots with residuals from the internal and external model fits
residuals_internal <- nfi_plots %>%
  mutate(residuals_internal = vol-predict(mod_internal, newdata = nfi_plots)) %>%
  select(residuals_internal, w, fire)


# Weighted mean of residuals is zero by construction!
weighted.mean(x = residuals_internal[residuals_internal$fire==1,"residuals_internal"], 
              w = residuals_internal[residuals_internal$fire==1,"w"]) %>%
  round(10)  # we don't care about rounding error





################
## ESTIMATION ##
################


################
# ONE-PHASE (NO AUXILIARY INFO)


Y1p <- one_phase_non_uniform(y = unlist(nfi_plots[nfi_plots$fire==1,"vol"]), 
                             w = unlist(nfi_plots[nfi_plots$fire==1,"w"]),
                             ci_method = "normal", type = "total", 
                             surface_area = fire_surface_area_ha)
print_estimates(Y1p)

# Same thing but only for mean density not total
one_phase_non_uniform(y = unlist(nfi_plots[nfi_plots$fire==1,"vol"]), 
                      w = unlist(nfi_plots[nfi_plots$fire==1,"w"]),
                      ci_method = "normal", type = "mean", 
                      surface_area = fire_surface_area_ha) %>% 
  print_estimates()

# Technically using "Student t" assumes normality which is a slight design-based violation
# but still recommended according to doi:10.1139/cjfr-2012-0381
one_phase_non_uniform(y = unlist(nfi_plots[nfi_plots$fire==1,"vol"]), 
                      w = unlist(nfi_plots[nfi_plots$fire==1,"w"]),
                      ci_method = "t", type = "total", 
                      surface_area = fire_surface_area_ha)$ci


# For fun use an intercept only model to reproduce results using sae_external()
mod_intercept <- lm(vol ~ 1, 
                   data = nfi_plots_not_in_fire, weights = w)
sae_external(mod_intercept, nfi_plots, nfi_plots$fire, exhaustive_mean, fire_surface_area_ha, 
             type = "total", alpha = 0.05, ci_method = "normal") %>% print_estimates()




################
# REG EXTERNAL (FIT INDEPENDENTLY OF SMALL AREA TO REDUCE VARIANCE TO RESIDUALS)


REG_external <- sae_external(ext_model_obj = mod_external,
                             data = nfi_plots,
                             small_area = nfi_plots$fire,
                             exhaustive_mean_df = exhaustive_mean,
                             surface_area = fire_surface_area_ha,
                             type = "total", alpha = 0.05, ci_method = "normal")

print_estimates(REG_external)


# Try "external model assumption" fitting everywhere and ignoring sampling variability on model coefficients
mod_external_by_assumption <- lm(vol ~ open + closed + logmean_in_closed_forest + gap_area_percent_in_closed_forest + fire, 
                                 data = nfi_plots, weights = w)
REG_external_by_assumption <- sae_external(ext_model_obj = mod_external_by_assumption,
                                           data = nfi_plots,
                                           small_area = nfi_plots$fire,
                                           exhaustive_mean_df = exhaustive_mean,
                                           surface_area = fire_surface_area_ha,
                                           type = "total", alpha = 0.05, ci_method = "normal")

print_estimates(REG_external_by_assumption)




################
# REG INTERNAL (FIT EVERYWHERE AND ACCOUNT FOR SAMPLING VARIABILITY WITH G-WEIGHT VARIANCE)
#   NOTE: Use the residual indicator technique to force mean residual to be zero in small area


REG_internal <- g_weight_internal(model_obj = mod_internal,
                                  data = nfi_plots,
                                  small_area = nfi_plots$fire,
                                  exhaustive_mean_df = exhaustive_mean,
                                  surface_area = fire_surface_area_ha,
                                  type = "total", alpha = 0.05, ci_method = "normal")

print_estimates(REG_internal)

# g_weight_internal() also calculates the variance by external model assumption which we can confirm with sae_external()
REG_external_by_assumption$variance
REG_internal$ext_variance


# We can also assess the bias avoided from including the indicator fire in the model
Y_synthetic <- g_weight_internal(model_obj = update(mod_internal, . ~ . - fire),
                                 data = nfi_plots,
                                 small_area = nfi_plots$fire,
                                 exhaustive_mean_df = exhaustive_mean,
                                 surface_area = fire_surface_area_ha,
                                 type = "total", alpha = 0.05, ci_method = "normal",
                                 allow_bias = TRUE)



# Results table for article
print_estimates(Y1p)
print_estimates(REG_external)
print_estimates(REG_internal)
print_estimates(Y_synthetic)


# The variance estimators are based on variants of the Hajek estimator of the total
# which is known to have a lower variance than the Horvitz-Thompson estimator of the total
# when the inclusion densities are negatively or weakly correlated with the target parameter.
# See p. 182 in Model Assisted Survey Sampling by Särndal, Swensson, and Wretman.
# We check this correlation below to show it is indeed negative.
cor.test(nfi_plots$vol, (1/nfi_plots$w))

