##################################################
## Project:         Fire in Aude department (France)
##
## Script purpose:  Make a plot of the approximate relative efficiency of two-phase versus one-phase
##                  estimation as a function of R squared of the assisting model
##
## Date:            2025-10-22
##
## Author:          Alexander Massey (alexander.massey@ign.fr)
##
## Notes:           It is preferable to consider the R^2 of a model restricted to the areas 
##                  where the target variable (e.g., timber volume) is concentrated, i.e., the forested part.
##                  Recall that model-assisted estimation reduces the variance of Y(x) to the variance of R(x). 
##                  Non-forest plots are often included out of necessity because the full forest boundary isn't known.
##                  However, these zero values are effectively reduced to (near-) zero residuals and 
##                  therefore contribute little to the overall variance reduction.
##
## Source:          page 83, Sampling Techniques For Forest Inventories by Daniel Mandallaz                 
##################################################


# Define R² values
r2 <- seq(0, 1, length.out = 500)

# Compute percent reduction in CI length
pct_red <- 100 * (1 - sqrt(1 - r2))

# Plot
plot(r2, pct_red, type = "l", lwd = 2, col = "darkgray",
     xlab = expression(R^2),
     ylab = "Percent reduction in CI length (%)",
     main = ("Relationship between R2 and expected reduction \n in confidence interval length \n using the model-assisted estimator"),
     xlim = c(0, 1), ylim = c(0, 100))
grid()

# Add illustrative points
r2_points <- c(0.25, 0.5, 0.61, 0.75, 0.85)
pct_points <- 100 * (1 - sqrt(1 - r2_points))
points(r2_points, pct_points, pch = 19, col = "black")

# Build an expression vector (not a list) of stacked labels:
labels_list <- lapply(seq_along(r2_points), function(i) {
  bquote(R^2 == .(r2_points[i]))
})
labels_expr <- as.expression(labels_list)  # convert to expression vector

# Define vertical offsets
y_offset_top <- 9    # for R^2 above the point
y_offset_bottom <- 4 # for percent below the point

# Place labels (R^2 above, percent below)
text(r2_points - 0.04, pct_points + y_offset_top,
     labels = labels_expr,
     cex = 0.8, col = "black")


# Add second line (percent)
text(r2_points - 0.04, pct_points + y_offset_bottom,
     labels = sprintf("%.1f%%", pct_points),
     cex = 0.8, col = "black")
