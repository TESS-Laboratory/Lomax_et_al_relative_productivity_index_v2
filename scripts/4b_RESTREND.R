## This script calculates RESTREND residuals for the updated GPP, precipitation
## and temperature datasets and outputs the results as a .tif file.

## Guy Lomax
## G.Lomax@exeter.ac.uk
## 18th December 2024

## 1. Setup ----

source("scripts/load.R")

YEARS <- 2000:2018

# nc <- availableCores() / 4
# 
# plan(multisession, workers = nc)

# Load RPI data for analysis years
rpi_rast <- rast("data/processed/raster/outputs/rpi_rast_v2.tif")

gpp <- subset(rpi_rast, str_detect(names(rpi_rast), "GPP"))
gpp_years <- subset(gpp, str_detect(names(gpp), paste(YEARS, collapse = "|")))

# Load ppt and T predictors for analysis years

dynamic_covariates <- rast("data/processed/raster/dynamic_covariates.tif") %>%
  crop(rpi_rast[[1]], mask = TRUE)

precipitation <- subset(dynamic_covariates, str_detect(names(dynamic_covariates), "precipitation"))
precipitation_years <- subset(precipitation, str_detect(names(precipitation), paste(YEARS, collapse = "|")))
tMean <- subset(dynamic_covariates, str_detect(names(dynamic_covariates), "temperature_2m"))
tMean_years <- subset(tMean, str_detect(names(tMean), paste(YEARS, collapse = "|")))

# Stack variables into combined raster

combined_raster <- c(gpp_years, precipitation_years, tMean_years)

##### RESTREND function

calc_restrend_resids <- function(combined_inputs, n_layers) {
  
  if (any(is.na(combined_inputs))) {
    return(rep(NA, n_layers))
  }
  
  if (length(combined_inputs) != 3 * n_layers) {
    stop("Error: combined raster has incorrect number of layers")
  }
  
  # Divide input vector into components
  gpp_values <- combined_inputs[1:n_layers]
  ppt_values <- combined_inputs[(n_layers + 1):(2 * n_layers)]
  tMean_values <- combined_inputs[(2 * n_layers + 1):(3 * n_layers)]
  
  # Fit model to each pixel
  
  # Use .lm.fit instead of lm() - much faster
  X <- cbind(1, ppt_values, tMean_values)
  fit <- .lm.fit(X, gpp_values)
    
  resids <- fit$residuals
    
  resids
  
}

tic()
restrend_rast <- terra::app(combined_raster, calc_restrend_resids, n_layers = nlyr(gpp_years))
x <- toc()

names(restrend_rast) <- paste0("resid.", YEARS)

writeRaster(restrend_rast, "data/processed/raster/outputs/restrend.tif", overwrite = TRUE)



