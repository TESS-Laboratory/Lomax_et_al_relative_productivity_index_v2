## Helper functions for evaluating model and raster performance metrics

#' @title Calculate performance
#' @description Calculates model performance metrics from vector of true vs.
#' predicted GPP. Returns a named vector of overall MAE, RMSE and R-squared values.
#' 
#' @usage calc_performance(truth, predicted, na.rm = FALSE)
#' 
#' @param truth numeric. A vector/data frame column of observed GPP values
#' @param predicted numeric. A vector/data frame column of model predicted GPP values
#' @param na.rm logical. Drop missing values where present

calc_performance <- function(truth, predicted, na.rm = FALSE) {
  
  error <- predicted - truth
  
  mae <- mean(abs(error), na.rm = na.rm)
  rmse <- sqrt(mean(error ^ 2, na.rm = na.rm))
  rsq <- 1 - sum((error ^ 2), na.rm = na.rm) / sum((truth - mean(truth, na.rm = na.rm)) ^ 2, na.rm = na.rm)
  
  c(mae = mae, rmse = rmse, rsq = rsq)
  
}

#' @title Calculate performance over raster
#' @description Calculates pixel-level performance metrics from rasters of true vs.
#' predicted GPP. Returns a three-layer raster with pixelwise MAE, RMSE and R-squared.
#' 
#' @usage calc_performance_raster(truth, predicted, na.rm = FALSE)
#' 
#' @param truth spatRaster of observed GPP values
#' @param predicted spatRaster of model predicted GPP values
#' @param na.rm logical. Drop missing values where present

calc_performance_raster <- function(truth, predicted, na.rm = FALSE) {
  
  annual_anomaly <- truth - mean(truth, na.rm = na.rm)
  pred_anomaly <- predicted - mean(predicted, na.rm = na.rm)
  error <- pred_anomaly - annual_anomaly
  
  mae <- mean(abs(error), na.rm = na.rm)
  rmse <- sqrt(mean(error ^ 2, na.rm = na.rm))
  rsq <- 1 - (sum(error ^ 2, na.rm = na.rm) / sum(annual_anomaly ^ 2, na.rm = na.rm))
  
  output <- c(mae, rmse, rsq)
  names(output) <- c("mae", "rmse", "rsq")
  output
  
}