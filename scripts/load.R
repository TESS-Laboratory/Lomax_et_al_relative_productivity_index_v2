load <- function() {
  # load all the packages you need for the project
  message("Loading packages...")
  suppressPackageStartupMessages(
    source("scripts/packages.R")
  )
  
  # read the functions from R directory
  message("Loading functions...")
  function_files <- list.files("scripts/functions", pattern = ".R$", full.names = TRUE)
  walk(function_files, source)
  message("Complete")
  
  # Generate folder structure
  dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/processed/csv", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/processed/raster/merit", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/processed/raster/dist_to_river", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/processed/raster/merit", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/processed/raster/rpi", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/processed/rds", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/processed/tmp", showWarnings = FALSE, recursive = TRUE)
  
  # set global options
  terra::terraOptions(
    memfrac = 0.8,
    memmax = NA,
    tempdir = "data/processed/tmp"
  )
  
  options(scipen = 999)
  
}

load()
