library(nctools)
library(fields)
library(maptools)
library(lubridate)

domain = "biscay-celtic"
lon = c(-20, 3)
lat = c(35, 61)
gridFile = NULL

# source("00.1-extract_clean.R") # add naming
# source("00.2-data_processing.R") # derived variables
# source("01.0-downscaling_observations.R") # quantiles from observed variables
# + possible interpolation (fishmip names) + grid
# source("01.1-downscaling_data.R") # interpolation + quantiles + databases
# source("01.2-downscaling_modelData.R")
# source("01.2-downscaling_model.R") # create models for the set of variables
# source("01.4-downscaling_prediction.R") # make predictions + ncdf files
# source("01.3-downscaling_modelDev-paper.R")
# source("01.3-downscaling_modelPred-paper.R")
# source("01.3-downscaling_modelInd-paper.R")
# source("01.3-downscaling-finalModels-paper.R")
# source("90.0-downscaling_paper-selectSpatialCovariate.R")
# source("90.5-downscaling_paper-finalModels.R")
# source("90.6-downscaling_paper-finalPredictions.R")
# source("90.7-downscaling_paper-finalIndicators.R")
beepr::beep(4)

domain = "humboldt-n"
lon = c(-93, -70)
lat = c(-20, 6)
gridFile = NULL

# source("00.1-extract_clean.R") # add naming
# source("00.2-data_processing.R") # derived variables
# source("01.0-downscaling_observations.R") # quantiles from observed variables
# + possible interpolation (fishmip names) + grid
# source("01.1-downscaling_data.R") # interpolation + quantiles + databases
# source("01.2-downscaling_modelData.R")
# source("01.2-downscaling_model.R") # create models for the set of variables
# source("01.3-downscaling_prediction.R") # make predictions + ncdf files
# source("01.3-downscaling_modelDev-paper.R")
# source("01.3-downscaling_modelPred-paper.R")
# source("01.3-downscaling_modelInd-paper.R")
# source("90.0-downscaling_paper-selectSpatialCovariate.R")
# source("90.5-downscaling_paper-finalModels.R")
# source("90.6-downscaling_paper-finalPredictions.R")
# source("90.7-downscaling_paper-finalIndicators.R")



