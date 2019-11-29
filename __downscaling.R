library(nctools)
library(fields)
library(maptools)
library(lubridate)
library(kali)

iDir = "/data/ISIMIP/ISIMIP2a/InputData/ocean/downscaling"
iDir = "ISIMIP2a/InputData/ocean/downscaling"

domain = "biscay-celtic"
lon = c(-20, 3)
lat = c(35, 61)
gridFile = NULL
bcsdModel = "m2.0"

source("90.9-downscaling_paper-RCP_predictions.R")

domain = "humboldt-n"
lon = c(-93, -70)
lat = c(-20, 6)
gridFile = NULL
bcsdModel = "m4.2"

# source("90.9-downscaling_paper-RCP_predictions.R")

for(i in 1:5) beepr::beep(4)

