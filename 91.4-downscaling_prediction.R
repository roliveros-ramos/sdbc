library(readr)
library(scam)
library(nctools)

source("99.9-auxiliar_functions.R")

# Read data ---------------------------------------------------------------

model = "gfdl-esm2m" # "*" for all models

variables =  c("npp_zint", "so_zs", "to_zb", "to_zs")
variables =  c("so_zs", "to_zb", "to_zs")
variables =  "to_zb"
useLog = c("to_zs"=FALSE, "to_zb"=FALSE, "so_zs"=FALSE, "npp_zint"=TRUE)
logConst = 1e-6

# variables = c("to_zs", "to_zb", "npp_zint", "so_zs")

fmlas = list()

scenarios = c("rcp26", "rcp45", "rcp60", "rcp85")
# scenarios = "rcp26"

inputDir = file.path("ISIMIP2a/InputData/ocean/downscaling", domain, "_input")
outputDir = file.path("ISIMIP2a/InputData/ocean/downscaling", domain)


for(ivar in variables) {

  kali::DateStamp("Downscaling for variable", ivar, "in progress...")
  xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
  xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
  xfiles = grep(x=xfiles, patt="quantile", value = TRUE)
  
  file = grep(x=xfiles, patt=model, value = TRUE)
  
  xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
  xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
  xfiles = grep(x=xfiles, patt="quantile", value = TRUE, invert=TRUE)
  pfiles = grep(x=xfiles, patt="2100", value = TRUE)
  hfiles = grep(x=xfiles, patt="2100", value = TRUE, invert=TRUE)
  
  pfile = grep(x=pfiles, patt=model, value = TRUE)
  hfile = grep(x=hfiles, patt=model, value = TRUE)
  
  pfiles = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
  pfiles = grep(x=pfiles, patt=ivar, value = TRUE)
  pfiles = grep(x=pfiles, patt=model, value = TRUE)
  pfiles = grep(x=pfiles, patt="quantile", value = TRUE, invert=TRUE)
  
  varid = strsplit(basename(file), split="_")[[1]][3]
  nenv = read_csv(file.path(inputDir, file), col_types = "dddddddi")
  penv = read_csv(file.path(inputDir, pfile), col_types = "dddddidddd")
  henv = read_csv(file.path(inputDir, hfile), col_types = "dddddid")
  
  slog = function(x) sign(x)*log(abs(x))
  
  nenv$ldc = slog(nenv$dc)
  nenv$lds = slog(nenv$ds)
  penv$ldc = slog(penv$dc)
  penv$lds = slog(penv$ds)
  henv$ldc = slog(henv$dc)
  henv$lds = slog(henv$ds)
  
  ovar = strsplit(ivar, split="_")[[1]][1]
  vars = names(nenv)
  modvar = grep(x=vars, pattern="_mod", value=TRUE)
  obsvar = gsub(x=modvar, pattern="_mod", replacement="")
  aux = setdiff(vars, c(modvar, obsvar, "qua"))
  
  varCheck = obsvar %in% vars
  
  nenv = nenv[!(nenv$lon>-5.4 & nenv$lat<43.2), ] # remove mediterranean sea
  
  if(obsvar=="sst" & mean(nenv[[modvar]], na.rm=TRUE)>200) 
    nenv[, modvar] = nenv[, modvar] - 273.15
  
  if(isTRUE(useLog[ivar])) {
    nenv[, modvar] = nenv[, modvar] + logConst
    nenv[, obsvar] = nenv[, obsvar] + logConst
  }
  
  # refit model
  
  nullFmla0 = sprintf("%s ~ te(by=%s, lon, lat, k=20) + te(lon, lat, k=15) + s(ldc, bs='cr', k=15)",
                     obsvar, modvar)
  nullFmla1 = sprintf("log(%s) ~ te(by=log(%s), lon, lat, k=20) + te(lon, lat, k=15) + s(ldc, bs='cr', k=15)",
                      obsvar, modvar)
  nullFmla = if(isTRUE(useLog[ivar])) nullFmla1 else nullFmla0
    
  kali::DateStamp("\tFitting model...") # automatize extraction of best formula? compare both covariates?
  fmla = if(is.null(fmlas[[ivar]])) nullFmla else fmlas[[ivar]]
  fmla = formula(fmla)
  final = bam(fmla, data=nenv)
  kali::DateStamp()
  
  
  # make predictions --------------------------------------------------------
  
  # future
  
  penv[, modvar] = NA_real_
  
  lon = sort(unique(penv$lon))
  lat = sort(unique(penv$lat))
  time = sort(as.numeric(as.character(unique(penv$time))))
  
  for(iScen in scenarios) {
    
    kali::DateStamp("\tPredicting scenarios", iScen)
    penv[, modvar] = penv[[iScen]]
    # remove mediterranean
    penv[(penv$lon>-5.4 & penv$lat<43.2), modvar] = NA
    
    if(isTRUE(useLog[ivar])) {
      penv[, modvar] = penv[, modvar] + logConst
    }
    
    pred = predict(final, newdata=penv, block.size=200000, type="response")
    
    dim(pred) = c(length(lon), length(lat), length(time))
    
    fileName = grep(pfiles, pattern=iScen, value=TRUE)
    
    kali::DateStamp("\n\tWriting file", fileName)
    
    if(isTRUE(useLog[ivar])) {
      
      lpred = exp(pred)-logConst
      lpred[lpred<0] = 0
      write_ncdf(lpred, filename = file.path(outputDir, fileName),
                 varid = ovar, dim = list(lon=lon, lat=lat, time=time),
                 dim.units=c("degrees_east", "degrees_north", "months"),
                 dim.longname = c("longitude", "latitude", "time"), unlim = "time")
      
    } else {
      
      write_ncdf(pred, filename = file.path(outputDir, fileName),
                 varid = ovar, dim = list(lon=lon, lat=lat, time=time),
                 dim.units=c("degrees_east", "degrees_north", "months"),
                 dim.longname = c("longitude", "latitude", "time"), unlim = "time")
      
    }
    
  }
  
  # historical
  
  lon = sort(unique(henv$lon))
  lat = sort(unique(henv$lat))
  time = sort(as.numeric(as.character(unique(henv$time))))
  
  iScen = "historical"
  kali::DateStamp("\tPredicting scenarios", iScen)
  henv[, modvar] = henv[[iScen]]
  # remove mediterranean
  henv[(henv$lon>-5.4 & henv$lat<43.2), modvar] = NA

  if(isTRUE(useLog[ivar])) {
    henv[, modvar] = henv[, modvar] + logConst
  }
  
  pred = predict(final, newdata=henv, block.size=200000, type="response")
  
  dim(pred) = c(length(lon), length(lat), length(time))
  
  fileName = grep(pfiles, pattern=iScen, value=TRUE)
  
  kali::DateStamp("\n\tWriting file", fileName)

  if(isTRUE(useLog[ivar])) {

    lpred = exp(pred)-logConst
    lpred[lpred<0] = 0    
    write_ncdf(lpred, filename = file.path(outputDir, fileName),
               varid = ovar, dim = list(lon=lon, lat=lat, time=time),
               dim.units=c("degrees_east", "degrees_north", "months"),
               dim.longname = c("longitude", "latitude", "time"), unlim = "time")
    
  } else {
    
    write_ncdf(pred, filename = file.path(outputDir, fileName),
               varid = ovar, dim = list(lon=lon, lat=lat, time=time),
               dim.units=c("degrees_east", "degrees_north", "months"),
               dim.longname = c("longitude", "latitude", "time"), unlim = "time")
    
  }
  
}

