library(readr)
library(scam)
library(foreach)
library(doSNOW)

source("99.9-auxiliar_functions.R")

model = "ipsl" # "*" for all models
ivar = "to_zs"
modvar = "sst_mod"

scenarios = c("rcp26", "rcp45", "rcp60", "rcp85")
# scenarios = "rcp26"

modelDir = "output/final"
outputDir = file.path(iDir, domain)

block.size = 200000

# Reading files -----------------------------------------------------------

inputDir = file.path(iDir, domain, "_input")

if(!dir.exists(outputDir)) dir.create(outputDir)

# quantile files
xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
xfiles = grep(x=xfiles, patt="quantile", value = TRUE)
file = grep(x=xfiles, patt=model, value = TRUE)
# future file
xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
pfiles = grep(x=xfiles, patt="2100", value = TRUE)
pfile = grep(x=pfiles, patt=model, value = TRUE)
# historical file
xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
xfiles = grep(x=xfiles, patt="1950", value = TRUE)
hfile = grep(x=xfiles, patt=model, value = TRUE)
# read files
varid = strsplit(basename(file), split="_")[[1]][3]
nenv = read_csv(file.path(inputDir, file), col_types = "dddddddd")
penv = read_csv(file.path(inputDir, pfile), col_types = "dddddddddd")
henv = read_csv(file.path(inputDir, hfile), col_types = "ddddddd")
# process files
slog = function(x) sign(x)*log(abs(x))
nenv$ldc = slog(nenv$dc)
nenv$lds = slog(nenv$ds)
penv$ldc = slog(penv$dc)
penv$lds = slog(penv$ds)
henv$ldc = slog(henv$dc)
henv$lds = slog(henv$ds)
names(henv)[7] = modvar
penv[, modvar] = NA_real_

nenv = nenv[!(nenv$lon>-5.4 & nenv$lat<43.2), ] # remove mediterranean sea
# penv[(penv$lon>-5.4 & penv$lat<43.2), modvar] = NA
henv[(henv$lon>-5.4 & henv$lat<43.2), modvar] = NA

if(mean(nenv$sst_mod, na.rm=TRUE)>200) nenv$sst_mod = nenv$sst_mod - 273.15

# Model predictions -------------------------------------------------------

if(!dir.exists(outputDir)) dir.create(outputDir)

modelPath = dir(path=modelDir, pattern=bcsdModel)
modelPath = file.path(modelDir, grep(x=modelPath, pattern=domain, value=TRUE))

mod = readRDS(modelPath)


outputFiles = dir(path=inputDir, pattern="nc4$")
outputFiles = grep(x=outputFiles, pattern=ivar, value=TRUE)

# historical data
iScen = "historical"
kali::DateStamp("\tPredicting ", iScen)
pred = predict(mod, newdata=henv, block.size=block.size)
dims = list(lon=unique(henv$lon), lat=unique(henv$lat), time=unique(henv$time))
dim(pred) = sapply(dims, length)
outputFile = grep(x=outputFiles, pattern = iScen, value=TRUE)
write_ncdf(pred, filename=file.path(outputDir, outputFile), 
           varid="sst", dim=dims, longname="Sea Surface Temperature", units="ºC")    

for(iScen in scenarios) {
    
    kali::DateStamp("\tPredicting scenarios", iScen)
    penv[, modvar] = penv[[iScen]]
    penv[(penv$lon>-5.4 & penv$lat<43.2), modvar] = NA
    
    # future data
    pred = predict(mod, newdata=penv, block.size=block.size)
    dims = list(lon=unique(penv$lon), lat=unique(penv$lat), time=unique(penv$time))
    dim(pred) = sapply(dims, length)
    outputFile = grep(x=outputFiles, pattern = iScen, value=TRUE)
    write_ncdf(pred, filename=file.path(outputDir, outputFile), 
               varid="sst", dim=dims, longname="Sea Surface Temperature", units="ºC")    

}


rm(mod)

beepr::beep(4)
