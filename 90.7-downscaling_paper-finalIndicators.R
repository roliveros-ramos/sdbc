library(readr)
library(scam)
library(foreach)
library(doSNOW)
library(ncdf4)
library(kali)

source("99.9-auxiliar_functions.R")

model = "ipsl" # "*" for all models
ivar = "to_zs"

modelDir = "output/final"
predDir = "output/final/predictions"
outputDir = "output/final/indicators"
figDir = "output/final/indicators/check"

dataSources = c(to_zs="avhrr", 
                npp="vgpm.s")

# Reading files -----------------------------------------------------------

inputDir = file.path(iDir, domain, "_input")

# quantile files
xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
xfiles = grep(x=xfiles, patt="quantile", value = TRUE)
file = grep(x=xfiles, patt=model, value = TRUE)
# future file
xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
xfiles = grep(x=xfiles, patt="2100", value = TRUE)
pfile = grep(x=xfiles, patt=model, value = TRUE)
# historical file
xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
xfiles = grep(x=xfiles, patt="1950", value = TRUE)
hfile = grep(x=xfiles, patt=model, value = TRUE)
# historical obs
xfiles = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=dataSources[ivar], value = TRUE)
ofile  = grep(x=xfiles, patt="quantile", value = TRUE, invert = TRUE)

# read files
varid = strsplit(basename(file), split="_")[[1]][3]
nenv = read_csv(file.path(inputDir, file), col_types = "dddddddd")
penv = read_csv(file.path(inputDir, pfile), col_types = "dddddddddd")
henv = read_csv(file.path(inputDir, hfile), col_types = "ddddddd")

nco = nc_open(file.path(inputDir, ofile))
oenv = ncvar_get(nco)
nc_close(nco)

# process files
slog = function(x) sign(x)*log(abs(x))
nenv$ldc = slog(nenv$dc)
nenv$lds = slog(nenv$ds)
penv$ldc = slog(penv$dc)
penv$lds = slog(penv$ds)
henv$ldc = slog(henv$dc)
henv$lds = slog(henv$ds)
# index files: historical from 198109
henv = henv[henv$time>12*(1981-1950)+(9-1), ]
names(henv)[6] = "sst_mod"
# index files: future from 208001
penv = penv[penv$time>12*(2080-2006)+(1-1), ]
names(penv)[9] = "sst_mod"
penv[, 6:8] = NULL

nenv = nenv[!(nenv$lon>-5.4 & nenv$lat<43.2), ] # remove mediterranean sea
penv$sst_mod[(penv$lon>-5.4 & penv$lat<43.2)] = NA
henv$sst_mod[(henv$lon>-5.4 & henv$lat<43.2)] = NA

if(mean(nenv$sst_mod, na.rm=TRUE)>200) nenv$sst_mod = nenv$sst_mod - 273.15

xpenv = penv$sst_mod
dims = list(lon=unique(penv$lon), lat=unique(penv$lat), time=unique(penv$time))
dim(xpenv) = sapply(dims, length)
xhenv = henv$sst_mod
dims = list(lon=unique(henv$lon), lat=unique(henv$lat), time=unique(henv$time))
dim(xhenv) = sapply(dims, length)

# the ESM without correction
esm = list(lon=dims$lon, lat=dims$lat, historical=xhenv, future=xpenv, quantile=nenv$sst_mod)
esm = addExtremes(esm)

# the observations
obs = list(lon=dims$lon, lat=dims$lat, historical=oenv, future=NULL, quantile=nenv$sst)
obs = addExtremes(obs)

yobs = density(nenv$sst)
ymod = density(nenv$sst_mod) 

# Model predictions -------------------------------------------------------

if(!dir.exists(outputDir)) dir.create(outputDir)
if(!dir.exists(figDir)) dir.create(figDir)

allModels = dir(path=modelDir, pattern=".*finalModel.*\\.rds")
allModels = file.path(modelDir, grep(x=allModels, pattern=domain, value=TRUE))

nullModel = grep(x=allModels, pattern="m0.0", value=TRUE)

ind00 = which(nenv$qua <= 0.00)
ind01 = which(nenv$qua <= 0.01)
ind05 = which(nenv$qua <= 0.05)
ind95 = which(nenv$qua >= 0.95)
ind99 = which(nenv$qua >= 0.99)
ind100 = which(nenv$qua >= 1)

indShelf = which(nenv$ds < 0)
indOpen  = which(nenv$ds > 0)

# null model (only one null model)
null = getPredData(nullModel, predDir, pattern="finalModel")
nullSE = (null$quantile - obs$quantile)^2
nullDev = sum(nullSE, na.rm=TRUE)
nullDevShelf = sum(nullSE[indShelf], na.rm=TRUE)
nullDevOpen  = sum(nullSE[indOpen], na.rm=TRUE)

nullDev00  = sum(nullSE[ind00], na.rm=TRUE)
nullDev01  = sum(nullSE[ind01], na.rm=TRUE)
nullDev05  = sum(nullSE[ind05], na.rm=TRUE)
nullDev95  = sum(nullSE[ind95], na.rm=TRUE)
nullDev99  = sum(nullSE[ind99], na.rm=TRUE)
nullDev100 = sum(nullSE[ind100], na.rm=TRUE)

esm_max = esm$max$historical - obs$max$historical
esm_min = esm$min$historical - obs$min$historical
esm_mean = esm$mean$historical - obs$mean$historical

presentBias = (esm$quantile - obs$quantile)
min_presentBias = presentBias[ind00]
max_presentBias = presentBias[ind100]

# zlim = range(c(esm_mean, esm_max, esm_min), na.rm=TRUE)
# biasPal = biasPalette(64, zlim=zlim, col=c("turquoise3", "red"))
# par(mfrow=c(1,3), mar=c(3,3,1,4), oma=c(1,1,1,1))
# image.map(esm$lon, esm$lat, esm_mean, col=biasPal, zlim=zlim)
# image.map(esm$lon, esm$lat, esm_min, col=biasPal, zlim=zlim)
# image.map(esm$lon, esm$lat, esm_max, col=biasPal, zlim=zlim)

# zlim = range(c(esm_max), na.rm=TRUE)
# biasPal1 = biasPalette(64, zlim=zlim)
# image.map(esm$lon, esm$lat, esm_max, col=biasPal1, zlim=zlim)

pb = txtProgressBar(style=3)
setTxtProgressBar(pb, value=0)

for(i in seq_along(allModels)) {

  indFile = gsub(x=basename(allModels[i]), pattern="finalModel", replacement="ind")
  indFile = gsub(x=indFile, pattern="rds$", replacement="txt")
  
  if(file.exists(file.path(outputDir, indFile))) next
  
  iMod = strsplit(basename(allModels[i]), "_")[[1]][3]
  thisModel = gsub(x=iMod, pattern="\\.rds", replacement = "")
  thisRep   = as.numeric(gsub(strsplit(iMod, "-")[[1]][2], pattern="rep", replacement=""))

  F1File = gsub(x=basename(allModels[i]), pattern="finalModel", replacement="check1")
  F1File = gsub(x=F1File, pattern="rds$", replacement="png")

  F2File = gsub(x=basename(allModels[i]), pattern="finalModel", replacement="check2")
  F2File = gsub(x=F2File, pattern="rds$", replacement="png")
  
  thisPred = getPredData(allModels[i], predDir, pattern="finalModel")

  # png(filename=file.path(figDir, F1File), width = 800, height = 640)
  # plot(thisPred)
  # dev.off()
  # png(filename=file.path(figDir, F2File), width = 600, height = 600)
  # density(thisPred, ymod=ymod, yobs=yobs, main=iMod)
  # dev.off()
  
  # indicators
  
  presentError = thisPred$quantile - obs$quantile # from quantiles
  futureBias  = esm$future - thisPred$future # ESM - DS
  
  bmax = thisPred$max$historical - obs$max$historical # pred error (max)
  bmin = thisPred$min$historical - obs$min$historical # pred error (min)
  
  cmax = esm$max$future - thisPred$max$future  # estimated future bias (max)
  cmin = esm$min$future - thisPred$min$future  # estimated future bias (min)
  
  ind = numeric()

# FUTURE INDICATORS -------------------------------------------------------
 
  # ind[01]: is the minimum estimated bias lower than minimum observed bias plus pred error?
  # lower than 1 mean true (for minimum tEMPERATURE) Are minimum TOO cold?
  ind[01] = abs(min(cmin, na.rm=TRUE))/(abs(min(esm_min, na.rm=TRUE)) + abs(min(bmin, na.rm=TRUE)))
  # ind[02]: is the maximum estimated bias lower than maximum observed bias plus pred error?
  # lower than 1 mean TRUE (for MAXIMUM TEMPERATURE). Are maximum TOO hot?
  ind[02] = abs(max(cmax, na.rm=TRUE))/(abs(max(esm_max, na.rm=TRUE)) + abs(max(bmax, na.rm=TRUE)))
  
  # total bias
  # relative change in bias future vs. present
  ind[03] = diff(range(futureBias, na.rm=TRUE))/diff(range(presentBias, na.rm=TRUE)) - 1
  # change in mean bias (future vs. present) 
  ind[04] = mean(futureBias, na.rm=TRUE) - mean(presentBias, na.rm=TRUE)
  # relative change in mean bias (future vs. present) 
  ind[05] = ind[04]/diff(range(presentBias, na.rm=TRUE))

  # bias for minima
  # relative change in bias future vs. present
  ind[06] = diff(range(cmin, na.rm=TRUE))/diff(range(min_presentBias, na.rm=TRUE)) - 1
  # change in mean bias (future vs. present) 
  ind[07] = mean(cmin, na.rm=TRUE) - mean(min_presentBias, na.rm=TRUE)
  # relative change in mean bias (future vs. present) 
  ind[08] = ind[07]/diff(range(min_presentBias, na.rm=TRUE))

  # bias for maxima
  # relative change in bias future vs. present
  ind[09] = diff(range(cmax, na.rm=TRUE))/diff(range(max_presentBias, na.rm=TRUE)) - 1
  # change in mean bias (future vs. present) 
  ind[10] = mean(cmax, na.rm=TRUE) - mean(max_presentBias, na.rm=TRUE)
  # relative change in mean bias (future vs. present) 
  ind[11] = ind[10]/diff(range(max_presentBias, na.rm=TRUE))
 
  # Bhattacharyya distance
  ind[12] = Bhattacharyya(futureBias, presentBias) # same bias distribution

# PRESENT INDICATORS ------------------------------------------------------

  # average magnitude of the absolute error 
  ind[13] = mean(abs(presentError), na.rm=TRUE) # average magnitude of error
  # range of the error 
  ind[14] = diff(range(presentError, na.rm=TRUE))
  # average bias in direction of the error
  ind[15] = mean(sign(presentError), na.rm=TRUE) # average direction of error
  # asymmetry indicator (0 for symmetric), detects extrapolation errors
  ind[16] = CTE(presentError[ind95],  presentError[ind05])
  ind[17] = CTE(presentError[ind99],  presentError[ind01])
  ind[18] = CTE(presentError[ind100], presentError[ind00])
  
  ind[19] = 1 - sum(presentError^2, na.rm=TRUE)/nullDev # deviance full dataset
  ind[20] = 1 - sum(presentError[indShelf]^2, na.rm=TRUE)/nullDevShelf # deviance inShelf
  ind[21] = 1 - sum(presentError[indOpen]^2, na.rm=TRUE)/nullDevOpen # deviance openwater
  
  ind[22] = 1 - sum(presentError[ind00]^2, na.rm=TRUE)/nullDev00 # deviance minima
  ind[23] = 1 - sum(presentError[ind01]^2, na.rm=TRUE)/nullDev01 # deviance 1%
  ind[24] = 1 - sum(presentError[ind05]^2, na.rm=TRUE)/nullDev05 # deviance 5%
  ind[25] = 1 - sum(presentError[ind95]^2, na.rm=TRUE)/nullDev95 # deviance 95%
  ind[26] = 1 - sum(presentError[ind99]^2, na.rm=TRUE)/nullDev99 # deviance 99%
  ind[27] = 1 - sum(presentError[ind100]^2, na.rm=TRUE)/nullDev100 # deviance maxima
  
  cat(c(thisModel, thisRep, ind), file=file.path(outputDir, indFile), sep=",")
  
  setTxtProgressBar(pb, value=i/length(allModels))
  
} # end for


beepr::beep(4)







