library(readr)
library(scam)
library(foreach)
library(doSNOW)

# source("99.9-auxiliar_functions.R")

RevoUtilsMath::setMKLthreads(n=1)

N = 24 # number of datasets
# n = 1e4 # sample size for each test dataset
model = "ipsl" # "*" for all models
ivar = "to_zs"

cores = 4
fitSCAM = TRUE

outputDir = "output"

nameCode = "model_%s-rep%02d_%s.rds"

# Reading files -----------------------------------------------------------

inputDir = file.path(iDir, domain, "_input")

if(!dir.exists(outputDir)) dir.create(outputDir)

xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
xfiles = grep(x=xfiles, patt="quantile", value = TRUE)

file = grep(x=xfiles, patt=model, value = TRUE)

xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
xfiles = grep(x=xfiles, patt="2100", value = TRUE)

pfile = grep(x=xfiles, patt=model, value = TRUE)

varid = strsplit(basename(file), split="_")[[1]][3]
nenv = read_csv(file.path(inputDir, file), col_types = "dddddddd")

# penv = read_csv(file.path(inputDir, pfile), col_types = "ddddddddd")

slog = function(x) sign(x)*log(abs(x))

nenv$ldc = slog(nenv$dc)
nenv$lds = slog(nenv$ds)

vars = names(nenv)
modvar = grep(x=vars, pattern="_mod", value=TRUE)
obsvar = gsub(x=modvar, pattern="_mod", replacement="")
aux = setdiff(vars, c(modvar, obsvar, "qua"))

varCheck = obsvar %in% vars

nenv = nenv[!(nenv$lon>-5.4 & nenv$lat<43.2), ] # remove mediterranean sea (43.2)

if(obsvar=="sst" & mean(nenv[[modvar]], na.rm=TRUE)>200) 
  nenv[, modvar] = nenv[, modvar] - 273.15


# Sampling data -----------------------------------------------------------

set.seed(123)

# sampling matrix (variance reduction)
n = nrow(nenv)
ind = matrix(sample(nrow(nenv)), nrow=ceiling(n/N), ncol=24)

# Model fitting -----------------------------------------------------------

# make cluster
cl = makeCluster(rep("localhost", cores))
registerDoSNOW(cl)

workDir = getwd()

# null model (only one null model)

modelName = "m0.0"
modFile = sprintf(nameCode, modelName, 0, domain)
if(!file.exists(file.path(outputDir, modFile))) {
  mod = gam(sst ~ offset(sst_mod), data=nenv)
  saveRDS(mod, file=file.path(outputDir, modFile))
}

# loop
foreach(i=seq_len(N), .inorder = FALSE, .packages = c("scam", "mgcv", "stats")) %dopar% {

  setwd(workDir)
  RevoUtilsMath::setMKLthreads(n=1)
  
  iDat = nenv[ind[, i], ]
  
  modelName = "m0.0"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam(sst ~ offset(sst_mod), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m0.1"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam(sst ~ offset(sst_mod) + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m0.2"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam(sst ~ sst_mod + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m1.0"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam( sst ~ s(sst_mod, bs="ad", k=40) + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  if(fitSCAM) {
    
    modelName = "m1.1"
    modFile = sprintf(nameCode, modelName, i, domain)
    if(!file.exists(file.path(outputDir, modFile))) {
      mod = scam(sst ~ s(sst_mod, bs="mpi", k=20)  + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat)
      saveRDS(mod, file=file.path(outputDir, modFile))
    }
    
    modelName = "m1.2"
    modFile = sprintf(nameCode, modelName, i, domain)
    if(!file.exists(file.path(outputDir, modFile))) {
      mod = scam(sst ~ s(sst_mod, bs="micx", k=20) + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat)
      saveRDS(mod, file=file.path(outputDir, modFile))
    }
    
    modelName = "m1.3"
    modFile = sprintf(nameCode, modelName, i, domain)
    if(!file.exists(file.path(outputDir, modFile))) {
      mod = scam(sst ~ s(sst_mod, bs="micv", k=20) + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat)
      saveRDS(mod, file=file.path(outputDir, modFile))
    }
    
  }  
  
  modelName = "m2.0"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam(sst ~ te(by=sst_mod, lon, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m2.1"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam(sst ~ te(by=sst_mod, lon, k=20) + te(by=sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m2.2"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam(sst ~ sst_mod + sst_mod:lon + sst_mod:lat + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m3.0"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam( sst ~ s(sst_mod, bs="ad", k=40) + ti(sst_mod, lon, k=20) + ti(sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  if(fitSCAM) {
    
    modelName = "m3.1"
    modFile = sprintf(nameCode, modelName, i, domain)
    if(!file.exists(file.path(outputDir, modFile))) {
      mod = scam(sst ~ s(sst_mod, bs="mpi", k=20) + ti(sst_mod, lon, k=20) + ti(sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat)
      saveRDS(mod, file=file.path(outputDir, modFile))
    }

    modelName = "m3.2"
    modFile = sprintf(nameCode, modelName, i, domain)
    if(!file.exists(file.path(outputDir, modFile))) {
      mod = scam(sst ~ s(sst_mod, bs="micx", k=20) + ti(sst_mod, lon, k=20) + ti(sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat)
      saveRDS(mod, file=file.path(outputDir, modFile))
    }
    
    modelName = "m3.3"
    modFile = sprintf(nameCode, modelName, i, domain)
    if(!file.exists(file.path(outputDir, modFile))) {
      mod = scam(sst ~ s(sst_mod, bs="micv", k=20) + ti(sst_mod, lon, k=20) + ti(sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat)
      saveRDS(mod, file=file.path(outputDir, modFile))
    }
    
  }
  
  modelName = "m3.4"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = try(gam(sst ~ te(sst_mod, lon, k=20) + te(sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat))
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m3.5"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = try(gam(sst ~ te(sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat))
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m3.6"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = try(gam(sst ~ te(sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat))
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m4.0"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam( sst ~ te(by=sst_mod, lon, lat, k=20) + s(sst_mod, bs="ad", k=40) + ti(sst_mod, lon, k=20) + ti(sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m4.1"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam( sst ~ te(by=sst_mod, lon, lat, k=20) + ti(sst_mod, lon, k=20) + ti(sst_mod, lat, k=20) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m4.2"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = gam( sst ~ te(by=sst_mod, lon, lat, k=20) + s(sst_mod, bs="ad", k=40) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  modelName = "m4.3"
  modFile = sprintf(nameCode, modelName, i, domain)
  if(!file.exists(file.path(outputDir, modFile))) {
    mod = try(gam(sst ~ te(lon, lat, sst_mod, k=10) + te(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML"))
    saveRDS(mod, file=file.path(outputDir, modFile))
  }
  
  NULL
  
}


stopCluster(cl)
beepr::beep(4)
