library(readr)
library(scam)

source("99.9-auxiliar_functions.R")

modelDir = "output"

model = "ipsl" # "*" for all models
ivar = "to_zs"

models = c("m0.0", "m0.2", "m1.0", "m2.0", "m3.0", "m4.0", "m4.2", "m1.3")

modelsHum = modelsBoB = models

fitSCAM = FALSE

# Reading files -----------------------------------------------------------

inputDir = file.path(iDir, domain, "_input")

xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
xfiles = grep(x=xfiles, patt="quantile", value = TRUE)

file = grep(x=xfiles, patt=model, value = TRUE)

xfiles = dir(path=inputDir, patt="\\.csv$", recursive=TRUE)
xfiles = grep(x=xfiles, patt=ivar, value = TRUE)
xfiles = grep(x=xfiles, patt="2100", value = TRUE)

varid = strsplit(basename(file), split="_")[[1]][3]
nenv = read_csv(file.path(inputDir, file), col_types = "dddddddd")

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

# Reading model files -----------------------------------------------------

files = dir(path=modelDir, pattern="rds$")
files = grep(x=files, pattern=domain, value=TRUE)

model = sapply(strsplit(x=files, split="-"), FUN="[", i=1)

mod1 = tapply(X = files, INDEX = model, FUN = "[", i=1)
mods = gsub(x=tapply(X = model, INDEX = model, FUN = unique),
            pattern="model_", replacement = "")

selectedModels = if(domain=="biscay-celtic") modelsBoB else modelsHum

nmod = as.numeric(gsub(x=selectedModels, pattern="m([0-9])\\..*", replacement="\\1"))

bs = as.numeric(sapply(strsplit(x=selectedModels, split="\\."), FUN="[", 2))


outputDir = file.path(modelDir, "final")
if(!dir.exists(outputDir)) dir.create(outputDir)

n = nrow(nenv)
xind0 = sample(nrow(nenv), size=0.90*n)
xind1 = sample(nrow(nenv), size=0.75*n)

sps = tapply(X = files, INDEX = model, FUN = readSps, path=modelDir)
names(sps) = names(mod1) = mods

# Fitting models ----------------------------------------------------------

for(i in seq_along(selectedModels)) {
  
  kali::DateStamp("Fitting", selectedModels[i])
  
  thisMod = selectedModels[i]
  
  if(thisMod!="m0.0") thisSp = apply(sps[[thisMod]], 2, median)
  
  refModel = readRDS(file.path(modelDir, mod1[thisMod]))
  modelType = class(refModel)[1]
  thisFormula = update(formula(refModel), ~ . - s(lds, bs="cr", k=15) + s(lds, bs="ad", k=40))
  
  if(modelType=="gam") {
    
    modFile = sprintf("%s_finalModel_%s.rds", domain, thisMod)
    if(!file.exists(file.path(outputDir, modFile))) {
      if(thisMod!="m0.0") {
        mod = bam(formula = thisFormula, data=nenv, method="fREML")
      } else {
        mod = gam(formula = thisFormula, data=nenv, method="REML")
      }
      saveRDS(mod, file=file.path(outputDir, modFile))
      rm(mod)
    }
    
  }
  
  if(modelType=="scam") {
    
    if(!fitSCAM) next
    
    gc()
    
    modFile = sprintf("%s_finalModel_%s.rds", domain, thisMod)
    

    # first pass with all data
    if(!file.exists(file.path(outputDir, modFile))) {
      # uses sp to speed up parameter estimation
      try({
        mod = scam(formula = thisFormula, data=nenv, sp=thisSp)
        saveRDS(mod, file=file.path(outputDir, modFile))
        rm(mod)
      })
    }
   
    gc()
    
    # second pass with 90% of data
    if(!file.exists(file.path(outputDir, modFile))) {
      # uses sp to speed up parameter estimation
      try({
        xenv = nenv[xind0, ]
        mod = scam(formula = thisFormula, data=xenv, sp=thisSp)
        save(mod, file=file.path(outputDir, modFile))
        rm(mod)
        rm(xenv)
        gc()
      })
    }
    
    # third pass with 75% of data
    if(!file.exists(file.path(outputDir, modFile))) {
      # uses sp to speed up parameter estimation
      try({
        xenv = nenv[xind1, ]
        mod = scam(formula = thisFormula, data=xenv, sp=thisSp)
        save(mod, file=file.path(outputDir, modFile))
        rm(mod)
        rm(xenv)
        gc()
      })
    }
    
  }
 
    gc()
   
}



