library(readr)
library(scam)
library(foreach)
library(doSNOW)

source("99.9-auxiliar_functions.R")

N = 12 # number of datasets
n = 1e4 # sample size for each test dataset
model = "ipsl" # "*" for all models
ivar = "to_zs"

cores = 4

modelDir = "output/final"
outputDir = "output/final/predictions"

# Reading files -----------------------------------------------------------

inputDir = file.path(iDir, domain, "_input")

if(!dir.exists(outputDir)) dir.create(outputDir)
if(!dir.exists(file.path(outputDir, "check"))) dir.create(file.path(outputDir, "check"))

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
# index files: historical from 198109
henv = henv[henv$time>12*(1981-1950)+(9-1), ]
names(henv)[7] = "sst_mod"
# index files: future from 208001
fenv = penv
penv = penv[penv$time>12*(2080-2006)+(1-1), ]
names(penv)[10] = "sst_mod"
penv[, 7:9] = NULL

nenv = nenv[!(nenv$lon>-5.4 & nenv$lat<43.2), ] # remove mediterranean sea
penv$sst_mod[(penv$lon>-5.4 & penv$lat<43.2)] = NA
henv$sst_mod[(henv$lon>-5.4 & henv$lat<43.2)] = NA

if(mean(nenv$sst_mod, na.rm=TRUE)>200) nenv$sst_mod = nenv$sst_mod - 273.15


# Model predictions -------------------------------------------------------

if(!dir.exists(outputDir)) dir.create(outputDir)

allModels = dir(path=modelDir, pattern=".*finalModel.*\\.rds")
allModels = file.path(modelDir, grep(x=allModels, pattern=domain, value=TRUE))

# make cluster
cl = makeCluster(rep("localhost", cores))
registerDoSNOW(cl)

workDir = getwd()

# null model (only one null model)

# loop
foreach(i=seq_along(allModels), .inorder = FALSE, .packages = c("scam", "mgcv", "stats", "nctools")) %dopar% {

  setwd(workDir)

  RevoUtilsMath::setMKLthreads(n=1)
  
  mod = readRDS(allModels[i])
  
  summFile = gsub(x=basename(allModels[i]), pattern="rds$", replacement="txt")
  
  if(!file.exists(file.path(outputDir, "check", summFile))) {
    sink(file = file.path(outputDir, "check", summFile))
    cat("\n", basename(allModels[i]), "\n")
    print(summary(mod))
    sink()
  }
  
  # historical data
  outputFile = gsub(x=basename(allModels[i]), pattern="finalModel", replacement="predH")
  outputFile = gsub(x=outputFile, pattern="rds$", replacement="nc4")
  if(!file.exists(file.path(outputDir, outputFile))) {
    pred = predict(mod, newdata=henv, block.size=120000)
    dims = list(lon=unique(henv$lon), lat=unique(henv$lat), time=unique(henv$time))
    dim(pred) = sapply(dims, length)
    write.ncdf(pred, filename=file.path(outputDir, outputFile), 
               varid="sst", dim=dims, longname="Sea Surface Temperature", units="ºC")    
  }

  # future data
  outputFile = gsub(x=basename(allModels[i]), pattern="finalModel", replacement="predF")
  outputFile = gsub(x=outputFile, pattern="rds$", replacement="nc4")
  if(!file.exists(file.path(outputDir, outputFile))) {
    pred = predict(mod, newdata=penv, block.size=120000)
    dims = list(lon=unique(penv$lon), lat=unique(penv$lat), time=unique(penv$time))
    dim(pred) = sapply(dims, length)
    write.ncdf(pred, filename=file.path(outputDir, outputFile), 
               varid="sst", dim=dims, longname="Sea Surface Temperature", units="ºC")    
  }
  
  # quantile data
  outputFile = gsub(x=basename(allModels[i]), pattern="finalModel", replacement="predQ")
  if(!file.exists(file.path(outputDir, outputFile))) {
    pred = predict(mod, newdata=nenv, block.size=120000)
    dimnames(pred) = NULL
    saveRDS(pred, file=file.path(outputDir, outputFile))    
  }
  
  rm(mod)
  
} # end foreach


stopCluster(cl)
beepr::beep(4)
