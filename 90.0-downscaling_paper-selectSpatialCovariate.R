library(readr)
library(scam)
library(RColorBrewer)

source("99.9-auxiliar_functions.R")


N = 12 # number of datasets
n = 1e5 # sample size for each test dataset
model = "ipsl" # "*" for all models
ivar = "to_zs"

cores = 4

RevoUtilsMath::setMKLthreads(n=cores)

outputDir = "output"

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

opath = file.path(outputDir, "selection")
if(!dir.exists(opath)) dir.create(opath, recursive = TRUE)


# Sampling data -----------------------------------------------------------

set.seed(123)

# sampling matrix (variance reduction)
ind = matrix(sample(nrow(nenv), 5*n*N, replace=TRUE), nrow=n, ncol=5*N)

# Model fitting -----------------------------------------------------------

i = 1
iDat = nenv[ind[,1], ]
# iDat = nenv

kali::DateStamp("Fitting models...")
null = gam(sst ~ offset(sst_mod), data=iDat)

kali::DateStamp("\tModel  1.0...")
mod1.0 = glm(sst ~ sst_mod + lon + lat + lon:lat + lds, data=iDat)
m1.0  = explainedDeviance(mod1.0, fixed=1, null=null)
rm()
kali::DateStamp("\tModel  2.0...")
mod2.0 = gam(sst ~ sst_mod + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
m2.0  = explainedDeviance(mod2.0, fixed=1, null=null)
kali::DateStamp("\tModel  3.0...")
mod3.0 = gam(sst ~ s(sst_mod, bs="cr", k=40) + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
m3.0  = explainedDeviance(mod3.0, fixed=1, null=null)
kali::DateStamp("\tModel  6.0a...")
mod6.0a = gam( sst ~ s(sst_mod, bs="cr", k=40) + ti(sst_mod, lat, k=20) + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
m6.0a = explainedDeviance(mod6.0a, fixed=1, null=null)
kali::DateStamp("\tModel  6.0b...")
mod6.0b = gam( sst ~ s(sst_mod, bs="cr", k=40) + ti(sst_mod, lds, k=20) + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
m6.0b = explainedDeviance(mod6.0b, fixed=1, null=null)
kali::DateStamp("\tModel  6.0c...")
mod6.0c = gam( sst ~ s(sst_mod, bs="cr", k=40) + ti(sst_mod, lon, k=20) + ti(lon, k=20) + ti(lat, k=30) + ti(lon, lat, k=15) + s(lds, bs="cr", k=15), data=iDat, method="REML")
m6.0c = explainedDeviance(mod6.0c, fixed=1, null=null)

beepr::beep(4)

kali::DateStamp("calculating deviance...")

ed = explDev(mod1.0, mod2.0, mod3.0, mod6.0a, mod6.0b, mod6.0c, null=null, sort=FALSE)


m1.0 = append(m1.0, rep(NA,3), 1)
m2.0 = append(m2.0, rep(NA,3), 1)
m3.0 = append(m3.0, rep(NA,3), 1)
m6.0a = append(append(m6.0a, NA, 1), NA, 3)
m6.0b = append(m6.0b, rep(NA,2), 1)
m6.0c = append(m6.0c, rep(NA,2), 2)

xnames = c(names(m6.0c)[2],names(m6.0a)[3], names(m6.0b)[4])
names(m6.0a)[2:4] = names(m6.0b)[2:4] = names(m6.0c)[2:4] = xnames
names(m1.0)[2:4] = paste("sst_mod", c("lon", "lat", "lds"), sep=":")

allM = cbind(m1.0, m2.0, m3.0, m6.0a, m6.0b, m6.0c)
allM[is.na(allM)] = 0
allD = t(t(allM)*ed)

write.csv(allM, file=file.path(opath, sprintf("%s_%s_explDevR.csv", domain, obsvar)))
write.csv(allD, file=file.path(opath, sprintf("%s_%s_explDev.csv", domain, obsvar)))

beepr::beep(4)
