library(nctools)
library(lubridate)

inputDir = file.path("ISIMIP2a/InputData/ocean/regional", domain)
outputDir = file.path("ISIMIP2a/InputData/ocean/downscaling", domain, "_input")

probs = unique(c(seq(0, 0.2, by=0.005), seq(0.2, 0.8, by=0.01), seq(0.8, 1, by=0.005)))

ncmask = nc_open(file.path(outputDir, "mask.nc4"))
mask = list(longitude = ncvar_get(ncmask, "lon"),
            latitude  = ncvar_get(ncmask, "lat"),
            mask      = ncvar_get(ncmask, "mask"))
nc_close(ncmask)

model = "gfdl-esm2m" # "*" for all models
variables = c("to_zs", "to_zb", "npp_zint", "so_zs") # "*" for all files, "to_zs" for SST
variables = "to_zb"
useLog = c(sst=FALSE, sss=FALSE, sbt=FALSE, mld=FALSE, npp=TRUE, do2=TRUE)

for(iVar in variables) {
  
  kali::DateStamp("Processing", sQuote(iVar), "...")
  
  files = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
  files = grep(x=files, patt=iVar, value = TRUE)
  files = grep(x=files, patt=model, value = TRUE)
  files = grep(x=files, patt="historical", value = TRUE)
  files = grep(x=files, patt="annual", value = TRUE, invert=TRUE)
  files = file.path(inputDir, files)
  
  origin = date("1950-01-01")
  starts = c(to=date("1981-09-01"), so=date("1980-01-01"), npp=date("1997-10-01"),
             do2=date("1958-01-01"), mld=date("1980-01-01"))
  
  timeStart = setNames(as.numeric(as.period(starts-origin)), names(starts))
  
  
  # Create quantile files for ESM data --------------------------------------
  
  
  for(file in files) {
    
    kali::DateStamp("\tProcessing file", file)
    
    varid  = strsplit(basename(file), split="_")[[1]][3]
    file = nc_regrid(filename=file, varid=varid, dim=c("longitude", "latitude"), 
                     log=useLog[iVar], new = mask, output=file.path(outputDir,basename(file)), 
                     extrap=TRUE)
    
    oFile = file.path(outputDir, gsub(x=basename(file), pattern=domain, 
                                      replacement=paste(domain, "quantile", sep="_")))
    oFile = gsub(x=oFile, pattern="195001|195901", replacement=format(starts[varid], "%Y%m"))
    nstart = interval(origin, starts[varid]) %/% months(1)
    
    kali::DateStamp("\t\tSubsetting time...")
    tmp = nc_subset(file, output = tempfile(), time=c(nstart, +Inf))
    
    kali::DateStamp("\t\tComputing quantiles...")
    nc_quantile(filename=tmp, varid = NA, MARGIN=c(1,2), 
                probs=probs, output=oFile, compression=9)
    
    file.remove(tmp)
    
  }
  
  
  # Future ------------------------------------------------------------------
  
  files = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
  files = grep(x=files, patt=iVar, value = TRUE)
  files = grep(x=files, patt=model, value = TRUE)
  files = grep(x=files, patt="historical", value = TRUE, invert=TRUE)
  files = grep(x=files, patt="annual", value = TRUE, invert=TRUE)
  files = file.path(inputDir, files)
  
  for(file in files) {
    
    kali::DateStamp("\tProcessing file", file)
    varid  = strsplit(basename(file), split="_")[[1]][3]
    file = nc_regrid(file=file, varid=varid, dim=c("longitude", "latitude"), log=useLog[iVar], 
                     new = mask, output=file.path(outputDir, basename(file)), 
                     extrap=TRUE)
    
  }
  
  
} # end of variables loop
