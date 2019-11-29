library(nctools)

inputDir = "ISIMIP2a/InputData/ocean"
patt = "*" # "*" for all files
model = "*" # "*" for all models
dx = 1
dy = 1

reDo = FALSE # re-extract all files

oldDim = c("LONGITUDE", "LATITUDE", "TIME", "DEPTH1_1", "zt_ocean", "DEPTHU1_1", "DEPTHV1_1", "DEPTHW1_1")
newDim = c("longitude", "latitude", "time", "depth", "depth", "depth", "depth", "depth")


# Extraction! -------------------------------------------------------------

xlat = lat + c(-1, +1)*dy
xlon = lon + c(-1, +1)*dx

pm = findPrimeMeridian(xlon)
harmo = if(pm=="center") "harmonized2" else "harmonized"
inputDir = file.path(inputDir, harmo)

files = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
files = grep(x=files, patt=patt, value = TRUE)
files = grep(x=files, patt=model, value = TRUE)
files = grep(x=files, patt="annual", value = TRUE, invert=TRUE)
files = file.path(inputDir, files)

rfiles = grep(x=files, patt="gfdl_reanalysis", value = TRUE)
files = grep(x=files, patt="gfdl_reanalysis", value = TRUE, invert=TRUE)

kali::DateStamp("Extracting data for", domain)

outputDir = file.path("regional", domain)

for(file in files) {

  kali::DateStamp("\tProcessing file", basename(file))
  iDir = dirname(file)
  oDir = gsub(x=iDir, patt=harmo, rep=outputDir)
  if(!dir.exists(oDir)) dir.create(oDir, recursive = TRUE)
  
  oFile = unlist(strsplit(x=basename(file), split="_"))
  oFile = paste(append(oFile, domain, after=4), collapse="_")
  ovar  = strsplit(oFile, split="_")[[1]][3]
  
  if(!file.exists(file.path(oDir, oFile)) | reDo) {
    
    tmp = nc_open(file)
    varid = names(tmp$var)[grep(x=tolower(names(tmp$var)), pattern=ovar)]
    nc_close(tmp)
    
    nc_subset(filename=file, varid = varid, newvarid = ovar, compression=9,
              output=file.path(oDir, oFile), latitude=xlat, longitude=xlon, 
              ignore.case = TRUE)
    
    nc_rename(filename=file.path(oDir, oFile), 
              oldnames = oldDim, newnames = newDim, 
              overwrite = TRUE, verbose = TRUE)
    
  }
  
}


for(file in rfiles) {
  
  kali::DateStamp("\tProcessing file", basename(file))
  iDir = dirname(file)
  oDir = gsub(x=iDir, patt=harmo, rep=outputDir)
  if(!dir.exists(oDir)) dir.create(oDir, recursive = TRUE)
  
  oFile = unlist(strsplit(x=basename(file), split="_"))
  oFile = paste(append(oFile, domain, after=4), collapse="_")
  oFile = gsub(x=oFile, pattern="1959", replacement="195901")
  oFile = gsub(x=oFile, pattern="2004", replacement="200412")
  ovar  = strsplit(oFile, split="_")[[1]][3]
  
  if(!file.exists(file.path(oDir, oFile)) | reDo) {
  
    
    tmp = nc_open(file)
    varid = names(tmp$var)[grep(x=tolower(names(tmp$var)), pattern=ovar)]
    nc_close(tmp)
    
    nc_subset(filename=file, varid = varid, newvarid = ovar, compression=9,
              output=file.path(oDir, oFile), latitude=xlat, longitude=xlon, 
              ignore.case = TRUE)
    
    nc_rename(filename=file.path(oDir, oFile), 
              oldnames = oldDim, newnames = newDim, 
              overwrite = TRUE, verbose = TRUE)
    
    # change time origin from 1959-01-01 to 1950-01-01 to match GFDL and IPSL
    
    nc = nc_open(file.path(oDir, oFile), write=TRUE)
    ncvar_put(nc, "time", ncvar_get(nc, "time") + 9*12)
    ncatt_put(nc, "time", "units", "months since 1950-01-01 00:00:00")
    nc_close(nc)
    
  }
    
}

kali::DateStamp("Done.")

