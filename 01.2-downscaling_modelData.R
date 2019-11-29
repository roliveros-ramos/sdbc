library(nctools)
library(mapdata)
library(readr)

# Input configuration -----------------------------------------------------

inputDir  = file.path("ISIMIP2a/InputData/ocean/downscaling", domain, "_input")
variables = NULL # NULL for all variables
models    = c("ipsl", "gfdl")

obsData = c(to_zs="avhrr", 
            so_zs="soda3.3.1.*sss",
            to_zb="soda3.3.1.*sbt",
            # mld="soda3.3.1.*mlt",
            npp_zint="vgpm.s")

modVars = c(to_zs="sst", 
            so_zs="sss",
            to_zb="sbt",
            # mld="mld",
            npp_zint="npp")

ncmask = nc_open(file.path(inputDir, "mask.nc4"))
aux = read_csv(file.path(inputDir, "auxiliar_info.csv"))
variables = if(is.null(variables)) names(obsData) else variables

# Model training data -----------------------------------------------------

kali::DateStamp("Creating model training databases...")
# TRAINING DATA
# create csv with data to fit the models
# one file per var and ESM, observations and auxiliar info added.

xfiles = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
xfiles = grep(x=xfiles, patt="quantile", value = TRUE)

for(ivar in variables) {
  
  kali::DateStamp("\tProcessing variable", ivar)
  
  modvar = modVars[ivar]
  obsFile = file.path(inputDir, grep(x=xfiles, pattern=obsData[ivar], value=TRUE))

  for(iModel in models) {

    files = grep(x=xfiles, patt=paste0("_", ivar), value = TRUE)
    files = grep(x=files, patt=iModel, value = TRUE)
    files = file.path(inputDir, files)
    
    for(file in files) {
      
      kali::DateStamp("\t\tProcessing", file, "...")
      oFile1 = gsub(x=basename(file), patt="\\.(nc|nc4)$", replacement=".RData")
      oFile2 = gsub(x=basename(file), patt="\\.(nc|nc4)$", replacement=".csv")
      
      nc = nc_open(file)
      msize = nc$var[[1]]$size
      xdims = ncvar_dim(nc, 1)
      dat = expand.grid(lon = ncvar_get(nc, xdims[1]),
                        lat = ncvar_get(nc, xdims[2]),
                        qua = ncvar_get(nc, xdims[3]))
      dat[[paste(modvar, "_mod", sep="")]] = as.numeric(ncvar_get(nc))
      nc_close(nc)
      
      nc = nc_open(obsFile)
      osize = nc$var[[1]]$size
      
      if(identical(msize, osize)) dat[[modvar]] = as.numeric(ncvar_get(nc))
      
      # merge with auxiliar data
      checkLon = identical(ncvar_get(ncmask, "lon"), ncvar_get(nc, "lon"))
      checkLat = identical(ncvar_get(ncmask, "lat"), ncvar_get(nc, "lat"))
      
      dat$index = seq_len(nrow(dat))
      
      if(checkLon & checkLat) dat = merge(dat, aux)
      
      dat = dat[order(dat$index), ]
      dat$index = NULL
      dat = dat[complete.cases(dat), ]
      
      save(dat, file=file.path(inputDir, oFile1))
      write_csv(dat, path=file.path(inputDir, oFile2))
      
      rm(dat)
      invisible(gc())
      
    } # end of files loop
    
  } # end of model loop
    
} # end of variable loop


# Prediction data - Future ------------------------------------------------

# PREDICTION DATA: FUTURE.
# Creates csv with ESM data and auxiliar information

kali::DateStamp("Creating future prediction database...")

xfiles = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
xfiles = grep(x=xfiles, patt="mask.nc4$", value = TRUE, invert=TRUE)
xfiles = grep(x=xfiles, patt="quantile", value = TRUE, invert=TRUE)
xfiles = grep(x=xfiles, patt="rcp", value = TRUE)


nc = nc_open(file.path(inputDir, xfiles[1]))
msize = nc$var[[1]]$size
dat0 = expand.grid(lon = ncvar_get(nc, "longitude"),
                   lat = ncvar_get(nc, "latitude"),
                   time = ncvar_get(nc, "time"))

checkLon = identical(ncvar_get(ncmask, "lon"), ncvar_get(nc, "longitude"))
checkLat = identical(ncvar_get(ncmask, "lat"), ncvar_get(nc, "latitude"))
nc_close(nc)

dat0$index = seq_len(nrow(dat0))

if(checkLon & checkLat) dat0 = merge(dat0, aux)

dat0 = dat0[order(dat0$index), ]
dat0$index = NULL

for(ivar in variables) {

  kali::DateStamp("\tProcessing variable", ivar)
  
  modvar = modVars[ivar]
  
  for(iModel in models) {
    
    dat = dat0
    files = grep(x=xfiles, patt=ivar, value = TRUE)
    files = grep(x=files, patt=iModel, value = TRUE)
    files = file.path(inputDir, files)
    
    mfile = sapply(lapply(strsplit(basename(files), split="_"), "[", i=-2), paste, collapse="_")[1]
    oFile1 = gsub(x=mfile, patt="\\.(nc|nc4)$", replacement=".RData")
    oFile2 = gsub(x=mfile, patt="\\.(nc|nc4)$", replacement=".csv")
    
    for(file in files) {
      
      kali::DateStamp("\t\tProcessing", file, "...")
      
      sce = strsplit(basename(file), split="_")[[1]][2]
      
      nc = nc_open(file)
      tmp = as.numeric(ncvar_get(nc))
      shift = if(modvar=="sst" & mean(tmp, na.rm=TRUE)>200) 273.15 else 0
      dat[[sce]] = tmp - shift
      nc_close(nc)
      
    } # end of files loop
    
    kali::DateStamp("\t\t\tSaving", oFile1, "...")
    save(dat, file=file.path(inputDir, oFile1))
    kali::DateStamp("\t\t\tSaving", oFile2, "...")
    write_csv(dat, path=file.path(inputDir, oFile2))
    
    rm(dat)
    invisible(gc())
    
  } # end of model loop
  
} # end of variable loop


# Prediction data - Historical --------------------------------------------

# PREDICTION DATA: HISTORICAL.
# Creates csv with ESM data and auxiliar information

xfiles = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
xfiles = grep(x=xfiles, patt="mask.nc4$", value = TRUE, invert=TRUE)
xfiles = grep(x=xfiles, patt="quantile", value = TRUE, invert=TRUE)
xfiles = grep(x=xfiles, patt="historical", value = TRUE)

kali::DateStamp("Creating historical prediction database...")

nc = nc_open(file.path(inputDir, xfiles[1]))
msize = nc$var[[1]]$size
dat0 = expand.grid(lon = ncvar_get(nc, "longitude"),
                   lat = ncvar_get(nc, "latitude"),
                   time = ncvar_get(nc, "time"))

checkLon = identical(ncvar_get(ncmask, "lon"), ncvar_get(nc, "longitude"))
checkLat = identical(ncvar_get(ncmask, "lat"), ncvar_get(nc, "latitude"))

nc_close(nc)

dat0$index = seq_len(nrow(dat0))

if(checkLon & checkLat) dat0 = merge(dat0, aux)

dat0 = dat0[order(dat0$index), ]
dat0$index = NULL

for(ivar in variables) {
  
  kali::DateStamp("\tProcessing variable", ivar)
  
  modvar = modVars[ivar]
  
  for(iModel in models) {
    
    dat = dat0
    files = grep(x=xfiles, patt=ivar, value = TRUE)
    files = grep(x=files, patt=iModel, value = TRUE)
    files = file.path(inputDir, files)
    
    mfile = sapply(lapply(strsplit(basename(files), split="_"), "[", i=-2), paste, collapse="_")[1]
    oFile1 = gsub(x=mfile, patt="\\.(nc|nc4)$", replacement=".RData")
    oFile2 = gsub(x=mfile, patt="\\.(nc|nc4)$", replacement=".csv")
    
    for(file in files) {
      
      kali::DateStamp("\t\tProcessing", file, "...")
      
      sce = strsplit(basename(file), split="_")[[1]][2]
      
      nc = nc_open(file)
      tmp = as.numeric(ncvar_get(nc))
      shift = if(modvar=="sst" & mean(tmp, na.rm=TRUE)>200) 273.15 else 0
      dat[[sce]] = tmp - shift
      nc_close(nc)
      
    }
    
    kali::DateStamp("\t\t\tSaving", oFile1, "...")
    save(dat, file=file.path(inputDir, oFile1))
    kali::DateStamp("\t\t\tSaving", oFile2, "...")
    write_csv(dat, path=file.path(inputDir, oFile2))
    
    rm(dat)
    invisible(gc())
    
  }
  
  
}

nc_close(ncmask)