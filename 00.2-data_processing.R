library(nctools)
library(lubridate)
library(fields)

inputDir = file.path("ISIMIP2a/InputData/ocean/regional", domain)

# NPP ---------------------------------------------------------------------

files = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
files = grep(x=files, patt="spp_zint", value = TRUE)
files = grep(x=files, patt="annual", value = TRUE, invert=TRUE)
files = file.path(inputDir, files)

for(iFile in files) {
  
  kali::DateStamp("Processing file", iFile)
  oFile = gsub(x=iFile, pattern="spp_zint", replacement="npp_zint")
  lFile = gsub(x=iFile, pattern="spp_zint", replacement="lpp_zint")
  oTemp = paste0(oFile, ".temp")
  
  file.copy(from=iFile, to=oTemp, overwrite = TRUE)
  nc  = nc_open(oTemp, write = TRUE)
  ncl = nc_open(lFile)
  spp = ncvar_get(nc)
  lpp = ncvar_get(ncl)
  # npp to mgC/m2/day (molC/m2/s)
  npp = 12000*86400*(spp+lpp) # mol->mg (12000), s->day (86400)
  ncvar_put(nc, vals=npp)
  ncatt_put(nc, varid=names(nc$var), attname = "units", attval = "mgC/m2/day")
  ncatt_put(nc, varid=names(nc$var), attname = "long_name", attval = "PP Phytoplankton")
  nc_close(nc)
  nc_close(ncl)
  nc_rename(oTemp, oldnames="spp", newnames="npp", output=oFile, overwrite=TRUE)
  file.remove(oTemp)
}


# DO2 ---------------------------------------------------------------------
# 22 umol.kg-1 = 0.5 ml.L-1 : oxy=o2[d=1,z=@loc:22]
# 22 umol/L = 22 mmol/m3 = 0.022 mol/m3, 

files = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
files = grep(x=files, patt="o2_z(all)?_", value = TRUE)
files = grep(x=files, patt="annual", value = TRUE, invert=TRUE)
files = file.path(inputDir, files)

loc = 0.022
MARGIN = "depth"

for(iFile in files) {

  kali::DateStamp("Processing file", iFile)
  oFile = gsub(x=iFile, pattern="o2_z(all)?_", replacement="do2_")
  nc_loc(iFile, varid=NA, MARGIN=MARGIN, loc=loc, output=oFile,
       name="do2")
  
}


# T15 ---------------------------------------------------------------------

files = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
files = grep(x=files, patt="to_z(all)?_", value = TRUE)
files = grep(x=files, patt="annual", value = TRUE, invert=TRUE)
files = file.path(inputDir, files)

loc = 15
MARGIN = "depth"

for(iFile in files) {
  
  kali::DateStamp("Processing file", iFile)
  
  ncx = nc_open(iFile)
  xtemp = mean(ncvar_get(ncx), na.rm=TRUE)
  shift = if(xtemp>200) 273.15 else 0
  nc_close(ncx)
  
  oFile = gsub(x=iFile, pattern="to_z(all)?_", replacement="dt15_")
  nc_loc(iFile, varid=NA, MARGIN=MARGIN, loc=loc+shift, output=oFile,
         name="dt15")
  
}


# MLD ---------------------------------------------------------------------


# SBT ---------------------------------------------------------------------

# Missing data for GFDL-ESM2M

files = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
files = grep(x=files, patt="to_z(all)?_", value = TRUE)
files = grep(x=files, patt="historical", value = TRUE)
files = grep(x=files, patt="gfdl-esm2m", value = TRUE)
files = grep(x=files, patt="annual", value = TRUE, invert=TRUE)
files = file.path(inputDir, files)

iFile = files[1]

kali::DateStamp("Processing file", iFile)


oFile = gsub(x=iFile, pattern="to_z(all)?_", replacement="to_zb_")

.bottom = function(x) {
  if(all(is.na(x))) return(NA)
  out = tail(x[!is.na(x)], 1)
  return(out)
}

nc_apply(iFile, varid = NA, MARGIN=c(1,2,4), FUN = .bottom, output = oFile)


