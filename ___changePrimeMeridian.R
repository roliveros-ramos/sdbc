library(nctools)

inputDir = "ISIMIP2a/InputData/ocean/harmonized"
patt = "*" # "*" for all files
model = "*" # "*" for all models
outputDir = "harmonized2"

# Extraction! -------------------------------------------------------------

files = dir(path=inputDir, patt="\\.nc4$", recursive=TRUE)
files = grep(x=files, patt=patt, value = TRUE)
files = grep(x=files, patt=model, value = TRUE)
files = grep(x=files, patt="annual", value = TRUE, invert=TRUE)
files = c(grep(x=files, patt="_zall_", value = TRUE, invert=TRUE),
          grep(x=files, patt="_zall_", value = TRUE))
files = c(grep(x=files, patt="_z_", value = TRUE, invert=TRUE),
          grep(x=files, patt="_z_", value = TRUE))
files = file.path(inputDir, files)

kali::DateStamp("Starting...")

for(file in files) {
  
  kali::DateStamp("\tProcessing file", basename(file))
  iDir = dirname(file)
  oDir = gsub(x=iDir, patt="harmonized", rep=outputDir)
  if(!dir.exists(oDir)) dir.create(oDir, recursive = TRUE)
  
  oFile = basename(file)
  ovar  = strsplit(oFile, split="_")[[1]][3]
  
  tmp = nc_open(file)
  varid = names(tmp$var)[grep(x=tolower(names(tmp$var)), pattern=ovar)]
  nc_close(tmp)
  
  if(!file.exists(file.path(oDir, oFile))) {

    try(nc_changePrimeMeridian(filename=file, output=file.path(oDir, oFile), 
                           varid=varid, primeMeridian = "center",
                           overwrite=TRUE, compression=9))
    
    invisible(gc())
    
  }
  
}

kali::DateStamp("Done.")

