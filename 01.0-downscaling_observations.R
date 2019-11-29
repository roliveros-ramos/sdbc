library(nctools)
library(mapdata)
library(lubridate)
library(readr)
library(kali)
library(raster)

inputDir = file.path("../envData/regional", domain)
outputDir = file.path("ISIMIP2a/InputData/ocean/downscaling", domain, "_input")

shelfFile  = "../spatial/shelfbreak/shelfbreak_200-0.0.1/shelfbreak_200.csv"
shelfShape = "../spatial/shelfbreak/shelfbreak_200-0.0.1/shelfbreak_200.shp"

coastFile = "world"

reDo = TRUE
reDoAux = FALSE
variables = "npp" # NULL for all dataSources

probs = unique(c(seq(0, 0.2, by=0.005), seq(0.2, 0.8, by=0.01), seq(0.8, 1, by=0.005)))

probs0 = seq(from=0, to=1, length=length(probs)) # test with no increase in tail sampling

if(!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

dataSources = c(sst="avhrr", 
                sss="soda3.3.1.*sss",
                sbt="soda3.3.1.*sbt",
                mld="soda3.3.1.*mlt",
                do2="roms-pisces.*do2",
                npp="vgpm.s")

# dataSources = c(sss="soda3.3.1.*sss") 

origin = date("1978-01-01 00:00:00") # origin for SODA and OI-SST
end = date("2005-12-31 23:59:59") # end of fishmip historical data

timeLimit = c(sst=as.numeric(end-origin), 
              sss=as.numeric(end-origin),
              sbt=as.numeric(end-origin),
              mld=as.numeric(end-origin),
              npp=2006,
              do2=(2006-1950)*12-1)

useLog = c(sst=FALSE, sss=FALSE, sbt=FALSE, mld=FALSE, npp=FALSE, do2=FALSE)
# use keepRange=TRUE if decreasing resolution when regriding
keepRange = c(sst=FALSE, sss=FALSE, sbt=FALSE, mld=FALSE, npp=TRUE, do2=FALSE)

# Grid calculation --------------------------------------------------------

if(is.null(gridFile)) {
  
  filename = grep(pattern=dataSources["sst"], x=dir(path=inputDir), value=TRUE)
  DateStamp("Grid file is NULL. Using", filename)
  x = nc_mask(filename=file.path(inputDir, filename),
              output=file.path(outputDir, "mask.nc4"))
  
} else {
  
  DateStamp("Using grid file ", gridFile)
  x = nc_mask(filename=gridFile,
              output=file.path(outputDir, "mask.nc4"))
  
}

res = 72
r = diff(range(x$lon))/diff(range(x$lat))
tiff(filename = file.path(outputDir, "mask.tif"), height=10*res, width=10*r*res)
par(mar=c(0,0,0,0), oma=c(3,3,1,1))
plot.new()
plot.window(xlim=range(x$lon), ylim=range(x$lat), asp=1, xaxs="i", yaxs="i")
image(x$lon, x$lat, x$mask, xlab="LONGITUDE", ylab="LATITUDE", 
      col="lightblue", axes=FALSE, add=TRUE)
map("worldHires", add=TRUE, interior=FALSE)
kali::map.axes2()
box()
dev.off()

ncmask = nc_open(file.path(outputDir, "mask.nc4"))
mask = list(lon  = ncvar_get(ncmask, "lon"),
            lat  = ncvar_get(ncmask, "lat"),
            mask = ncvar_get(ncmask, "mask"))
nc_close(ncmask)

# Quantile calculation ----------------------------------------------------

DateStamp("Processing observations for domain", sQuote(domain))
variables = if(is.null(variables)) names(dataSources) else variables
for(iVar in variables) {
  
  kali::DateStamp("Processing", sQuote(iVar), "...")
  
  files = unlist(sapply(dataSources[iVar], x=dir(path=inputDir), value=TRUE, FUN=grep))
  
  if(length(files)==0) next 
  
  for(file in files) {
    
    oFile = file.path(outputDir, gsub(x=file, pattern=domain, 
                                      replacement=paste(domain, "quantile", sep="_")))
    oFile = gsub(x=oFile, pattern="([-|_])[0-9].*\\.", 
                 replacement=sprintf("\\1%s.", format(end, "%Y%m")))
    
    # test with simple quantile sampling
    oFile2 = file.path(outputDir, gsub(x=file, pattern=domain, 
                                      replacement=paste(domain, "quantile2", sep="_")))
    oFile2 = gsub(x=oFile2, pattern="([-|_])[0-9].*\\.", 
                 replacement=sprintf("\\1%s.", format(end, "%Y%m")))
    
    iFile = file.path(outputDir, gsub(x=file, pattern="([-|_])[0-9].*\\.", 
                                      replacement=sprintf("\\1%s.", format(end, "%Y%m"))))
    
    tmp = nc_subset(file.path(inputDir, file), output = tempfile(), time=c(-Inf, timeLimit[iVar]))
   
    if(!file.exists(iFile) | isTRUE(reDo)) {
      kali::DateStamp("\tInterpolating...")
      iFile = nc_regrid(file=tmp, varid=NA, dim=1:2, log=useLog[iVar], new = mask, 
                        output=iFile, extrap=TRUE, keepRange = keepRange[iVar])
      kali::DateStamp("\tFinished.")
    }
    
    if(!file.exists(oFile) | isTRUE(reDo)) {
      kali::DateStamp("\tCalculating quantiles...")
      nc_quantile(filename=iFile, varid = NA, MARGIN=c(1,2), probs=probs, output=oFile, compression=9)
      
    } 
    
    # if(!file.exist(oFile2)) 
    #   nc_quantile(filename=iFile, varid = NA, MARGIN=c(1,2), probs=probs0, output=oFile2, compression=9)
    
    
    file.remove(tmp)
    
  }
  
}


# Auxiliar data calculation -----------------------------------------------

if(!file.exists(file.path(outputDir, "auxiliar_info.csv")) | reDoAux) {
  
  kali::DateStamp("Processing auxiliar information...")
  
  aux = expand.grid(lon = x$lon,
                    lat = x$lat)
  
  shelf = read.csv(shelfFile)[,-1]
  coast = as.data.frame(setNames(map(coastFile, interior=FALSE, plot=FALSE)[c("x", "y")], 
                                 c("lon", "lat")))
  coast = coast[coast$lon<180.1, ]
  coast = coast[complete.cases(coast), ]
  
  sh = shapefile(shelfShape)
  
  aux$ds = getDistance(data=aux, ref=shelf)
  aux$dc = getDistance(data=aux, ref=coast)
  
  pos = SpatialPoints(cbind(lon=aux[, "lon"], lat=aux[, "lat"]), 
                       proj4string=CRS(proj4string(sh)))
  
  inShelf = ifelse(is.na(over(pos, sh)), 1, -1)
  aux$ds  = aux$ds*inShelf
  aux$mask = as.numeric(x$mask)
  
  write_csv(aux, path=file.path(outputDir, "auxiliar_info.csv"))
  
}


beepr::beep(4)
