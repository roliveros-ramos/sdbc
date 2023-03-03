
# explDev = function(...) {
#   
#   names = as.character(match.call()[-1])
#   
#   .explDev = function(object) 
#     (object$null.deviance - object$deviance)/object$null.deviance
#   
#   out = sapply(list(...), FUN=.explDev)
#   names(out) = names
#   out = sort(out, decreasing = TRUE)
#   return(out)
# }



explDev3 = function(..., newdata, y, null) {
  
  names = grep(as.character(match.call(expand.dots = FALSE)[-1]), 
               pattern="pairlist", value=TRUE)
  names = strsplit(gsub(x=names, pattern=".*\\((.*)\\)", rep="\\1"), split=", ")[[1]]
  
  pnull = predict(null, type="response", newdata=newdata)
  nulldev = sum(null$family$dev.resids(y=y, mu=pnull, wt = 1))
  kali::DateStamp("Calculating predictions...")
  preds = getPredictions(..., newdata=newdata)
  .getDev = function(mu) null$family$dev.resids(y=y, mu=mu, wt = 1)
  kali::DateStamp("Calculating deviances...")
  devs = (nulldev - colSums(apply(preds, 2, FUN=.getDev)))/nulldev 
  names(devs) = names
  return(devs)
}

explDev4 = function(..., newdata, y) {

  null = ..1  
  y = newdata[, y]
  names = grep(as.character(match.call(expand.dots = FALSE)[-1]), 
               pattern="pairlist", value=TRUE)
  names = strsplit(gsub(x=names, pattern=".*\\((.*)\\)", rep="\\1"), split=", ")[[1]]
  
  preds = getPredictions(..., newdata=newdata)

  .getDev = function(mu) null$family$dev.resids(y=y, mu=mu, wt = 1)
  devs = apply(preds, 2, FUN=.getDev) 
  colnames(devs) = names
  return(as.data.frame(devs))
}

calcDev = function(x, null="m0", sort=TRUE) {
  out = colSums(x)
  out = (out[null] - out)/out[null]
  out = out[names(out)!=null]
  out = sort(out, decreasing = TRUE)
  return(out)
}


.getDev = function(mu, y, null) null$family$dev.resids(y=y, mu=mu, wt = 1)

explDev2 = function(models, null=NULL, full=FALSE, sort=TRUE) {
  
  names = names(models)
  if(is.null(names)) names = seq_along(models)
  
  if(inherits(null, "glm")) null = null$null.deviance
  
  .explDev = function(object, null=NULL) {
    if(is.null(null)) null = object$null.deviance 
    (null - object$deviance)/null
  }
  
  .explDev2 = function(object, null=NULL) {
    if(is.null(null)) null = object$null.deviance 
    int = coef(object)["(Intercept)"]
    return(c(int, null, object$deviance, (null - object$deviance)/null))
  }
  
  if(!isTRUE(full)) {
    out = sapply(models, FUN=.explDev, null=null)
    names(out) = names
    if(isTRUE(sort)) out = sort(out, decreasing = TRUE)
  } else {
    out = t(sapply(models, FUN=.explDev2, null=null))
    rownames(out) = names
    colnames(out) = c("Intercept", "null", "deviance", "expl.deviance")
    out = as.data.frame(out)
    if(isTRUE(sort)) out = out[order(out$expl.deviance), ]
  }
  
  return(out)
  
}




mainEffect = function(..., full=FALSE, digits=3) {
  
  names = grep(as.character(match.call(expand.dots = FALSE)[-1]), 
               pattern="pairlist", value=TRUE)
  names = strsplit(gsub(x=names, pattern=".*\\((.*)\\)", rep="\\1"), split=", ")[[1]]
  
  .explDev = function(object, null) {
    x = predict(object, type="terms")
    of = if(is.null(object$offset)) 0 else mean(object$offset)
    y = sum(c(colMeans(x)[1], 0), na.rm=TRUE)
    return(y)
  }
  
  .explDev2 = function(object) {
    y = round(colMeans(predict(object, type="terms")), digits)
  }
  
  if(!isTRUE(full)) {
    out = sapply(list(...), FUN=.explDev)
    names(out) = names
  } else {
    out = lapply(list(...), FUN=.explDev2)

  }
  
  return(out)
  
}



getConstant = function(..., full=FALSE, digits=3) {
  
  names = grep(as.character(match.call(expand.dots = FALSE)[-1]), 
               pattern="pairlist", value=TRUE)
  names = strsplit(gsub(x=names, pattern=".*\\((.*)\\)", rep="\\1"), split=", ")[[1]]
  
  .explDev = function(object, null) {
    x = predict(object, type="terms")
    of = if(is.null(object$offset)) 0 else mean(object$offset)
    y = sum(c(of, attr(x, "constant")), na.rm=TRUE)
    return(y)
  }
  
  .explDev2 = function(object) {
    x = predict(object, type="terms")
    y = c(mean(object$offset), attr(x, "constant"))
  }
  
  if(!isTRUE(full)) {
    out = sapply(list(...), FUN=.explDev)
    names(out) = names
  } else {
    out = lapply(list(...), FUN=.explDev2)
    
  }
  
  return(out)
  
}


getPredictions = function(..., newdata) {
  
  names = grep(as.character(match.call(expand.dots = FALSE)[-1]), 
               pattern="pairlist", value=TRUE)
  names = strsplit(gsub(x=names, pattern=".*\\((.*)\\)", rep="\\1"), split=", ")[[1]]
  
  .predict = function(object, null) {
    y = predict(object, type="response", newdata=newdata)
    return(y)
  }
  
  out = sapply(list(...), FUN=.predict)
  colnames(out) = names
  out = as.data.frame(out)
  
  return(out)
  
}




getDeviance = function(y, mu, family) {
  if (is.character(family))
    family <- eval(parse(text = family))
  if (is.function(family))
    family <- family()
  if (is.null(family$family))
    stop("family not recognized")
  wts = rep(1, length(y))
  out = family$dev.resids(y, mu, wts)
  return(out)
}


.getCol = function(x, n=9) as.numeric(cut(x, seq(floor(min(x)), 
                                                 ceiling(max(x)), length=n+1)))

.getLevels = function(x, n=9) cut(x, seq(floor(min(x)), 
                                         ceiling(max(x)), length=n+1))


.plot = function(fml, data, var, lim, n=4) {
  plot(fml, data=data, pch=".", col=.getCol(data[[var]], n=n), xlim=lim, ylim=lim)
  abline(0,1, lty=3)
  mtext(var, 3, adj=0.05, cex=1.5, line=-2)
  return(invisible())
} 

.boxplot = function(data, by, n=9) {
  
  modvar = grep(x=names(data), pattern="_mod", value=TRUE)
  obsvar = gsub(x=modvar, pattern="_mod", replacement="")
  test = all(c(modvar, obsvar) %in% names(data))
  if(!test) stop("Variables not found.")
  tmp = data.frame(bias = data[[obsvar]] - data[[modvar]])
  tmp$grp = .getLevels(data[[by]], n=n)
  boxplot(bias ~ grp, data=tmp, pars=list(outpch=".", outcol="blue"))
  abline(h=0, lty=3, col="red")
  mtext(by, 1, adj=0.95, cex=1.5, line=-2)
  return(invisible())
  
}

getQBias = function(model, obs, probs, FUN=median) {
  bias = -(obs - model) # observed bias
  q = sapply(bias, FUN = function(x) tapply(x, probs, FUN=FUN))
  qq = matrix(NA, nrow=length(quants), ncol=ncol(q))
  colnames(qq) = colnames(q)
  rownames(qq) = quants
  qq[rownames(q), ] = q
  return(qq)
}

# explainedDeviance -------------------------------------------------------

.explainedDeviance2 = function(object, fixed=NULL, null=NULL, newdata=NULL, obs=NULL, ...) {
  
  cat(date(), "\n")
  
  if(!is.null(newdata) & is.null(obs)) stop("'obs' must be provided to use newdata.")
  
  if(is.null(null)) null = object
  if(inherits(null, "glm")) null = null$null.deviance
  
  term = stats::terms(object$formula, keep.order=TRUE)
  
  sn = attr(term,"term.labels")
  ch = setNames(lapply(sn, FUN=function(...) c(TRUE, FALSE)), sn)
  ch[fixed] = TRUE
  mm = do.call(expand.grid, ch)
  ind = apply(mm, 1, any)
  mm = mm[ind, ]
  
  .update = function(index, term) nfmla = stats::formula(term[index])
  
  .mexplDev = function(object, null) (null - object$deviance)/null
  
  fmlas = apply(mm, 1, FUN=.update, term=term)

  xcall = object$call
  models = list()
  models[[1]] = object
  
  for(i in seq_along(fmlas)[-1]) {
    pb <- txtProgressBar(style=3)
    setTxtProgressBar(pb, i/length(fmlas))
    icall = xcall
    icall$formula = stats::as.formula(fmlas[[i]], env=environment())
    icall$data = quote(object$model)
    models[[i]] = eval(icall)
  }
  
  mm2 = as.data.frame(sapply(mm, as.numeric))
  # if(attr(term,"intercept")==1) mm2$'(Intercept)' = 1 
  
  mm2$dev = sapply(models, FUN=.mexplDev, null=null)
  
  edm = stats::glm(dev ~ . + 0 , data=mm2) 
  
  ed = stats::coef(edm)/sum(stats::coef(edm))  
  
  pb <- txtProgressBar(style=3)
  setTxtProgressBar(pb, 1)
  cat("\n")
  cat(date(), "\n")
  
  return(ed)
  
}




# Statistical distance ----------------------------------------------------

Bhattacharyya = function(x, y) {
  mux = mean(x, na.rm=TRUE)
  muy = mean(y, na.rm=TRUE)
  sd2x = var(x, na.rm=TRUE)
  sd2y = var(y, na.rm=TRUE)
  t1 = sd2x/sd2y + sd2y/sd2x + 2
  t2 = ((mux - muy)^2)/(sd2x + sd2y)
  return(0.25*log(0.25*t1) + 0.25*t2)
}

Db = Bhattacharyya

.mySummary = function(x, sim="sst_mod") {
  n = nrow(x)
  sst = range(x$sst, na.rm=TRUE)
  sst_mod = range(x[, sim], na.rm=TRUE)
  dc = range(x$dc, na.rm=TRUE)
  ds = range(x$ds, na.rm=TRUE)
  ps = sum(x$ds<0)/n
  return(c(n, sst, sst_mod, dc, ds, ps))
}

.getPredFiles = function(model, path, pattern="model") {
  outputFile = gsub(x=basename(model), pattern=pattern, replacement="predH")
  f1 = gsub(x=outputFile, pattern="RData$", replacement="nc4")
  outputFile = gsub(x=basename(model), pattern=pattern, replacement="predF")
  f2 = gsub(x=outputFile, pattern="RData$", replacement="nc4")
  f3 = gsub(x=basename(model), pattern=pattern, replacement="predQ")
  return(c(historical=file.path(path, f1), future=file.path(path, f2), 
           quantile=file.path(path, f3)))
}

addExtremes = function(output, probs=c(0.05, 0.95)) {
  
  probs = unique(round(probs, 2))
  # min
  output$min = list()
  output$min$historical = if(!is.null(output$historical)) apply(output$historical, 1:2, min) else NULL
  output$min$future     = if(!is.null(output$future)) apply(output$future,     1:2, min) else NULL
  output$min$range = suppressWarnings(range(pretty(unlist(output$min))))
  # mean
  output$mean = list()
  output$mean$historical = if(!is.null(output$historical)) apply(output$historical, 1:2, mean) else NULL
  output$mean$future     = if(!is.null(output$future)) apply(output$future,     1:2, mean) else NULL
  output$mean$range = suppressWarnings(range(pretty(unlist(output$mean))))
  # max
  output$max = list()
  output$max$historical = if(!is.null(output$historical)) apply(output$historical, 1:2, max) else NULL
  output$max$future     = if(!is.null(output$future)) apply(output$future,     1:2, max) else NULL
  output$max$range = suppressWarnings(range(pretty(unlist(output$max))))
  
  if(!is.null(probs)) {

    for(i in probs) {
      iName = sprintf("p%02.0f", 100*i)
      output[[iName]] = list()
      output[[iName]]$historical = if(!is.null(output$historical)) apply(output$historical, 1:2, quantile, probs=i, na.rm=TRUE) else NULL
      output[[iName]]$future     = if(!is.null(output$future)) apply(output$future, 1:2, quantile, probs=i, na.rm=TRUE) else NULL
      output[[iName]]$range = suppressWarnings(range(pretty(unlist(output[[iName]]))))
    }
    
  }
  
  class(output) = c("downscaling.results", class(output))
  
  return(output)
}

getPredData = function(model, path, pattern="model") {
  files = .getPredFiles(model, path, pattern=pattern)
  output = list()
  nc = nc_open(files["historical"])
  output$lon = ncvar_get(nc, "lon")
  output$lat = ncvar_get(nc, "lat")
  output$historical = ncvar_get(nc)
  nc_close(nc)
  nc = nc_open(files["future"])
  output$future = ncvar_get(nc)
  nc_close(nc)
  attach(files["quantile"])
  output$quantile = pred
  detach()
  # rm(pred)
  dimnames(output$quantile) = NULL
  
  output = addExtremes(output)

  return(output)
}

getObsData = function(domain) {
  
  inputDir = file.path("ISIMIP2a/InputData/ocean/downscaling", domain, "_input")
  obsFile = file.path(inputDir, sprintf("avhrr-only-v2-%s.198109-200512.nc4", domain))
  nc = nc_open(obsFile)
  on.exit(nc_close(nc))
  x = ncvar_get(nc)
  out = list(lon = ncvar_get(nc, "lon"),
             lat = ncvar_get(nc, "lat"),
             historical=ncvar_get(nc))
  out = addExtremes(out)
  return(out)
  
}


getESMData = function(domain, offset=0, scenario="historical", filter=TRUE, input=NULL,
                      path="ISIMIP2a/InputData/ocean/downscaling", model="ipsl-cm5a-lr",
                      var="to_zs", addHistorical=TRUE) {
  if(is.null(input)) input = ""
  if(scenario!="historical") filter = FALSE
  inputDir = file.path(path, domain, input)
  patt = sprintf("%s_%s_%s_%s_monthly.*\\.nc4", model, scenario, var, domain)
  esmFile = file.path(inputDir, dir(path=inputDir, pattern=patt))
  nc = nc_open(esmFile)
  on.exit(nc_close(nc))
  x = ncvar_get(nc)
  time = ncvar_get(nc, "time")
  dimNames = ncvar_dim(nc)[[1]]
  ind = if(isTRUE(filter)) (time > 12*(1981-1950)+(9-1)) else TRUE
  
  out = list(lon = ncvar_get(nc, dimNames[1]),
             lat = ncvar_get(nc, dimNames[2]))
  out$historical = if(scenario=="historical") x[, , ind] - offset else NULL
  out$future     = if(scenario!="historical") x[, , ind] - offset else NULL
  if(scenario!="historical" & addHistorical) {
    patt = sprintf("%s_%s_%s_%s_monthly.*\\.nc4", model, "historical", var, domain)
    esmFile = file.path(inputDir, dir(path=inputDir, pattern=patt))
    nc0 = nc_open(esmFile)
    on.exit(nc_close(nc0), add = TRUE)
    x = ncvar_get(nc0)
    out$historical = x - offset
  }
  out = addExtremes(out)
  return(out)
}


plot.downscaling.results = function(x, mask=1, zlim=NULL, ...) {
  
  opar = par(no.readonly = TRUE)
  on.exit(par(opar))
  
  mask = as.numeric(mask)
  
  par(mfrow=c(2,3), mar=c(3,3,1,4), oma=c(1,1,1,1))
  zmin  = if(is.null(zlim)) x$min$range else zlim
  zmean = if(is.null(zlim)) x$mean$range else zlim
  zmax  = if(is.null(zlim)) x$max$range else zlim
  image.map(x$lon, x$lat, x$min$historical*mask, zlim=zmin, main="MIN (HIST)")
  image.map(x$lon, x$lat, x$mean$historical*mask, zlim=zmean, main="MEAN (HIST)")
  image.map(x$lon, x$lat, x$max$historical*mask, zlim=zmax, main="MAX (HIST)")
  
  image.map(x$lon, x$lat, x$min$future*mask, zlim=zmin, main="MIN (FUT)")
  image.map(x$lon, x$lat, x$mean$future*mask, zlim=zmean, main="MEAN (FUT)")
  image.map(x$lon, x$lat, x$max$future*mask, zlim=zmax, main="MAX (FUT)")
  
  return(invisible())
  
}

density.downscaling.results = function(x, yobs, ymod, main="", ...) {
  
  opar = par(no.readonly = TRUE)
  on.exit(par(opar))
  
  xd = density(x$quantile)
  xd$y = xd$y/max(xd$y, na.rm=TRUE)
  yobs$y = yobs$y/max(yobs$y, na.rm=TRUE)
  ymod$y = ymod$y/max(ymod$y, na.rm=TRUE)
  
  xrange = range(pretty(c(yobs$x, ymod$x, xd$x)))
  
  plot.new()
  plot.window(xlim=xrange, ylim=c(0,1))
  polygon(yobs, col=kali::makeTransparent("blue", 0.3), border=NA)
  # lines(yobs, col="blue", lwd=2)
  lines(ymod, col="red", lwd=2)
  lines(xd, lwd=2)
  mtext(c("OBSERVED", "ESM", "BIAS CORRECTION"), 3, 
        adj=0.95, line=-(2:4), col=c(4,2,1), font=2)
  mtext(main, 3, adj=0.05, line=-2, cex=1.5)
  axis(1)
  axis(2, las=1)
  box()
  
  return(invisible())
  
}

.mytable = function(x, k) {
  out = numeric(k)
  names(out) = 1:k
  y = table(x)
  if(length(y)>0) {
    out[as.numeric(names(y))] = y
    out = out/sum(out)
  }
  return(out)
}


CTE = function(x, y) {
  x = abs(mean(x, na.rm=TRUE))
  y = abs(mean(y, na.rm=TRUE))
  return((x - y)/(x + y))
}

corplot = function(x, start=1) {
  n = ncol(x)
  image.plot((start-1) + 1:n, (start-1) + 1:n, x, col=biasPalette(),
             xlab="", ylab="")
  grid(nx=n, ny=n)
  return(invisible())
}

removedModels = function(ind, data, maxBias=1) {
  x = data[, ind]
  remove = which(abs(x) > maxBias)
  remModel = data[remove, 1]
  remValue = 100*data[remove, ind]
  msg = sprintf("%s (%0.1f%%)", remModel, remValue)
  msg = paste0(msg, collapse=", ")
  msg = paste0("Removing models: ", msg)
  message(msg)
  return(invisible(remove))
}

readSps = function(files, path) {
  files = file.path(path, files)
  .readSps = function(file) {
    attach(file)
    on.exit(detach())
    if(is.null(mod$sp)) return(NULL)
    return(mod$sp)
  }
  t(sapply(files, FUN=.readSps))
}

loadModel = function(file) {
  attach(file)
  on.exit(detach())
  return(mod)
}



transform_indicators = function(x, g3, max.value=5) {
  
  x$ind3 = pmin(abs(x$ind3), max.value)
  x$ind6 = pmin(abs(x$ind6), max.value)
  x$ind9 = pmin(abs(x$ind9), max.value)
  
  x$ind5 = abs(x$ind5) # 0-10%
  x$ind8 = abs(x$ind8)
  x$ind11 = abs(x$ind11)
  x$ind14 = x$ind14
  x$ind16 = abs(x$ind16)
  x$ind18 = abs(x$ind18)
  x[, g3] = (1 - x[, g3]) # change to "unexplained deviance" 0-40%
  
  return(x)
  
}


calculateBias = function(sim, obs) {
  check = identical(sim$lon, obs$lon) & identical(sim$lat, obs$lat)
  if(!check) stop("spatial dimensions do not match.")
  out = list(lon=sim$lon, lat=sim$lat)
  out$min = sim$min$historical - obs$min$historical
  out$max = sim$max$historical - obs$max$historical
  out$mean = sim$mean$historical - obs$mean$historical
  out$p05 = sim$p05$historical - obs$p05$historical
  out$p95 = sim$p95$historical - obs$p95$historical
  out$range = range(out$min, out$max, out$mean, na.rm=TRUE)
  return(out)
}

calculateTrends = function(x, t.window=120, trendOnly=TRUE) {
  xh = if(!is.null(x$historical)) apply(x$historical, 3, mean, na.rm=TRUE) else NULL
  xf = if(!is.null(x$future)) apply(x$future, 3, mean, na.rm=TRUE) else NULL
  mf = ts(c(xh, xf), freq=12, start=c(1950,1))
  sf = stl(x=mf, s.window = "periodic", t.window = t.window)
  if(isTRUE(trendOnly)) return(sf$time.series[, "trend"])
  return(sf$time.series)
}
