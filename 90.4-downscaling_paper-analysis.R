library(readr)
library(ncdf4)
library(kali)
library(colorful)
library(RColorBrewer)

source("99.9-auxiliar_functions.R")

inputDir = "indicators"

# Reading files -----------------------------------------------------------

files = dir(path=inputDir, pattern=domain)

allInd = paste(sapply(file.path(inputDir, files), FUN=readLines), collapse="\n")
allInd = read.csv(text=allInd, header=FALSE, stringsAsFactors = FALSE)
names(allInd) = c("model", "rep", paste0("ind", seq_len(ncol(allInd)-2)))

modInfo = gsub(x=allInd$model, pattern="m", rep="")
modInfo = do.call(rbind, strsplit(modInfo, split="\\."))
modInfo = as.data.frame(apply(modInfo, 2, as.numeric))
names(modInfo) = c("model", "bs")


y0 = aggregate(allInd[, -(1:2)], by = list(modInfo$model), FUN=mean)
y1 = aggregate(allInd[, -(1:2)], by = list(modInfo$bs), FUN=mean)
yx = aggregate(allInd[, -(1:2)], by = list(allInd$rep), FUN=mean)
y2 = aggregate(allInd[, -(1:2)], by = list(allInd$model), FUN=mean)
y2b = aggregate(allInd[, -(1:2)], by = list(allInd$model), FUN=max)
y2c = aggregate(allInd[, -(1:2)], by = list(allInd$model), FUN=quantile, prob=0.75)

print(y0, digits=2)
print(y1, digits=2)
print(y2, digits=2)

# View(y2[, c(1, 1 + c(3,6,9, 4,7,10))])
# corplot(cor(y2[,-1]))
# corplot(cor(y2[,1 + 1:12]))
# corplot(cor(y2[,1 + 1:12], method = "spearman"))
# corplot(cor(y2[,1 + 13:18]), start=13)
# corplot(cor(y2[,1 + 13:18], method="spearman"), start=13)
# corplot(cor(y2[,1 + 19:27]), start=19)
# corplot(cor(y2[,1 + 19:27], method="spearman"), start=19)

# remove 'bad' models, overfitting

maxBias = 1

# indMean = removedModels("ind3", y2, maxBias = 2)
# indMin  = removedModels("ind6", y2, maxBias = 2)
# indMax  = removedModels("ind9", y2, maxBias = 2)

indMean = removedModels("ind3", y2c, maxBias = 1)
indMin  = removedModels("ind6", y2c, maxBias = 1)
indMax  = removedModels("ind9", y2c, maxBias = 1)

overfit = y2[unique(c(indMean, indMin, indMax)), 1]
overfit = "m0.0"

checkInd = y2[!(y2$Group.1 %in% overfit), ]
indicators = allInd[!(allInd$model %in% overfit), ]

# indicator 3 and 4: normalize by value of null model

i1 = c(6, 3, 9) #c(6,3,9,11,5,8)
i2 = c(13, 14, 18) #c(13:14,16,18)
i3 = c(19, 20, 22, 27) #c(19:22,24,25,27)

indLabels = c(paste0("A", seq_along(i1)),
              paste0("B", seq_along(i2)),
              paste0("C", seq_along(i3)))

chosenOnes = paste0("ind", c(i1, i2, i3))

g1 = 0 + seq_along(i1)        # FUTURE BIAS
g2 = max(g1) + seq_along(i2)  # RESIDUALS
g3 = max(g2) + seq_along(i3)  # DEVIANCE

included = TRUE 
x = indicators[included, chosenOnes]

x$ind3 = pmin(abs(x$ind3), 5)
x$ind6 = pmin(abs(x$ind6), 5)
x$ind9 = pmin(abs(x$ind9), 5)

# x$ind5 = abs(x$ind5) # 0-10%
# x$ind8 = abs(x$ind8)
# x$ind11 = abs(x$ind11)
x$ind14 = pmin(x$ind14, 20)
# x$ind16 = abs(x$ind16)
x$ind18 = abs(x$ind18)

x[, g3] = pmax(0, unlist(x[, g3])) # remove negative deviance
x[, g3] = (1 - x[, g3]) # change to "unexplained deviance" 0-40%

x$ind27 = pmin(x$ind27, 10)

nx = rbind(apply(x, 2, quantile, prob=0.05),
           apply(x, 2, median),
           apply(x, 2, quantile, prob=0.95))

mx0 = aggregate(x, by = list(indicators$model[included]), FUN=quantile, prob=0.95)
# mx = aggregate(x, by = list(indicators$model[-1]), FUN=mean)
names(mx0)[1] = "model"

mind = TRUE 
mx = mx0[mind, ]

labs = gsub(names(mx[,-1]), pattern="ind", replacement="")

# ylim0 = c(rep(5,3), 0.5, 12, 1, rep(0.4, length(g3)))
# ylim0 = rbind(0, ylim0)

ylim1 = c(rep(5,3), 1, 20, 1, rep(1, length(g3)))
ylim1 = rbind(0, ylim1)
ylim0 = ylim1
ylim = if(domain=="biscay-celtic") ylim0 else ylim1

theta1 = 60/length(g1) + head(seq(0, 120, length=length(g1)+1), length(g1))
theta2 = 60/length(g2) + head(seq(0, 120, length=length(g2)+1), length(g2))
theta3 = 60/length(g3) + head(seq(0, 120, length=length(g3)+1), length(g3))

theta = 30 + c(theta1, 120 + theta2, 240 + theta3)

cols = c("black", "red", "darkgreen", "blue", "cyan4", "magenta", "orange")
cols = brewer.pal(8, name = "Dark2")
mfrow = c(4,5)
par(mfrow=mfrow, mar=c(0,0,1,0), oma=c(1,1,1,1))
for(i in seq_len(nrow(mx))) {
  icol = cols[1+as.numeric(gsub(mx[i,1], pattern="m([0-9].*)\\..*", replacement="\\1"))]
  px = as.numeric(mx[i,-1])
  px1 = px2 = px3 = numeric(length(px))
  px1[g1] = px[g1]
  px2[g2] = px[g2]
  px3[g3] = px[g3]
  
  labs = if(i==1) indLabels else NA 
  spiderplot(px, ylim, fill=NA, type="n", alpha=0.5,
             labels=labs, at=c(0.25, 0.50, 0.75), theta=theta, rmin=0)
  spider(px1, ylim, fill=icol, type="n", alpha=0.75,
         labels=labs, at=c(0.25, 0.50, 0.75), theta=theta, rmin=0)
  spider(px2, ylim, fill=icol, type="n", alpha=0.75,
         labels=labs, at=c(0.25, 0.50, 0.75), theta=theta, rmin=0)
  spider(px3, ylim, fill=icol, type="n", alpha=0.75,
         labels=labs, at=c(0.25, 0.50, 0.75), theta=theta, rmin=0)
  mtext(mx[i,1], 3, line=-0.5, col=icol)
  segments(x0=0, x1=0, y0=0, y1=-1, lty=2)
  segments(x0=0, x1=cos(pi/6), y0=0, y1=sin(pi/6), lty=2)
  segments(x0=0, x1=cos(pi-pi/6), y0=0, y1=sin(pi-pi/6), lty=2)
}
# spiderplot(px, ylim, fill=NA, type="n", alpha=0.5,
           # labels=indLabels, at=c(0.25, 0.50, 0.75), theta=theta, rmin=0)
