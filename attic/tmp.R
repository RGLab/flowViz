## load a flowFrame
library(flowViz)
data(GvHD)
fcs1 <- GvHD[[1]]



############################################################################
##                            rectangelGates first                        ##
############################################################################

## This should cause a warning because we need to guess the
## plotting parameters
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650),
              "SSC-H" = c(400, 575))
plot(fcs1, c("FSC-H", "SSC-H"))
glines(rg)
glines(rg, col=4, verbose=FALSE)
gpolygon(rg, col=2)
gpolygon(rg, col=3, verbose=FALSE)
gpoints(rg, fcs1)
gpoints(rg, fcs1, verbose=FALSE, pch=20)

## This should fail because we don't know how to match to the plot
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650))
glines(rg)
gpolygon(rg, col=2)
gpoints(rg, fcs1)
## This should fail because we need a character of lengh 2 as second argument
glines(rg, "a")
gpolygon(rg, "a")
gpoints(rg, fcs1, "a")
## This should fail because the parameters don't match
glines(rg, c("a", "b"))
gpolygon(rg, c("a", "b"))
gpoints(rg, fcs1, c("a", "b"))


## This should create a region gate in the y dimension and also know
## how to deal with infinite values
plot(fcs1, c("FSC-H", "SSC-H"))
rg <- rectangleGate(filterId="Rectangle", "SSC-H" = c(300, 600))
glines(rg, c("FSC-H", "SSC-H"))
glines(rg, c("FSC-H", "SSC-H"), verbose=FALSE, col=4)
gpolygon(rg, c("FSC-H", "SSC-H"), col=2)
gpolygon(rg, c("FSC-H", "SSC-H"), col=3, verbose=FALSE)
gpoints(rg, fcs1, c("FSC-H", "SSC-H"))
gpoints(rg, fcs1, c("FSC-H", "SSC-H"), verbose=FALSE, pch=20)
plot(fcs1, c("FSC-H", "SSC-H"))
rg <- rectangleGate(filterId="Rectangle", "SSC-H" = c(-Inf, 650))
glines(rg, c("FSC-H", "SSC-H"))
glines(rg, c("FSC-H", "SSC-H"), verbose=FALSE, col=4)
gpolygon(rg, c("FSC-H", "SSC-H"), col=2)
gpolygon(rg, c("FSC-H", "SSC-H"), col=3, verbose=FALSE)
gpoints(rg, fcs1, c("FSC-H", "SSC-H"))
gpoints(rg, fcs1, c("FSC-H", "SSC-H"), verbose=FALSE, pch=20)

## This should create a region gate in the y dimension and also know
## how to deal with infinite values
plot(fcs1, c("FSC-H", "SSC-H"))
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(300, 600))
glines(rg, c("FSC-H", "SSC-H"))
glines(rg, c("FSC-H", "SSC-H"), verbose=FALSE, col=4)
gpolygon(rg, c("FSC-H", "SSC-H"), col=2)
gpolygon(rg, c("FSC-H", "SSC-H"), col=3, verbose=FALSE)
gpoints(rg, fcs1, c("FSC-H", "SSC-H"))
gpoints(rg, fcs1, c("FSC-H", "SSC-H"), verbose=FALSE, pch=20)
plot(fcs1, c("FSC-H", "SSC-H"))
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(-Inf, 250))
glines(rg, c("FSC-H", "SSC-H"))
glines(rg, c("FSC-H", "SSC-H"), verbose=FALSE, col=4)
gpolygon(rg, c("FSC-H", "SSC-H"), col=2)
gpolygon(rg, c("FSC-H", "SSC-H"), col=3, verbose=FALSE)
gpoints(rg, fcs1, c("FSC-H", "SSC-H"))
gpoints(rg, fcs1, c("FSC-H", "SSC-H"), verbose=FALSE, pch=20)


## This should create a warning that we don't need a filter result
## but still work
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650),
              "SSC-H" = c(400, 575))
res <- filter(fcs1, rg)
plot(fcs1, c("FSC-H", "SSC-H"))
glines(rg, res)
glines(rg, res, col=3, verbose=FALSE)
gpolygon(rg, res, col=2)
gpolygon(rg, res, col=4, verbose=FALSE)
gpoints(rg, fcs1, filterResult=res)
gpoints(rg, fcs1, filterResult=res, verbose=FALSE, pch=20)

## This should work with the usual warnings
plot(fcs1, c("FSC-H", "SSC-H"))
glines(res)
glines(res, col=3, verbose=FALSE)
gpolygon(res, col=2)
gpolygon(res, col=4, verbose=FALSE)
gpoints(res, fcs1)
gpoints(res, fcs1, verbose=FALSE, pch=20)


############################################################################
##                            now polygonGates                            ##
############################################################################

## This should cause a warning because we need to guess the
## plotting parameters
sqrcut <- matrix(c(300,400,600,500, 200, 50, 100, 300, 400,70),
                 ncol=2)
colnames(sqrcut) <- c("FSC-H","SSC-H")
pg <- polygonGate(filterId="nonDebris", boundaries= sqrcut)
plot(fcs1, c("FSC-H", "SSC-H"))
glines(pg)
glines(pg, col=3, verbose=FALSE)
gpolygon(pg, col=2)
gpolygon(pg, col=4, verbose=FALSE)
gpoints(pg, fcs1)
gpoints(pg, fcs1, verbose=FALSE, pch=20)


## This should create a warning that we don't need a filter result
## but still work
res <- filter(fcs1, pg)
plot(fcs1, c("FSC-H", "SSC-H"))
glines(pg, res)
glines(pg, res, verbose=FALSE, col=3)
gpolygon(pg, res, col=2)
gpolygon(pg, res, col=4, verbose=FALSE)


## This should work with the usual warnings
plot(fcs1, c("FSC-H", "SSC-H"))
glines(res)
glines(res, col=3, verbose=FALSE)
gpolygon(res, col=2)
gpolygon(res, col=4, verbose=FALSE)
gpoints(res, fcs1)
gpoints(res, fcs1, pch=20, verbose=FALSE)


############################################################################
##                               norm2Filter                              ##
############################################################################

## This should fail because we need either the raw data or a filter result
## to plot norm2Filters
nf <- norm2Filter(filterId = "BVNorm", "FSC-H", "SSC-H", scale=2)
plot(fcs1, c("FSC-H", "SSC-H"))
glines(nf)
glines(nf, c("FSC-H", "SSC-H"))
gpolygon(nf, col=2)
gpolygon(nf, c("FSC-H", "SSC-H"), col=2)


## This should cause a warning because we need to guess the
## plotting parameters
res <- filter(fcs1, nf)
glines(nf, res)
glines(nf, res, verbose=FALSE, col=3)
gpolygon(nf, res, col=2)
gpolygon(nf, res, col=4, verbose=FALSE)
gpoints(nf, fcs1, filterResult=res)
gpoints(nf, fcs1, filterResult=res, verbose=FALSE, pch=20)
plot(fcs1, c("FSC-H", "SSC-H"))
glines(nf, fcs1)
glines(nf, fcs1, col=3, verbose=FALSE)
gpolygon(nf, fcs1, col=2)
gpolygon(nf, fcs1, col=4, verbose=FALSE)
gpoints(nf, fcs1)
gpoints(nf, fcs1, verbose=FALSE, pch=20)



## This should work with the usual warnings
plot(fcs1, c("FSC-H", "SSC-H"))
glines(res)
glines(res, verbose=FALSE, col=3)
gpolygon(res, col=2)
gpolygon(res, col=4, verbose=FALSE)
gpoints(res, fcs1)
gpoints(res, fcs1, verbose=FALSE, pch=20)



############################################################################
##                               curv2Filter                              ##
############################################################################
## This should fail because we need either the raw data or a filter result
## to plot norm2Filters
cf <- curv2Filter(filterId = "BVCurv", "FSC-H", "SSC-H")
plot(fcs1, c("FSC-H", "SSC-H"))
glines(cf)
glines(nf, c("FSC-H", "SSC-H"))
gpolygon(nf, col=2)
gpolygon(nf, c("FSC-H", "SSC-H"), col=2)


## This should cause a warning because we need to guess the
## plotting parameters
res <- filter(fcs1, cf)
glines(cf, fcs1)
glines(cf, fcs1, col=3, verbose=FALSE)
gpolygon(cf, fcs1, col=2)
gpolygon(cf, fcs1, col=4, verbose=FALSE)
gpoints(cf, fcs1)
gpoints(cf, fcs1, verbose=FALSE, pch=20)
plot(fcs1, c("FSC-H", "SSC-H"))
glines(cf, res)
glines(cf, res, col=3, verbose=FALSE)
gpolygon(cf, res)
gpolygon(cf, res, col=4, verbose=FALSE)
gpoints(cf, fcs1, filterResult=res)
gpoints(cf, fcs1, filterResult=res, verbose=FALSE, pch=20, col=2)


## This should work with the usual warnings
plot(fcs1, c("FSC-H", "SSC-H"))
glines(res)
glines(res, col=3, verbose=FALSE)
gpolygon(res, col=2)
gpolygon(res, col=4, verbose=FALSE)
gpoints(res, fcs1)
gpoints(res, fcs1, verbose=FALSE, pch=20, col=2)



############################################################################
##                               curv1Filter                              ##
############################################################################
## This should fail because we need either the raw data or a filter result
## to plot norm2Filters
c1f <- curv1Filter(filterId = "BVCurv", "SSC-H")
res <- filter(fcs1, c1f)
plot(fcs1, c("FSC-H", "SSC-H"))
glines(c1f, c("FSC-H", "SSC-H"))
gpolygon(c1f, c("FSC-H", "SSC-H"))

## This should fail because we can't match the plotting parameters
glines(kf, fcs1)
gpolygon(c1f, fcs1)

## This should work
glines(c1f, fcs1, channels=c("FSC-H", "SSC-H"))
glines(res, channels=c("FSC-H", "SSC-H"), col=4)
gpolygon(c1f, fcs1, channels=c("FSC-H", "SSC-H"))
gpolygon(res, channels=c("FSC-H", "SSC-H"), col=4)
gpoints(res, fcs1, c("FSC-H", "SSC-H"), pch=20)


############################################################################
##                               kmeansFilter                             ##
############################################################################
## All these should fail because we don't know how to plot lines or
## polygons for kmeansFilter
kf <- kmeansFilter("kmfilt", "FSC-H" = c("Low", "High"))
res <- filter(fcs1, kf)
plot(fcs1, c("FSC-H", "SSC-H"))
glines(kf, c("FSC-H", "SSC-H"))
glines(kf, c("FSC-H", "SSC-H"), verbose=FALSE)
gpolygon(kf, c("FSC-H", "SSC-H"))
gpolygon(kf, c("FSC-H", "SSC-H"), verbose=FALSE)

## This should fail because we can't match the plotting parameters
gpoints(kf, fcs1)
gpoints(kf, fcs1, verbose=FALSE, pch=20)

## This should work
gpoints(kf, fcs1, c("FSC-H", "SSC-H"))
gpoints(kf, fcs1, c("FSC-H", "SSC-H"), verbose=FALSE, pch=20, col=2)
plot(fcs1, c("FSC-H", "SSC-H"))
gpoints(kf, fcs1, c("FSC-H", "SSC-H"), filterResult=res)
gpoints(kf, fcs1, c("FSC-H", "SSC-H"), filterResult=res, verbose=FALSE,
        pch=20, col=2)






## need to fix ellipsoidGates first
## m <- matrix(c(300, 200, 500, 200),ncol = 2)
## colnames(m) <- c("FSC-H", "SSC-H")
## eg <- ellipsoidGate(filterId= "ellipsoidGateIII", m, distance=400)
## plot(fcs1, c("FSC-H", "SSC-H"))
## lines(pe)
## polygon(eg, col=2)




