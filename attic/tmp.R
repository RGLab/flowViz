## load a flowFrame
library(flowViz)
library(ncdfFlow)
library(grid)
library(IDPmisc)
library(hexbin)
library(latticeExtra)
library(gridExtra)
lapply(list.files("~/rglab/workspace/flowViz/R",full=T),source)
data(GvHD)
fcs1 <- GvHD[[1]]
fcs2 <- GvHD[[1]]
colnames(fcs2)[1] <- "noname"
fs <- GvHD[c(1,2,8,9)]
fs <- ncdfFlowSet(fs)
g1 <- rectangleGate("SSC-H"=c(400,Inf))
g2 <- rectangleGate("SSC-H"=c(900,950))
overlay <- Subset(fs[1],g2)
dev.off()
x11()
grid.arrange(
    xyplot(`FSC-H`~`SSC-H`,fs[1]
              ,xbin=64
              ,smooth=F 
            ,filter = sapply(sampleNames(fs),function(x)filters(list(g1,g2)),simplify=F)
#              ,filter = g1
#    , overlay =overlay
    ,stats=T
      )
  ,xyplot(`FSC-H`~`SSC-H`, fs[1]
      ,xbin=64
                ,smooth=F 
      , par.settings = list(overlay.symbol = list(fill = "black")
      )
      ,filter = g1, overlay =overlay
      )
  )
  
aa <- trellis.par.get()
bb <- flowViz.par.get()

system.time(
densityplot(~`SSC-H`
                ,fs
              ,stack=T
#              ,filter = sapply(sampleNames(fs),function(x)filters(list(g1,g2)),simplify=F)
              ,filter =g1
              ,fitGate = F
      ,stats =T
          )

)
pData(GvHD)
## a wrapper that catches and reports errors
wrap <- function(x)
{
    tmp <- try(x)
    if(!is(tmp, "try-error")){
        print("done")
        if(!is.null(tmp))
            warning("Return value should be NULL.", call.=FALSE)
    }
}
    

## plot a gate with all possible combinations of arguments
plotGates <- function(gate, data=fcs1, verbose=TRUE)
{
    opar <- par(ask=TRUE)
    oo <- options(warn=1)
    on.exit({par(opar); options(oo)})
    yp <- xp <- p <- parameters(gate)
    xp[1] <- "falseX"
    yp[2] <- "falseY"
    fres <- try(filter(data, gate))
    ## first without a second argument
    print("Only dispatch on filter")
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(gate, verbose=verbose, lwd=2, lty=3))
    ## Now with optional channel argument
    print("Only dispatch on filter, but with optional channels argument")
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(gate, channels=c("FSC-H", "SSC-H"), lty=3,
             lwd=2, verbose=verbose))
    ## Now the second argument are the channels
    print("Dispatch on filter and channel names")   
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(gate, c("FSC-H", "SSC-H"), verbose=verbose,
                col=2, lwd=2, lty=3))
    ## Some of the channels are not in the filter
    print("Dispatch on filter and channel names with y channel wrong")  
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(gate, yp, verbose=verbose, col=2, lwd=2, lty=3))
    ## And with optional channel argument added on top
    print("Dispatch on filter with optional channels argument")  
    plot(fcs1, c("SSC-H", "FSC-H"))
    wrap(fun(gate, c("FSC-H", "SSC-H"), verbose=verbose,
                col=2, lwd=2, lty=3, channels=c("SSC-H", "FSC-H")))
    ## Some of the channels are not in the filter
    print("Dispatch on filter and channel names with y channel wrong and with optional channels argument")
    plot(fcs1, c("SSC-H", "FSC-H"))
    wrap(fun(gate, yp, verbose=verbose, col=2, lwd=2, lty=3, channels=xp))
    ## Now the second argument is a filterResult
    print("Dispatch on filter and filterResult")   
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(gate, fres, verbose=verbose,
                lwd=2, lty=3))
    ## And with optional channel argument added on top
    print("Dispatch on filter and filterResult with optional channels argument")   
    plot(fcs1, c("SSC-H", "FSC-H"))
    wrap(fun(gate, fres, verbose=verbose,
                col=5, lwd=2, lty=3, channels=c("SSC-H", "FSC-H")))
    ## Now the second argument is a flowFrame
    print("Dispatch on filter and flowFrame")   
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(gate, data, verbose=verbose,
                col=5, lwd=2, lty=3))
    ## And with optional channel argument added on top
    print("Dispatch on filter and flowFrame with optional channels argument")   
    plot(fcs1, c("SSC-H", "FSC-H"))
    wrap(fun(gate, data, verbose=verbose,
                lwd=2, lty=3, channels=c("SSC-H", "FSC-H")))
}


## plot points in a gate with all possible combinations of arguments
plotPoints <- function(gate, data=fcs1, verbose=TRUE)
{
    opar <- par(ask=TRUE)
    oo <- options(warn=1)
    on.exit({par(opar); options(oo)})
    yp <- xp <- p <- parameters(gate)
    xp[1] <- "falseX"
    yp[2] <- "falseY"
    fres <- try(filter(data, gate))
    ## first without a third argument
    print("Only dispatch on filter")
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(gate, fcs1, verbose=verbose, lwd=2, lty=3))
    ## Now with the channel argument
    print("Only dispatch on filter, but with optional channels argument")
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(gate, fcs1, channels=c("FSC-H", "SSC-H"), lty=3,
             lwd=2, verbose=verbose))
    ## Some of the channels are not in the filter
    print("Dispatch on filter and channel names with y channel wrong")  
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(gate, fcs1, yp, verbose=verbose, col=2, lwd=2, lty=3))
    ## Now the first argument is a filterResult
    print("Dispatch on filter and filterResult")   
    plot(fcs1, c("FSC-H", "SSC-H"))
    wrap(fun(fres, fcs1, verbose=verbose,
                lwd=2, lty=3))
    ## And with channel argument added on top
    print("Dispatch on filter and filterResult with channels argument")   
    plot(fcs1, c("SSC-H", "FSC-H"))
    wrap(fun(fres, fcs1, verbose=verbose,
                col=5, lwd=2, lty=3, channels=c("SSC-H", "FSC-H")))
}



############################################################################
##                            rectangelGates first                        ##
############################################################################
## lines
fun <- glines
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650),
              "SSC-H" = c(400, 575))
plotGates(rg)
plotGates(rg, verbose=FALSE)
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650),
              "SSC-H" = c(400, 575), "FL2-H" = c(400, 500))
plotGates(rg)
plotGates(rg, verbose=FALSE)
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650))
plotGates(rg)
plotGates(rg, verbose=FALSE)

## polygons
fun <- gpolygon
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650),
              "SSC-H" = c(400, 575))
plotGates(rg)
plotGates(rg, verbose=FALSE)
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650),
              "SSC-H" = c(400, 575), "FL2-H" = c(400, 500))
plotGates(rg)
plotGates(rg, verbose=FALSE)
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650))
plotGates(rg)
plotGates(rg, verbose=FALSE)

## points
fun <- gpoints
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650),
              "SSC-H" = c(400, 575))
plotPoints(rg)
plotPoints(rg, verbose=FALSE)
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650),
              "SSC-H" = c(400, 575), "FL2-H" = c(400, 500))
plotPoints(rg)
plotPoints(rg, verbose=FALSE)
rg <- rectangleGate(filterId="Rectangle", "FSC-H" = c(250, 650))
plotPoints(rg)
plotPoints(rg, verbose=FALSE)


############################################################################
##                                quadGates                               ##
############################################################################
## lines
fun <- glines
qg <- quadGate(filterId="nonDebris", "FSC-H"=500, "SSC-H"=800)
plotGates(qg)
plotGates(qg, verbose=FALSE)

## polygons
fun <- gpolygon
plotGates(qg)
plotGates(qg, verbose=FALSE)

## points
fun <- gpoints
plotPoints(qg)
plotPoints(qg, verbose=FALSE)


############################################################################
##                               polygonGates                             ##
############################################################################
## lines
fun <- glines
sqrcut <- matrix(c(300,400,600,500, 200, 50, 100, 300, 400,70),
                 ncol=2)
colnames(sqrcut) <- c("FSC-H","SSC-H")
pg <- polygonGate(filterId="nonDebris", boundaries= sqrcut)
plotGates(pg)
plotGates(pg, verbose=FALSE)

## polygons
fun <- gpolygon
plotGates(pg)
plotGates(pg, verbose=FALSE)

## points
fun <- gpoints
plotPoints(pg)
plotPoints(pg, verbose=FALSE)


############################################################################
##                               norm2Filter                              ##
############################################################################
## lines
fun <- glines
nf <- norm2Filter(filterId = "BVNorm", "FSC-H", "SSC-H", scale=2)
plotGates(nf)
plotGates(nf, verbose=FALSE)

## polygons
fun <- gpolygon
plotGates(nf)
plotGates(nf, verbose=FALSE)

## points
fun <- gpoints
plotPoints(nf)
plotPoints(nf, verbose=FALSE)


############################################################################
##                               curv2Filter                              ##
############################################################################
## lines
fun <- glines
c2f <- curv2Filter(filterId = "BVCurv", "FSC-H", "SSC-H")
plotGates(c2f)
plotGates(c2f, verbose=FALSE)

## polygons
fun <- gpolygon
plotGates(c2f)
plotGates(c2f, verbose=FALSE)

## points
fun <- gpoints
plotPoints(c2f)
plotPoints(c2f, verbose=FALSE)


############################################################################
##                               curv1Filter                              ##
############################################################################
## lines
fun <- glines
c1f <- curv1Filter(filterId = "BVCurv", "SSC-H")
plotGates(c1f)
plotGates(c1f, verbose=FALSE)

## polygons
fun <- gpolygon
plotGates(c1f)
plotGates(c1f, verbose=FALSE)

## points
fun <- gpoints
plotPoints(c1f)
plotPoints(c1f, verbose=FALSE)


############################################################################
##                               kmeansFilter                             ##
############################################################################
## lines
fun <- glines
kf <- kmeansFilter("kmfilt", "FSC-H" = c("Low", "High"))
plotGates(kf)

## polygons
fun <- gpolygon
plotGates(kf)

## points
fun <- gpoints
plotPoints(kf)
plotPoints(kf, verbose=FALSE)

############################################################################
##                               hexbin
############################################################################
xyplot(`FSC-H` ~ `SSC-H`, GvHD[1:3], smooth=F,xbin=128)

############################################################################
##                               conditioning lattice
############################################################################

library(flowViz)
library(IDPmisc)
data(GvHD)

fs<-GvHD[c(1,2,9,10)]

xyplot(`SSC-H` ~ `FSC-H`|Patient:Visit:name ,data =fs)
xyplot(`SSC-H` ~ `FSC-H`|Patient:name ,data =fs)

xyplot(`SSC-H` ~ `FSC-H` ,data =fs)
exprs(fs[[3]])<-exprs(fs[[3]])[0,,drop=F]
xyplot(`SSC-H` ~ `FSC-H` ,data =fs[[3]])

xyplot(Grade~factor(name)|Patient+Visit,data=pData(fs))

require(flowWorkspace)
library(flowStats)
library(flowClust)
lapply(list.files("~/rglab/workspace/flowViz/R",full=T),source)

load("~/rglab/workspace/HIMCLyoplate/Gottardo/flowCAP/wf.rda")
png("~/rglab/workspace/HIMCLyoplate/Gottardo/pipeline/analysis/Tcell/flowVizFixed.png",width=800,height=600)
plotGate(wf)
dev.off()



overlay(wf[1:2],gate=4,overlay.gate.indices=c(4,5))

overlay <- function(gs,gate=9,overlay.gate.indices=c(9,40),trans=TRUE,grid=65,h=c(0.1,0.1),nlines=25,wh=1,...){
	#gs is a gating set, and we want to generate one plot at a time..
	browser()
	gh <- gs[[wh]]
	p <- getParent(gh,gate)
	d <- getData(gh,p)
	g <- getGate(gh,gate)
	cn <- colnames(g@boundaries)
	d <- exprs(d)[,cn]
	
	x <- d[,1]
	y <- d[,2]
	
	ogi <- and(do.call(cbind,lapply(overlay.gate.indices,function(i)flowWorkspace:::getIndices(gh,i))))
	a <- and(do.call(cbind,lapply(c(overlay.gate.indices,gate),function(i)flowWorkspace:::getIndices(gh,i))))
	b<-flowWorkspace:::getIndices(gh,p)
	prop <- prop.table(table(a[b]))[2]
	cl <- kde2d(x=x,y=y,h=h,n=grid)
	if(trans){
		cl$z <- sqrt(cl$z)
	}
	contour(cl,n=nlines,...)
	d <- exprs(getData(gh)[,cn])[ogi,,drop=FALSE]
	points(x=d[,1],y=d[,2],col="red",cex=3,pch='.')
	
	polygon(getGate(gh,gate)@boundaries,border="blue")
	text(3.5,2,signif(prop,3))
}
or <- function(x,...){
	if(ncol(x)==1){
		x
	}else if(ncol(x)==2){
		x[,1]|x[,2]
	}else{
		Recall(x[,-1L])
	}
}
and <- function(x,...){
	if(ncol(x)==1){
		x
	}else if(ncol(x)==2){
		x[,1]&x[,2]
	}else{
		Recall(x[,-1L])
	}
}
