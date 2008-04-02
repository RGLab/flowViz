## Plot values of a flowSet against the time domain.
## TO DO: Eventualy, we might want individual methods for flowSets and for
## flowFrames, although it might be worthwhile to reimplement using
## lattice graphics first...


## run over a cytoFrame, bin the values according to the time domain
## and compute a location measure locM as well as variances varM for each bin.
## The result of this function will be the input to the plotting functions
## and the basis for the quality score.
prepareSet <- function(x, parm, binSize, locM=median, varM=mad){
    cn <- colnames(x)
    wcn <- grep("Time|time", cn)
    if(length(wcn)==0)
        stop("No time domain recording for this data")
    else
        ti <- wcn[1]
    xx <- x[, ti]
    ord <- order(xx)
    xx <- xx[ord]
    yy <- x[ord, parm]
    lenx <- length(xx)
    nrBins <- floor(lenx/binSize)
    ## time parameter is already binned or very sparse events
    if(length(unique(xx)) < nrBins){
        tmpy <- split(yy, xx)
        yy <- sapply(tmpy, locM, na.rm=TRUE)
        xx <- unique(xx)
        var <- sapply(tmpy, varM, na.rm=TRUE)
        binSize <- 1
    }else{
        ## bin values in nrBins bins
        if(lenx > binSize){
            cf <- c(rep(1:nrBins, each=binSize),
                    rep(nrBins+1, lenx-nrBins*binSize))
            stopifnot(length(cf) == lenx)
            tmpx <- split(xx,cf)
            tmpy <- split(yy,cf)
            yy <- sapply(tmpy, locM, na.rm=TRUE)
            xx <- sapply(tmpx, mean, na.rm=TRUE)
            var <- sapply(tmpy, varM, na.rm=TRUE)
        }else{
            ## very little events
            warning("Low number of events", call.=FALSE)
            tmpy <- split(yy,xx)
            yy <- sapply(tmpy, locM, na.rm=TRUE)
            xx <- unique(xx)
            var <- sapply(tmpy, varM, na.rm=TRUE)
            binSize <- 1
        }
    }
    return(list(smooth=cbind(xx,yy), variance=var, binSize=binSize))
}



## Call above function for each flowFrame in a flowSet and produce the
## plots, according to "type" (see below for details). Also compute a
## quality value modified by the setting of varCut, which is basically
## the sum of average positive distances of means for each bin from the global
## mean confidence interval, scaled by the variance cutoff:
##       sum(z[z>0])/length(z)/varCut
timeLinePlot <- function(x, channel, type=c("stacked", "scaled", "native"),
                         col, ylab=names(x), binSize, varCut=1, ...)
{
    if(!length(channel)==1)
        stop("'channel' must be character scalar")
    if(!channel %in% colnames(x[[1]]))
        stop(channel, " is not a valid channel in this flowSet.")
    ## making sure the sample names fit the axis annotation
    sampleNames(x) <- truncNames(sampleNames(x))
    ## lets fix ourselves some nioce colors
    if(missing(col)){
        require(RColorBrewer)
        colp <- brewer.pal(8, "Dark2")
        col <- colorRampPalette(colp)(length(x))
        set.seed(1000)
        col <- sample(col) 
    }else{
        if(length(col)!=1 || length(col)!=length(x))
            stop("'col' must be color vector of length 1 or same length ",
                 "as the flowSet")
    }
    ## a reasonable default for the bin size
    if(missing(binSize))
        binSize <- min(max(1, floor(median(fsApply(x, nrow)/100))), 500)
    opar <- par(c("mar", "mgp", "mfcol", "mfrow", "las"))
    on.exit(par(opar))
    ## bin the data and compute location and variance
    timeData <-  fsApply(x, prepareSet, channel, binSize=binSize,
                         use.exprs=TRUE, simplify=FALSE)
    type <- match.arg(type)
    mr <- range(x[[1]])[,channel]
    mr[1] <- max(mr[1], 0)
    ## standardize to compute maeningful scale-free QA scores
    med <- sapply(timeData, function(z) median(z$smooth[,2], na.rm=TRUE))
    gvars <- sapply(timeData, function(x) mean(x$variance))
    stand <-  mapply(function(z,m,v) abs(z$smooth[,2]-m)/(v*varCut), timeData,
                 med, gvars)
    ## create teh plot
    if(length(med)==1){
        nativePlot(timeData, p=channel, range=mr, col="darkblue", med=med,
                   varCut=varCut, ...)
        
        return((sum(stand[stand>0])/length(stand))/varCut)
    }
    layout(matrix(1:2), heights=c(0.8, 0.2))
    switch(type,
           "scaled"=scaledPlot(timeData, p=channel, range=mr, col=col,
           med=med, varCut=varCut, ...),
           "stacked"=stackedPlot(timeData, p=channel, range=mr, col=col,
           ylab=ylab, med=med, varCut=varCut, ...),
           "native"=nativePlot(timeData, p=channel, range=mr, col=col,
           med=med, varCut=varCut, ...),
           stop("Unknown type"))
    qaScore <- sapply(stand, function(z) sum(z[abs(z)>1])/length(z))*100
    par(mar=c(5,3,0,3), las=2)
    on.exit(par(opar))
    top <- 2
    barplot(qaScore, axes=FALSE, col=col, cex.names=0.7,
            ylim=c(0, min(c(top, max(qaScore)))), border=col, space=0.2)
    wh <- which(qaScore >= top)
    points((wh+wh*0.2)-0.5, rep(top-(top/12), length(wh)), pch=17, col="white",
           cex=0.7)
    attr(qaScore, "binSize") <- binSize
    return(qaScore)
}



## truncData
truncNames <- function(names){
    nc <- nchar(names)
    names[nc>11] <- paste(substr(names[nc>11], 1, 8), "...", sep="")
    return(names)
}

## lighten up colors
desat <- function(col, by=50)
{
    rgbcol <- col2rgb(col)
    rgbcol <- rgbcol + by
    rgbcol[rgbcol>255] <- 255
    return(rgb(rgbcol[1,], rgbcol[2,], rgbcol[3,], max=255))
}



## align values around 0 and plot
scaledPlot <- function(y, p, main=paste("time line for", p),
                       range, col, med, lwd, varCut, ...){
    par(mar=c(1,2.5,3,2.5), mgp=c(1.5,0.5,0))
    yy <- mapply(function(z, m) data.frame(x=z$smooth[,1], y=z$smooth[,2]-m),
                 y, med, SIMPLIFY=FALSE)
    var <- sapply(y, function(z) mean(z$var))
    maxX <- max(sapply(yy, function(z) max(z[,1], na.rm=TRUE)), na.rm=TRUE)
    xlim <- c(0, maxX)
    maxY <- max(sapply(yy, function(z) max(abs(range(z[,2], na.rm=TRUE)), 
                                          na.rm=TRUE)), na.rm=TRUE)
    minRange <- max(c(diff(range)/20, var)) 
    ylim <- c(min(-maxY, -minRange), max(maxY, minRange))
    if(missing(lwd))
        lwd <- 2
    plot(yy[[1]], xlab="time", type="n", col=col[1], xaxt="n", yaxt="n",
         lwd=lwd, xlim=xlim, ylim=ylim, main=main, ylab="", ...)
    if(varCut>0){
        xl <- par("usr")[1:2]
        xl <- xl + c(1,-1)*(diff(xl)*0.01)
        rect(xl[1], max(var)*varCut, xl[2], -max(var)*varCut,
             col=desat("gray", by=30), border=NA)
    }
    abline(h=0, col="darkgray")
    for(j in 1:length(y))
        lines(yy[[j]], col=col[j], lwd=lwd)
}


## plot stacked values for each flowFrame
stackedPlot <- function(y, p, main=paste("time line for", p),
                        range, col, ylab, med, lwd, varCut, ...){
    par(mar=c(1,5,3,3), mgp=c(2,0.5,0), las=1)
    var <- sapply(y, function(z) median(z$var))
    actualRange <- max(c(diff(range)/10, sapply(y, function(x)
                        diff(range(x$smooth[,2], na.rm=TRUE))), var*2))*1.01
    stacks <- ((1:length(y))-1) * actualRange
    yy <- mapply(function(z, m, s) data.frame(x=z$smooth[,1],
                                              y=z$smooth[,2]-m+s), y,
                 med, stacks, SIMPLIFY=FALSE)
    maxX <- max(sapply(yy, function(z) max(z[,1], na.rm=TRUE)))
    xlim <- c(0, maxX)
    ylim <- range(sapply(yy, function(z) range(z[,2], na.rm=TRUE)))
    if(missing(ylab) | is.null(ylab))
       ylab <- names(y)
    if(missing(lwd))
        lwd <- 2
    plot(yy[[1]], xlab="", ylab="", type="n", xaxt="n", 
         lwd=lwd, xlim=xlim, ylim=ylim, main=main, yaxt="n", ...)
    xl <- par("usr")[1:2]
    xl <- xl + c(1,-1)*(diff(xl)*0.01)
    if(length(ylab)>1)
        axis(2, stacks, ylab, cex.axis=0.8)
    for(j in 1:length(y)){
        if(varCut>0)
            rect(xl[1], mean(yy[[j]]$y)-var[[j]]*varCut, xl[2],
                 mean(yy[[j]]$y)+var[[j]]*varCut,
                 col=desat("gray", by=30), border=NA)
        lines(yy[[j]], col=col[j], lwd=lwd)   
    }
}



## plot values in the "native" dimensions
nativePlot <- function(y, p, main=paste("time line for", p),
                       range, col, med, lwd, varCut, ...){
    par(mar=c(1,2.5,3,2.5), mgp=c(1.5,0.5,0))
    var <- sapply(y, function(z) mean(z$var))
    actualRange <- max(c(diff(range)/10, sapply(y, function(x)
                        diff(range(x$smooth[,2], na.rm=TRUE))),
                         max(var)*2.1))*1.01
    maxX <- max(sapply(y, function(z) max(z$smooth[,1], na.rm=TRUE)),
                na.rm=TRUE)
    xlim <- c(0, maxX)
    m <- mean(sapply(y, function(z) mean(z$smooth[,2])))
    maxY <-  max(c(sapply(y,function(z)
                          max(abs(range(z$smooth[,2],na.rm=TRUE)))),
                   m+max(var)*1.05, m-diff(range)/20))
    minY <-  min(c(sapply(y,function(z)
                        min(abs(range(z$smooth[,2],na.rm=TRUE)))),
                   m-max(var)*1.05, m-diff(range)/20))
    ylim <- c(minY, maxY)
    if(missing(lwd))
        lwd <- 2                  
    plot(y[[1]]$smooth, xlab="", ylab="", type="n", col=col[1], yaxt="n", 
         lwd=lwd, xlim=xlim, ylim=ylim, main=main, xaxt="n", ...)
    if(varCut>0 && length(var)==1){
        xl <- par("usr")[1:2]
        xl <- xl + c(1,-1)*(diff(xl)*0.01)
        
        rect(xl[1], m-max(var)*varCut, xl[2], m+max(var)*varCut,
             col=desat("gray", by=30), border=NA)
    }
    for(j in 1:length(y))
        lines(y[[j]]$smooth, col=col[j], lwd=lwd, ...)
}

