## Plot values of a flowSet against the time domain.
## TO DO: Eventualy, we might want individual methods for flowSets and for
## flowFrames, although it might be worthwhile to reimplement using
## lattice graphics first...


## run over a cytoFrame, bin the values according to the time domain
## and compute medians for each bin. The resultof this function will
## be the input to the plotting functions
prepareSet <- function(x, parm){
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
    bs <- 500
    if(lenx>bs){
        bc <- floor(lenx/bs)
        cf <- c(rep(1:bc, each=bs), rep(bc+1, lenx-bc*bs))
        stopifnot(length(cf) == lenx)
        tmpx <- split(xx,cf)
        tmpy <- split(yy,cf)
        yy <- sapply(tmpy, median, na.rm=TRUE)
        xx <- sapply(tmpx, mean, na.rm=TRUE)
    }else{
        warning("Low number of events", call.=FALSE)
        tmpy <- split(yy,xx)
        yy <- sapply(tmpy, median, na.rm=TRUE)
        xx <- unique(xx)
    }
    
    return(cbind(xx,yy))
}


## wrapper function to produce the plots
timeLinePlot <- function(x, channel, type=c("stacked", "scaled", "native"), col,
                         ylab=names(x), ...){
    if(!channel %in% colnames(x[[1]]))
        stop(channel, " is not a valid channel in this flowSet.")
    if(missing(col)){
        require(RColorBrewer)
        colp <- brewer.pal(12, "Paired")
        col <- colorRampPalette(colp)(length(x))
        set.seed(1000)
        col <- sample(col)				 
    }else{
        if(length(col)!=1 || length(col)!=length(x))
            stop("'col' must be color vector of length 1 or same length ",
                 "as the flowSet")
    }
    opar <- par(c("mar", "mgp", "mfcol", "mfrow", "las"))
    on.exit(par(opar))
    timeData <-  fsApply(x, prepareSet, channel, use.exprs=TRUE, simplify=FALSE)
    type <- match.arg(type)
    mr <- range(x[[1]])[,channel]
    mr[1] <- max(mr[1], 0)
    med <- sapply(timeData, function(z) median(z[,2], na.rm=TRUE))
    if(length(med)==1){
        nativePlot(timeData, p=channel, range=mr, col=col, med=med, ...)
        return(sum(abs(timeData[[1]][,2]-med), na.rm=TRUE)/nrow(timeData[[1]]))
    }
    layout(matrix(1:2), heights=c(0.8, 0.2))
    switch(type,
           "scaled"=scaledPlot(timeData, p=channel, range=mr, col=col,
           med=med, ...),
           "stacked"=stackedPlot(timeData, p=channel, range=mr, col=col,
           ylab=ylab, med=med, ...),
           "native"=nativePlot(timeData, p=channel, range=mr, col=col,
           med=med, ...),
           stop("Unknown type"))
    ad <- mapply(function(z,m) sum(abs(z[,2]-m), na.rm=TRUE)/nrow(z),
                 timeData, med)
    par(mar=c(5,3,0,3), las=2)
    on.exit(par(opar))
    barplot(ad, axes=FALSE, col=col, cex.names=0.8,
            ylim=c(0, mr[2]/25), border=col)
    box(col="darkgray")
    return(ad)
}


## align values around 0 and plot
scaledPlot <- function(y, p, main=paste("time line for", p),
                       range, col, med, ...){
    par(mar=c(3.5,2.5,3,2.5), mgp=c(1.5,0.5,0))
    y <- mapply(function(z, m) data.frame(x=z[,1], y=z[,2]-m), y, med,
                SIMPLIFY=FALSE)
    maxX <- max(sapply(y, function(z) max(z[,1], na.rm=TRUE)), na.rm=TRUE)
    xlim <- c(0, maxX)
    maxY <- max(sapply(y, function(z) max(abs(range(z[,2], na.rm=TRUE)), 
                                          na.rm=TRUE)), na.rm=TRUE)
    minRange <- diff(range)/4
    ylim <- c(min(-maxY, -minRange), max(maxY, minRange))
    plot(y[[1]], xlab="time", type="l", col=col[1], 
         lwd=1, xlim=xlim, ylim=ylim, main=main, ylab="", ...)
    abline(h=0, col="darkgray")
    for(j in 2:length(y))
        lines(y[[j]], col=col[j], lwd=1)
}


## plot stacked values for each flowFrame
stackedPlot <- function(y, p, main=paste("time line for", p),
                        range, col, ylab, med, ...){
    par(mar=c(4,5,3,3), mgp=c(2,0.5,0), las=1)
    stacks <- length(y):1 * (diff(range)/10)
    y <- mapply(function(z, m, s) data.frame(x=z[,1], y=z[,2]-m+s), y, med, 
			stacks, SIMPLIFY=FALSE)
    maxX <- max(sapply(y, function(z) max(z[,1], na.rm=TRUE)), na.rm=TRUE)
    xlim <- c(0, maxX)
    ylim <- range(sapply(y, function(z) range(z[,2], na.rm=TRUE)), na.rm=TRUE)
    if(missing(ylab) | is.null(ylab))
       ylab <- names(y)
    plot(y[[1]], xlab="time", ylab="", type="l", col=col[1], 
         lwd=1, xlim=xlim, ylim=ylim, main=main, yaxt="n", ...)
    if(length(ylab)>1)
        axis(2, stacks, ylab, cex.axis=0.8)
    for(j in 2:length(y))
        lines(y[[j]], col=col[j], lwd=1)
}


## plot values in the "native" dimensions
nativePlot <- function(y, p, main=paste("time line for", p),
                       range, col, med, ...){
    par(mar=c(3.5,2.5,3,2.5), mgp=c(1.5,0.5,0))
    maxX <- max(sapply(y, function(z) max(z[,1], na.rm=TRUE)), na.rm=TRUE)
    xlim <- c(0, maxX)
    maxY <-  max(sapply(y, function(z) max(abs(range(z[,2], na.rm=TRUE)), 
                                           na.rm=TRUE)), na.rm=TRUE)
    minRange <- diff(range)/4
    ylim <- c(0, max(minRange, maxY))
    plot(y[[1]], xlab="time", ylab="", type="l", col=col[1], 
         lwd=1, xlim=xlim, ylim=ylim, main=main, ...)
    if(length(y)>1)
        for(j in 2:length(y))
            lines(y[[j]], col=col[j], lwd=1)
}

