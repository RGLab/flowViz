## Plot values of a flowSet against the time domain.
## TO DO: Eventualy, we might want individual methods for flowSets and for
## flowFrames, although it might be worthwhile to reimplement using
## lattice graphics first...



## Call flowCore:::prepareSet for each flowFrame in a flowSet and produce the
## plots, according to "type" (see below for details). Also compute a
## quality value modified by the setting of varCut, which is basically
## the sum of average positive distances of means for each bin from the global
## mean confidence interval, scaled by the variance cutoff:
##       sum(z[z>0])/length(z)/varCut
## For type "frequency", this is the sum of events per time tick for
## bins where the total number of events is > 2 times the expected average
## number
## A spearman rank correlation coefficient is return as an attribute.
timelineplot <- function(x, channel, type=c("stacked", "scaled", "native",
                                     "frequency"),
                         col, ylab=names(x), binSize, varCut=1, ...)
{
    ## Sanity checking up front
    if(!length(channel)==1)
        stop("'channel' must be character scalar")
    if(!channel %in% colnames(x[[1]]))
        stop(channel, " is not a valid channel in this flowSet.")
    if(tolower(channel) == "time")
        stop("Argument 'channel' can not be the time channel")

    ## Making sure the sample names fit on the plot as axis annotation
    sampleNames(x) <- truncNames(sampleNames(x))

    ## Lets fix ourselves some nice colors
    if(missing(col) | is.null(col)){
#        require(RColorBrewer)
        colp <- brewer.pal(8, "Dark2")
        col <- colorRampPalette(colp)(length(x))
        set.seed(1000)
        col <- sample(col) 
    }else{
        if(length(col)!=1 || length(col)!=length(x))
            stop("'col' must be color vector of length 1 or same length ",
                 "as the flowSet")
    }

    ## Bin the data and compute local variances and locations
    time <- flowCore:::findTimeChannel(xx= exprs(x[[1]]))
    timeData <- fsApply(x, flowCore:::prepareSet, parm=channel, time=time,
                        binSize=binSize, use.exprs=FALSE, simplify=FALSE)
    opar <- par(c("mar", "mgp", "mfcol", "mfrow", "las"))
    on.exit(par(opar))
    type <- match.arg(type)
    mr <- range(x[[1]])[,channel]
    mr[1] <- max(mr[1], 0)
    
    ## Standardize to compute meaningful scale-free QA scores
    med <- sapply(timeData, function(z) median(z$smooth[,2], na.rm=TRUE))
    lm <- length(med)
    gvars <- sapply(timeData, function(x) mean(x$variance))
    stand <-  mapply(function(z,m,v) abs(z$smooth[,2]-m)/(v*varCut), timeData,
                 med, gvars, SIMPLIFY=FALSE)
    tvals <- lapply(timeData, function(x) x$smooth[,1])
    
    ## Create the plot, either one of the 4 possible types. For flowFrames
    ## we always use the native scaling.
    if(lm==1 && type != "frequency"){
        nativePlot(timeData, p=channel, range=mr, col="darkblue",
                   varCut=varCut, ...)
        return((sum(unlist(stand)[unlist(stand)>0])/length(stand))/varCut)
    }else if(lm>1)
        layout(matrix(1:2), heights=c(0.8, 0.2))
    switch(type,
           "scaled"=scaledPlot(timeData, p=channel, range=mr, col=col,
           med=med, varCut=varCut, ...),
           "stacked"=stackedPlot(timeData, p=channel, range=mr, col=col,
           ylab=ylab, med=med, varCut=varCut, ...),
           "native"=nativePlot(timeData, p=channel, range=mr, col=col,
           varCut=varCut, ...),
           "frequency"={
               freqPlot(timeData, p=channel, col=col, varCut=varCut,
                        ylab=ylab, ...)
               stand <- lapply(timeData, function(x)
                               x$frequencies[,2] /
                               (mean(x$frequencies[,2])*varCut)-1)
               tvals <- lapply(timeData, function(x) x$frequencies[,1])
           },
           stop("Unknown type"))
    
    ## the QA score and correlation coefficient
    qaScore <- computeQAScore(stand)
    corr <- mapply(cor, tvals, stand, method="spearman", use="pairwise")

    ## a legend indicating the problematic samples
    if(lm>1){
        par(mar=c(5,3,0,3), las=2)
        on.exit(par(opar))
        top <- 2
        barplot(qaScore, axes=FALSE, col=col, cex.names=0.7,
                ylim=c(0, min(c(top, max(qaScore, na.rm=TRUE)))), border=col,
                space=0.2)
        wh <- which(qaScore >= top)
        points((wh+wh*0.2)-0.5, rep(top-(top/12), length(wh)), pch=17, col="white",
               cex=0.7)
    }

    ## The return value with attributes attached.
    attr(qaScore, "binSize") <- binSize
    attr(qaScore, "correlation") <- corr
    attr(qaScore, "raw") <- stand
    return(invisible(qaScore))
}


## Compute the QA score from the standardized values
computeQAScore <- function(stand, varCutoff=1)
    sapply(stand, function(z) sum(z[abs(z) > varCutoff])/length(z))*100

## Truncate the names to make them fit on the plot
truncNames <- function(names){
    nc <- nchar(names)
    names[nc>11] <- paste(substr(names[nc>11], 1, 8), "...", sep="")
    if(any(duplicated(names))){
        ns <- split(names, names)
        is <- split(seq_along(names), names)
        ns<- lapply(ns, function(x)
                    if(length(x)>1) paste(x, seq_along(x), sep="_") else x)
        names <- unlist(ns)[unlist(is)]
    }
    return(names)
}



## Align all values around 0 and plot
scaledPlot <- function(y, p, main=paste("time line for", p),
                       range, col, med, lwd, varCut, ...){
    par(mar=c(1,2.5,3,2.5), mgp=c(1.5,0.5,0))
    yy <- mapply(function(z, m) data.frame(x=z$smooth[,1], y=z$smooth[,2]-m),
                 y, med, SIMPLIFY=FALSE)
    var <- sapply(y, function(z) mean(z$var, na.rm=TRUE))
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
             col=desat("lightgray", by=30), border=NA)
    }
    abline(h=0, col="darkgray")
    for(j in 1:length(y))
        lines(yy[[j]], col=col[j], lwd=lwd)
}


## Plot values for each flowFrame and stack individual results
stackedPlot <- function(y, p, main=paste("time line for", p),
                        range, col, ylab, med, lwd, varCut, ...){
    par(mar=c(1,5,3,3), mgp=c(2,0.5,0), las=1)
    var <- sapply(y, function(z) mean(z$var, na.rm=TRUE))
    actualRange <- max(c(diff(range)/10, sapply(y, function(x)
                        diff(range(x$smooth[,2], na.rm=TRUE))), var*2))*1.01
    stacks <- ((length(y):1)-1) * actualRange
    yy <- mapply(function(z, m, s) data.frame(x=z$smooth[,1],
                                              y=z$smooth[,2]-m+s), y,
                 med, stacks, SIMPLIFY=FALSE)
    maxX <- max(sapply(yy, function(z) max(z[,1], na.rm=TRUE)))
    xlim <- c(0, maxX)
    ylim <- range(c(sapply(yy, function(z) range(z[,2], na.rm=TRUE))),
                  median(yy[[1]]$y)+var[[1]], median(yy[[length(yy)]]$y) -
                                        var[length(yy)])
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
    if(varCut>0)
    for(j in 1:length(y))
        rect(xl[1], mean(yy[[j]]$y)-var[[j]]*varCut, xl[2],
             mean(yy[[j]]$y)+var[[j]]*varCut,
             col=desat("lightgray", by=30), border=NA)
    for(j in 1:length(y))
        lines(yy[[j]], col=col[j], lwd=lwd)   
}



## Plot values in the "native" dimensions
nativePlot <- function(y, p, main=paste("time line for", p),
                       range, col, lwd, varCut, ...){
    par(mar=c(1,2.5,3,2.5), mgp=c(1.5,0.5,0))
    var <- sapply(y, function(z) mean(z$var, na.rm=TRUE))
    actualRange <- max(c(diff(range)/10, sapply(y, function(x)
                        diff(range(x$smooth[,2], na.rm=TRUE))),
                         max(var)*2.1))*1.01
    maxX <- max(sapply(y, function(z) max(z$smooth[,1], na.rm=TRUE)),
                na.rm=TRUE)
    xlim <- c(0, maxX)
    m <- mean(sapply(y, function(z)
                 {
                     sel <- z$smooth[,2]>range[1] & z$smooth[,2]<range[2]
                     if(length(sel))
                         mean(z$smooth[sel,2])
                     else
                         z$smooth[1,2]
                 }))
    maxY <-  max(c(sapply(y,function(z)
                          max(z$smooth[,2],na.rm=TRUE)),
                   m+max(var)*1.05, m-diff(range)/20))
    minY <-  min(c(sapply(y,function(z)
                        min(z$smooth[,2],na.rm=TRUE)),
                   m-max(var)*1.05, m-diff(range)/20))
    ylim <- c(minY, maxY)
    if(missing(lwd))
        lwd <- 2                  
    plot(y[[1]]$smooth, xlab="", ylab="", type="n", col=col[1], yaxt="n", 
         lwd=lwd, xlim=xlim, ylim=ylim, main=main, xaxt="n", ...)
    if(varCut>0 && length(var)==1){
        xl <- par("usr")[1:2]
        xl <- xl + c(1,-1)*(diff(xl)*0.01)
        
        rect(xl[1], m-var*varCut, xl[2], m+var*varCut,
             col=desat("lightgray", by=30), border=NA)
    }
    for(j in 1:length(y))
        lines(y[[j]]$smooth, col=col[j], lwd=lwd, ...)
}



## Plot frequency values for each flowFrame
freqPlot <- function(y, p, main="time line frequencies",
                     col, ylab, lwd, varCut, ...){
    par(mar=c(1,5,3,3), mgp=c(2,0.5,0), las=1)
    var <- 1
    stX <- lapply(y, function(x) x$frequencies[,1])
    stY <- lapply(y, function(x)
                  x$frequencies[,2] / mean(x$frequencies[,2])-1) 
    actualRange <- max(diff(range(stY)), var*varCut*2)*1.01
    stacks <- ((length(y):1)-1) * actualRange
    stYY <- mapply(function(x,s) x+s, stY, stacks, SIMPLIFY=FALSE)
    if(!is.list(stYY))
        stYY <- list(stYY)
    xlim <- c(0, max(unlist(stX)))
    ylim <- range(c(stYY), stacks+var*varCut, stacks-var*varCut)
    if(missing(ylab) | is.null(ylab))
       ylab <- names(y)
    if(missing(lwd))
        lwd <- 2
    plot(stX[[1]], stYY[[1]], xlab="", ylab="", type="n", xaxt="n", 
         lwd=lwd, xlim=xlim, ylim=ylim, main=main, yaxt="n", ...)
    xl <- par("usr")[1:2]
    xl <- xl + c(1,-1)*(diff(xl)*0.01)
    if(length(ylab)>1)
        axis(2, stacks, ylab, cex.axis=0.8)
    if(varCut>0)
    for(j in 1:length(y))
        rect(xl[1], stacks[j]-var*varCut, xl[2],
             stacks[j]+var*varCut,
             col=desat("lightgray", by=30), border=NA)
    for(j in 1:length(y))
        lines(stX[[j]], stYY[[j]], col=col[j], lwd=lwd)
}

## A method for flowSets
setMethod("timeLinePlot",
          signature(x="flowSet", channel="character"),
          function(x, channel, type=c("stacked", "scaled", "native",
                               "frequency"),
                   col=NULL, ylab=sampleNames(x), binSize, varCut=1, ...)
      {
          ## a reasonable default for the bin size
          if(missing(binSize))
              binSize <- min(max(1, floor(median(fsApply(x, nrow)/100))), 500)
          type <- match.arg(type)
          timelineplot(x, channel, binSize=binSize, col=col,
                       varCut=varCut, type=type, ylab = ylab, ...)
      })


## A method for flowFrames. We coerce into a flowSet and use that method
setMethod("timeLinePlot",
          signature(x="flowFrame", channel="character"),
          function(x, channel, ...)
      {
          timeLinePlot(as(x, "flowSet"), channel=channel, ...)
      })

## An error if no channel is specified
setMethod("timeLinePlot",
          signature(x="ANY", channel="missing"),
          function(x, channel, ...)
          stop("Argument 'channel' is missing", call.=FALSE))

