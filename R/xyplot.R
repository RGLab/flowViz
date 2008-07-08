##############################################################################
##                             flowFrame methods                            ##
##############################################################################

## Default xyplot for flowFrames without any formula. This tries to guess the
## 'Time' parameter from the colnames of the expression matrix and creates
## timeline plots for each parameter in a stacked layout. The method casts
## an error if it isn't able to guess the 'Time' parameter. There are three
## options for plotting:
##    discretize (the default): bin y values in time interval and plot median
##                              values
##    smooth: standard smoothScatter plot of time vs parameter (without
##            points added by default)
##    none of the above: standard type argument for xyplots, defaults to 'l'
setMethod("xyplot",
          signature(x="flowFrame", data="missing"),
          function(x, data, xlab=time, ylab="", time="Time",
                   layout, panel=panel.xyplot.flowframe.time,
                   prepanel=prepanel.xyplot.flowframe.time,
                   type="discrete", ...)
      {
          ## guess the time parameter
          expr <- exprs(x)
          if (!(time %in% colnames(expr)))
              stop("Name of the Time (X) variable must be specified as the ",
                   "'time' argument")
          time.x <- expr[, time]
          ## set up fake data frame for dispatch
          fakedf <- data.frame(channel=setdiff(colnames(expr), time),
                               time=1)
          ## set up stacked layout
          if (missing(layout)) layout <- c(1, ncol(expr) - 1)
          xyplot(channel ~ time | channel, data=fakedf, type=type,
                 prepanel=prepanel, panel=panel, layout=layout,
                 time.x=time.x, expr=expr, xlab=xlab, ylab=ylab,
                 default.scales=list(y=list(relation="free", rot=0)),
                 ...)
      })


## Dedicated prepanel function to set up dimensions
prepanel.xyplot.flowframe.time <- 
    function(x, y, time.x, expr, ...)
{
    xx <- time.x
    yy <- expr[, as.character(y)]
    list(xlim=range(xx, finite=TRUE), ylim=range(yy, finite=TRUE),
         dx=diff(xx), dy=diff(yy))
}

## Dedicated panel function to do the timeline plotting with
## different options
panel.xyplot.flowframe.time <- 
    function(x, y, time.x, expr, type="l", nrpoints=0, binSize=100, ...)
{
    xx <- time.x
    yy <- expr[, as.character(y)]
    plotType("gtime", c("Time",  as.character(y)))
    if (type == "smooth"){
        ## smoothScatter plot of time vs. parameter (not sure how useful
        ## that is...)
        panel.smoothScatter(xx, yy, nrpoints=nrpoints, ...)
        return(invisible(NULL))
    }else if(type == "discrete"){
        ## discetize y values into bins acording to time intervals and
        ## compute median values for each interval  
        ord <- order(xx)
        xx <- xx[ord]
        yy <- yy[ord]
        lenx <- length(xx)
        nrBins <- floor(lenx/binSize)
        ## time parameter is already binned or very sparse events
        if(length(unique(xx)) < nrBins){
            yy <- sapply(split(yy, xx), median)
            xx <- unique(xx)
        }else{
            ## bin values in nrBins bins
            if(lenx > binSize){
                cf <- c(rep(1:nrBins, each=binSize),
                        rep(nrBins+1, lenx-nrBins*binSize))
                stopifnot(length(cf) == lenx)
                tmpx <- split(xx,cf)
                tmpy <- split(yy,cf)
                yy <- sapply(tmpy, median, na.rm=TRUE)
                xx <- sapply(tmpx, mean, na.rm=TRUE)
            }else{
                ## very little events
                warning("Low number of events", call.=FALSE)
                tmpy <- split(yy,xx)
                yy <- sapply(tmpy, median, na.rm=TRUE)
                xx <- unique(xx)
            }
        }
        panel.xyplot(xx, yy, type="l", ...)
    } else{
        ## regular xyplot using lines
        panel.xyplot(xx, yy, type=type, ...)
    }
}



## xyplot for flowFrames with a formula.  We'll make this very simple;
## the drawback being that the expression matrix will be copied,
## the upshot being that all the fancy xyplot formula stuff will be valid.
setMethod("xyplot",
          signature(x="formula", data="flowFrame"),
          function(x, data, smooth=TRUE, panel=panel.xyplot.flowframe,
                   gpar=flowViz.par(), xlim, ylim, ...)
      {
          ## deparse the formula structure
          channel.y <- x[[2]]
          channel.x <- x[[3]]
          if (length(channel.x) == 3)
              channel.x <- channel.x[[2]]
          channel.x.name <- expr2char(channel.x)
          channel.y.name <- expr2char(channel.y)
          ## find useful xlim and ylim defaults
          if(missing(xlim)){
              xlim <- unlist(range(data, channel.x.name))
              xd <- diff(xlim)/20
              xlim <- xlim + c(-1,1)*xd
          }
          if(missing(ylim)){
              ylim <- unlist(range(data, channel.y.name))
              yd <- diff(ylim)/20
              ylim <- ylim + c(-1,1)*yd
          }
          xyplot(x, data=as.data.frame(exprs(data)), smooth=smooth, 
                 panel=panel, frame=data, channel.x.name=channel.x.name,
                 channel.y.name=channel.y.name, gpar=gpar, xlim=xlim,
                 ylim=ylim, ...)
      })


## Dedicated panel function to be abble to add filters on the plot
panel.xyplot.flowframe <- function(x,y, frame, filter=NULL, smooth=TRUE,
                                   channel.x.name,  channel.y.name,
                                   pch=".", gpar, margin=TRUE, ...)
{
    if (smooth){
        ## We remove margin events before computing the density and
        ## indicate those by lines
        if(margin){
            r <- range(frame, c(channel.x.name, channel.y.name))
            l <- length(x)
            inc <- apply(r,2,diff)/1e5
            dots <- list(...)
            nb <- if("nbin" %in% names(dots)) rep(dots$nbin,2) else rep(64,2)
            selxL <- x > r[2,channel.x.name]-inc[1]
            selxS <- x < r[1,channel.x.name]+inc[1]
            selyL <- y > r[2,channel.y.name]-inc[2]
            selyS <- y < r[1,channel.y.name]+inc[2]
            allsel <- !(selxL | selxS | selyL | selyS)
            panel.smoothScatter(x[allsel], y[allsel],
                                range.x=list(r[,1], r[,2]), ...)
            addMargin(r[1,channel.x.name], y[selxS], r, l, nb[1])
            addMargin(r[2,channel.x.name], y[selxL], r, l, nb[1], b=TRUE)
            addMargin(x[selyS], r[1,channel.y.name], r, l, nb[2])
            addMargin(x[selyL], r[2,channel.y.name], r, l, nb[2], b=TRUE)
        }else{
            panel.smoothScatter(x, y, ...)
        }
        plotType("gsmooth", c(channel.x.name, channel.y.name))
        if(!is.null(filter)){
            glpolygon(filter, frame,
                      channels=c(channel.x.name, channel.y.name),
                      verbose=FALSE, gpar=gpar, ...)
        }
    }else{
        panel.xyplot(x, y, pch=pch, ...)
        plotType("gpoints", c(channel.x.name, channel.y.name))
        if(!is.null(filter)){
            glpoints(filter, frame,
                     channels=c(channel.x.name, channel.y.name),
                     verbose=FALSE, pch=pch, col="red",
                     gpar=gpar, ...)
        }
    }
}




##############################################################################
##                            flowSet methods                               ##
##############################################################################

## xyplot method for flowSets with formula. 
setMethod("xyplot",
          signature(x="formula", data="flowSet"),
          function(x, data, xlab, ylab, as.table=TRUE,
                   prepanel=prepanel.xyplot.flowset,
                   panel=panel.xyplot.flowset, pch = ".",
                   smooth=TRUE, filter=NULL,
                   gpar=NULL, ...)
      {
          ## ugly hack to suppress warnings about coercion introducing
          ## NAs (needs to be `undone' inside prepanel and panel
          ## functions):
          pd <- pData(data)
          uniq.name <- createUniqueColumnName(pd)
          pd[[uniq.name]] <- factor(sampleNames(data))
          ## deparse the formula structure
          channel.y <- x[[2]]
          channel.x <- x[[3]]
          if (length(channel.x) == 3)
          {
              channel.x <- channel.x[[2]]
              x[[3]][[2]] <- as.name(uniq.name)
              x[[2]] <- NULL
          }
          else
          {
              x[[3]] <- as.name(uniq.name)
              x[[2]] <- NULL
          }
          channel.x.name <- expr2char(channel.x)
          channel.y.name <- expr2char(channel.y)
          channel.x <- as.expression(channel.x)
          channel.y <- as.expression(channel.y)
          ## axis annotation
          if (missing(xlab)) xlab <- channel.x.name
          if (missing(ylab)) ylab <- channel.y.name
          ## use densityplot method with dedicated panel and prepanel
          ## functions to do the actual plotting
          densityplot(x, data=pd, prepanel=prepanel, panel=panel,
                      frames=data@frames, channel.x=channel.x,
                      channel.y=channel.y, channel.x.name=channel.x.name,
                      channel.y.name=channel.y.name, filter=filter,
                      as.table=as.table, xlab=xlab, ylab=ylab,
                      pch=pch, smooth=smooth, gpar=gpar, ...)
      })



## Dedicated prepanel function to set up dimensions
prepanel.xyplot.flowset <- 
    function(x,  frames, channel.x, channel.y, ...)
{
    if (length(nm <- as.character(x)) > 1)
        stop("must have only one flow frame per panel")
    if (length(nm) == 1)
    {
        xx <- evalInFlowFrame(channel.x, frames[[nm]])
        yy <- evalInFlowFrame(channel.y, frames[[nm]])
        list(xlim=range(xx, finite=TRUE), ylim=range(yy, finite=TRUE),
             dx=diff(xx), dy=diff(yy))
    }
    else list()
}



## Dedicated panel function to do the plotting and add gate boundaries
panel.xyplot.flowset <-
    function(x, frames, channel.x, channel.y, channel.x.name, channel.y.name, 
             filter=NULL, pch=".", smooth=TRUE, gpar, margin=TRUE, ...)
{
    x <- as.character(x)
    if (length(x) > 1) stop("must have only one flow frame per panel")
    if (length(x) < 1) return()
    nm <- x
    if(is(filter, "filter")){
        filter <- list(filter)
        names(filter) <- nm
    }
    if(!is.null(filter) && (!(nm %in% names(filter) ||
                              !is(filter[[nm]] ,"filter"))))
        stop("'filter' must inherit from class 'filter' or a be list of ",
             "such objects.", call.=FALSE) 
    xx <- flowViz:::evalInFlowFrame(channel.x, frames[[nm]])
    yy <- flowViz:::evalInFlowFrame(channel.y, frames[[nm]])

    if(smooth) {
        ## We remove margin events before computing the density and
        ## indicate those by lines
        if(margin){
            r <- range(frames[[nm]], c(channel.x.name, channel.y.name))
            l <- length(xx)
            inc <- apply(r,2,diff)/1e5
            dots <- list(...)
            nb <- if("nbin" %in% names(dots)) rep(dots$nbin,2) else rep(64,2)
            selxL <- xx > r[2,channel.x.name]-inc[1]
            selxS <- xx < r[1,channel.x.name]+inc[1]
            selyL <- yy > r[2,channel.y.name]-inc[2]
            selyS <- yy < r[1,channel.y.name]+inc[2]
            allsel <- !(selxL | selxS | selyL | selyS)
            panel.smoothScatter(xx[allsel], yy[allsel],
                                range.x=list(r[,1], r[,2]), ...)
            addMargin(r[1,channel.x.name], yy[selxS], r, l, nb[1])
            addMargin(r[2,channel.x.name], yy[selxL], r, l, nb[1], b=TRUE)
            addMargin(xx[selyS], r[1,channel.y.name], r, l, nb[2])
            addMargin(xx[selyL], r[2,channel.y.name], r, l, nb[2], b=TRUE)
        }else
        panel.smoothScatter(xx, yy, ...)
        plotType("gsmooth", c(channel.x.name, channel.y.name))
        if(!is.null(filter)){
            if(is.null(gpar))
                gpar <- flowViz.par()
            glpolygon(filter[[nm]], frames[[nm]],
                      channels=c(channel.x.name, channel.y.name),
                      verbose=FALSE, gpar=gpar, ...)
        }
    }
    else{
        panel.xyplot(xx, yy, pch=pch, ...)
        plotType("gpoints", c(channel.x.name, channel.y.name))
        col = brewer.pal(9, "Set1")
         if(!is.null(filter)){
             glpoints(filter[[nm]], frames[[nm]],
                      channels=c(channel.x.name, channel.y.name),
                      verbose=FALSE, gpar=gpar, pch=pch, ...)
         }
    }
}



addMargin <- function(x, y, r, total, nb, len=200, b=FALSE)
{
    if(length(x) & length(y)){
        dx <-  diff(r[,1])
        dy <- diff(r[,2])
        lenx <- dx/len
        leny <- dy/len
        addX <- if(b) dx/nb[1]*1.1 else 0
        addY <- if(b) dy/nb[2]*1.1 else 0
        n <- 100
        colvec <- c(NA, colorRampPalette(c("lightgray", "black"))(99)) 
        if(length(x)==1){
            hh <- hist(y, n=n, plot=FALSE)
            xoff <- dx/(nb*2)-addX
            col <- colvec[pmin(100, as.integer(hh$counts/total*5000)+1)]
            panel.segments(rep(x-xoff, n), hh$mids, rep(x-xoff-lenx, n),
                           hh$mids, col=col, lwd=2)
        }else{
            hh <- hist(x, n=n, plot=FALSE)
            yoff <- dy/(nb*2)-addY
            col <- colvec[pmin(100, as.integer(hh$counts/total*5000)+1)]
            panel.segments(hh$mids, rep(y-yoff, n), hh$mids,
                           rep(y-yoff-leny, n), col=col, lwd=2)
        }
    }
}








