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
                   layout, panel=panel.xyplot.flowframe,
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
                 prepanel=prepanel.xyplot.flowframe.time,
                 panel=panel.xyplot.flowframe.time, layout=layout,
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

## Dedicated panel function to do the plotting with different options
panel.xyplot.flowframe.time <- 
    function(x, y, time.x, expr, type="l", nrpoints=0, binSize=100, ...)
{
    xx <- time.x
    yy <- expr[, as.character(y)]
    
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
          function(x, data, smooth=TRUE, 
                   panel=if (smooth) panel.smoothScatter else panel.xyplot,
                   ...)
      {
          xyplot(x, data=as.data.frame(exprs(data)), smooth=smooth, 
                 panel=panel, ...)
      })







##############################################################################
##                            flowSet methods                               ##
##############################################################################

## xyplot method for flowSets with formula. 
setMethod("xyplot",
          signature(x="formula", data="flowSet"),
          function(x, data, xlab, ylab, as.table=TRUE,
                   prepanel=prepanel.xyplot.flowset,
                   panel=panel.xyplot.flowset, pch = ".",
                   smooth=TRUE, filter=NULL, filterResults=NULL,
                   displayFilter=TRUE, ...)
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
                      filterResults=filterResults, displayFilter=displayFilter,
                      as.table=as.table, xlab=xlab, ylab=ylab,
                      pch=pch, smooth=smooth, ...)
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
             filter=NULL, filterResults=NULL, pch=".", smooth=TRUE, ...)
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
        panel.smoothScatter(xx, yy, ...)
        if(!is.null(filter)){
            glpolygon(filter[[nm]], frames[[nm]],
                      channels=c(channel.x.name, channel.y.name),
                      verbose=FALSE, ...)
        }
    }
    else{
        panel.xyplot(xx, yy, pch=pch, ...)
        col = brewer.pal(9, "Set1")
         if(!is.null(filter)){
             glpoints(filter[[nm]], frames[[nm]],
                      channels=c(channel.x.name, channel.y.name),
                      verbose=FALSE, pch=pch, col=col, ...)
         }
    }
}












