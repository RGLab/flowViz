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
          signature=signature(x="flowFrame",
                              data="missing"),
          definition=function(x,
                              data,
                              time,
                              xlab,
                              ylab="",
                              layout,
                              prepanel=prepanel.xyplot.flowframe.time,
                              panel=panel.xyplot.flowframe.time,
                              type="discrete",
                              ...)
      {
          ## guess the time parameter
          expr <- exprs(x)
          if(missing(time))
              time <- flowCore:::findTimeChannel(expr)
          if(!(time %in% colnames(expr)))
              stop("Invalid name of variable (", time, ") recording the ",
                   "\ntime domain specified as 'time' argument.", call.=FALSE)
          if(missing(xlab))
              xlab <- time
          ## set up fake data frame for dispatch
          fakedf <- data.frame(channel=setdiff(colnames(expr), time),
                               time=1)
          ## set up stacked layout
          if (missing(layout)) layout <- c(1, ncol(expr) - 1)
          xyplot(channel ~ time | channel, data=fakedf, type=type,
                 prepanel=prepanel, panel=panel, layout=layout, frame=x,
                 time=time, xlab=xlab, ylab=ylab,
                 default.scales=list(y=list(relation="free", rot=0)), ...)
      })


## Prepanel function to set up dimensions. We want to use the instrument measurement
## range instead of the absolute range of the data for the y dimension. The x-dimension
## is constant since we have the same time recordings for each channel. We also record
## the data ranges in the internal state environment for further use.
prepanel.xyplot.flowframe.time <- 
    function(x, y, frame, time, ...)
{
    yc <- as.character(y)
    expr <- exprs(frame)
    xx <- expr[, time]
    yy <- expr[, yc]
    xlim <- range(xx, finite=TRUE)
    ylim <- unlist(range(frame, yc))
    plotLims(xlim, ylim)
    return(list(xlim=xlim, ylim=ylim))
}

## Panel function to do the timeline plotting with several different options
panel.xyplot.flowframe.time <- 
    function(x, y, frame, time, type="discrete", nrpoints=0, binSize=100, ...)
{
    y <- as.character(y)
    expr <- exprs(frame)
    xx <- expr[, time]
    yy <- expr[, y]
    ## We record in the state environment which type of plot we produce
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
          signature=signature(x="formula",
                              data="flowFrame"),
          definition=function(x,
                              data,
                              smooth=TRUE,
                              prepanel=prepanel.xyplot.flowframe,
                              panel=panel.xyplot.flowframe,
                              ...)
      {
          ## par.settings will not be passed on to the panel functions, so
          ## we have to fetch it from ... and stick the gate relevant stuff
          ## back it in there manually
          gp <- list(...)[["par.settings"]]
          
          ## deparse the formula structure
          ## FIXME: shouldn't all.vars be helpful here?!?
          channel.y <- x[[2]]
          channel.x <- x[[3]]
          if (length(channel.x) == 3)
              channel.x <- channel.x[[2]]
          channel.x.name <- expr2char(channel.x)
          channel.y.name <- expr2char(channel.y)
          ## call data.frame xyplot method with our panel function
          xyplot(x, data=as.data.frame(exprs(data)), smooth=smooth, 
                 prepanel=prepanel, panel=panel, frame=data,
                 channel.x.name=channel.x.name,
                 channel.y.name=channel.y.name,
                 gp=gp, ...)
      })



## Prepanel function to set up dimensions. We want to use the instrument measurement
## range instead of the absolute range of the data. We also record the data ranges
## in the internal state environment for further use.
prepanel.xyplot.flowframe <- 
    function(frame, channel.x.name, channel.y.name, ...)
{
    ranges <- range(frame)
    xlim <- if(!length(grep("`", channel.x.name, fixed=TRUE))){
        tmp <- ranges[, channel.x.name]
        xd <- diff(tmp)/15
        tmp + c(-1,1)*xd
    }else NULL
    ylim <- if(!length(grep("`", channel.y.name, fixed=TRUE))){
        tmp <- ranges[, channel.y.name]
        yd <- diff(tmp)/15
        tmp + c(-1,1)*yd
    }else NULL
    plotLims(xlim, ylim)
    return(list(xlim=xlim, ylim=ylim))
}

## Panel function that allows us to add filters on the plot. The actual plotting
## is done by either the panel.smoothScatter or the default lattice panel.xyplot
## function
panel.xyplot.flowframe <- function(x,
                                   y,
                                   frame,
                                   filter=NULL,
                                   smooth=TRUE,
                                   margin=TRUE,
                                   outline=FALSE,
                                   channel.x.name,
                                   channel.y.name,
                                   pch=gpar$flow.symbol$pch,
                                   alpha=gpar$flow.symbol$alpha,
                                   cex=gpar$flow.symbol$cex,
                                   col=gpar$flow.symbol$col,
                                   gp,
                                   ...)
{
    ## graphical parameter defaults
    gpar <- flowViz.par.get()
    if(!is.null(gp))
        gpar <- lattice:::updateList(gpar, gp)
    if(is.null(gpar$gate$cex))
        gpar$gate$cex <- cex
    if(is.null(gpar$gate$pch))
        gpar$gate$pch <- pch
    ## Whenever we have a function call in the formula we might no longer be the
    ## original scale in which the gate was defined, and we have no clue how to
    # plot it
    validName <- !(length(grep("\\(", channel.x.name)) ||
                   length(grep("\\(", channel.y.name)))
    if (smooth){
        ## We remove margin events before passing on the data to panel.smoothScatter
        ## and after plotting indicate those events by grayscale lines on the plot
        ## margins
        if(margin){
            r <- range(frame, c(channel.x.name, channel.y.name))
            l <- length(x)
            inc <- apply(r, 2, diff)/1e5
            dots <- list(...)
            nb <- if("nbin" %in% names(dots)) rep(dots$nbin, 2) else rep(64, 2)
            selxL <- x > r[2,channel.x.name]-inc[1]
            selxS <- x < r[1,channel.x.name]+inc[1]
            selyL <- y > r[2,channel.y.name]-inc[2]
            selyS <- y < r[1,channel.y.name]+inc[2]
            allsel <- !(selxL | selxS | selyL | selyS)
            panel.smoothScatter(x[allsel], y[allsel], range.x=list(r[,1], r[,2]),
                                ...)
            addMargin(r[1,channel.x.name], y[selxS], r, l, nb)
            addMargin(r[2,channel.x.name], y[selxL], r, l, nb, b=TRUE)
            addMargin(x[selyS], r[1,channel.y.name], r, l, nb)
            addMargin(x[selyL], r[2,channel.y.name], r, l, nb, b=TRUE)
        }else{
            panel.smoothScatter(x, y, ...)
        }
        plotType("gsmooth", c(channel.x.name, channel.y.name))
        if(!is.null(filter) & validName){
            glpolygon(filter, frame,
                      channels=c(channel.x.name, channel.y.name),
                      verbose=FALSE, gpar=gpar, strict=FALSE, ...)
        }
    }else{
        if(!is.null(filter) && validName){
            if(!is(filter, "filterResult"))
                filter <- filter(frame, filter)
            rest <- Subset(frame, !filter)
            x <- exprs(rest[,channel.x.name])
            y <- exprs(rest[,channel.y.name])
            panel.xyplot(x, y, col=col, cex=cex, pch=pch, alpha=alpha, ...)
            
            glpoints(filter, frame,
                     channels=c(channel.x.name, channel.y.name),
                     verbose=FALSE, gpar=gpar, strict=FALSE, ...)
            if(outline)
                glpolygon(filter, frame,
                          channels=c(channel.x.name, channel.y.name),
                          verbose=FALSE, gpar=gpar, names=FALSE,
                          strict=FALSE)
        }else{
            panel.xyplot(x, y, col=col, cex=cex, pch=pch, alpha=alpha, ...)
        }
        plotType("gpoints", c(channel.x.name, channel.y.name))
    }
}




##############################################################################
##                            flowSet methods                               ##
##############################################################################
## xyplot method for flowSets with formula. 
setMethod("xyplot",
          signature=signature(x="formula",
                              data="flowSet"),
          definition=function(x,
                              data,
                              smooth=TRUE,
                              filter=NULL,
                              as.table=TRUE,
                              prepanel=prepanel.xyplot.flowset,
                              panel=panel.xyplot.flowset,
                              xlab=channel.x.name,
                              ylab=channel.y.name,
                              par.settings=NULL,
                              ...)
      {
          ## no conditioning variable, we chose 'name' as default
          if (length(x[[3]]) == 1){
              tmp <- x[[3]]
              x[[3]] <- (~dummy | name)[[2]]
              x[[3]][[2]] <- tmp
          }
          ## par.settings will not be passed on to the panel functions, so
          ## we have to fetch it from ... and stick the gate relevant stuff
          ## back it in there manually
          gp <- par.settings
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
          ## use densityplot method with dedicated panel and prepanel
          ## functions to do the actual plotting
          densityplot(x, data=pd, prepanel=prepanel, panel=panel,
                      frames=data@frames, channel.x=channel.x,
                      channel.y=channel.y, channel.x.name=channel.x.name,
                      channel.y.name=channel.y.name, xlab=xlab, ylab=ylab,
                      smooth=smooth, gp=gp, as.table=as.table, filter=filter,
                      par.settings=par.settings, ...)
      })



## Prepanel function to set up dimensions. We want to use the instrument measurement
## range instead of the absolute range of the data. We also record the data ranges
## in the internal state environment for further use.
prepanel.xyplot.flowset <- 
    function(x, frames, channel.x.name, channel.y.name, ...)
{
    if (length(nm <- as.character(x)) > 1)
        stop("must have only one flow frame per panel")
    if (length(nm) == 1){
        ranges <- range(frames[[nm]])
        xlim <- if(!length(grep("`", channel.x.name, fixed=TRUE))){
            tmp <- ranges[, channel.x.name]
            xd <- diff(tmp)/15
            tmp + c(-1,1)*xd
        }else NULL
        ylim <- if(!length(grep("`", channel.y.name, fixed=TRUE))){
            tmp <- ranges[, channel.y.name]
            yd <- diff(tmp)/15
            tmp + c(-1,1)*yd
        }else NULL
        plotLims(xlim, ylim)
        return(list(xlim=xlim, ylim=ylim))
    }
}



## FIXMES:
##   - How can we cleanly deparse the formula to get to the channel
##     names?
##   - The Formula interface allows for arbitrary functions to be
##     called on the data before plotting. If that happens, we have no
##     clue what to do with the gates, since they are still defined in
##     the old coordinate system.
##   - The same is true for the range parameters stored in the flowFrame.
##     These only make sense in the original coordinates. Do we want to
##     transform them, and if so, how can we do that? Can we substitute
##     the flow data by a vector containing the ranges in the formula? 
## Panel function that allows us to add filters on the plot. The actual plotting
## is done by panel.xyplot.flowframe
panel.xyplot.flowset <- function(x,
                                 frames,
                                 filter=NULL,
                                 channel.x,
                                 channel.y,
                                 ...)
{
    nm <- as.character(x)
    if (length(nm) < 1) return()
    ## 'filter' either has to be a single filter, or a list of filters matching
    ## the flowSet, or a filterResultList.
    if(!is.null(filter)){
        if(!is.list(filter)){
            if(is(filter, "filter")){
                filter <- lapply(seq_along(nm), function(x) filter)
                names(filter) <- nm
            }
        }else if(!is(filter, "filterResultList"))
            filter <- as(filter, "filterResultList")
        if(!(nm %in% names(filter) || !is(filter[[nm]] ,"filter"))){
            warning("'filter' must either be a filterResultList, a single\n",
                    "filter object or a named list of filter objects.",
                    call.=FALSE)
            filter <- NULL
        }
    }
    x <- flowViz:::evalInFlowFrame(channel.x, frames[[nm]])
    y <- flowViz:::evalInFlowFrame(channel.y, frames[[nm]])
    panel.xyplot.flowframe(x, y, frame=frames[[nm]], filter=filter[[nm]], ...)
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
        colvec <- colorRampPalette(c("lightgray", "black"))(100) 
        if(length(x)==1){
            hh <- hist(y, n=n, plot=FALSE)
            sel <- hh$counts > 0
            xoff <- dx/(nb[1]*2)-addX
            col <- colvec[pmin(100, as.integer(hh$counts[sel]/total*5000)+1)]
            panel.segments(rep(x-xoff, n), hh$mids[sel], rep(x-xoff-lenx, n),
                           hh$mids[sel], col=col, lwd=3, lineend=2)
        }else{
            hh <- hist(x, n=n, plot=FALSE)
            sel <- hh$counts > 0
            yoff <- dy/(nb[2]*2)-addY
            col <- colvec[pmin(100, as.integer(hh$counts[sel]/total*5000)+1)]
            panel.segments(hh$mids[sel], rep(y-yoff, n), hh$mids[sel],
                           rep(y-yoff-leny, n), col=col, lwd=3, lineend=2)
        }
    }
}





## ==========================================================================
## Plot a view object. For everything but gates we have the modified data
## available, hence we can plot directly
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("xyplot",
          signature(x="formula", data="view"),
          function(x, data, ...) xyplot(x, Data(data), ...))


setMethod("xyplot",
          signature(x="view", data="missing"),
          function(x, data, ...) xyplot(Data(x), ...))



## ==========================================================================
## Plot a gateView object. Essentially, this is calling the plot
## method on the parent data of the view, adding gates if appropriate.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("xyplot",
          signature(x="formula", data="gateView"),
          function(x, data, filter=NULL, par.settings, ...)
      {
          fres <-
              if(!is.null(filter)) filter else get(action(data)@filterResult)
          ## deparse the formula structure
          channel.y <- x[[2]]
          channel.x <- x[[3]]
          if (length(channel.x) == 3)
              channel.x <- channel.x[[2]]
          channel.x.name <- flowViz:::expr2char(channel.x)
          channel.y.name <- flowViz:::expr2char(channel.y)
          thisData <- Data(parent(data))
          if(!is.null(fres) && all(c(channel.x.name, channel.y.name) %in%
                                   unique(unlist(parameters(fres))))){
              l <- max(2, if(is(fres, "filterResultList"))
                       length(fres[[1]]) else length(fres))
              n <- if(is(fres, "filterResultList")) names(fres[[1]]) else
              names(fres)
              col <- rep("#00000030", l)
              names(col) <- n
              pop <- data@frEntry
              col[pop] <- "transparent"
              if(missing(par.settings))
                  par.settings <- list(gate=list(fill=col, col="#00000040"))
              xyplot(x, thisData, filter=fres, par.settings=par.settings, ...)
          }else{
              xyplot(x, thisData,  ...)
          }
      })



