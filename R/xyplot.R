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
#' @export 
#' @rdname xyplot
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
              time <- flowCore:::findTimeChannel(expr, strict = TRUE)
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
#' @export 
#' @rdname xyplot
prepanel.xyplot.flowframe.time <- 
    function(x, y, frame, time, xlim, ylim, ...)
{
    yc <- as.character(y)
    expr <- exprs(frame)  
    
    if(missing(xlim)){
      #use default calculation
      xx <- expr[, time]
      xlim <- range(xx, finite=TRUE)
    }
    
    if(missing(ylim)){
      #use default calculation
      yy <- expr[, yc]
      ylim <- unlist(range(frame, yc))
    }


    
    return(list(xlim=xlim, ylim=ylim))
}

## Panel function to do the timeline plotting with several different options
#' @param time A character string giving the name of the data column recording
#' time. If not provided, we try to guess from the available parameters.
#' @param nrpoints The number of points plotted on the smoothed plot in sparse
#'          regions. This is only listed here because we use a different default. See
#'              \code{\link[lattice:panel.smoothScatter]{panel.smoothScatter}} for details.
#' @param type type of rendering; see
#'                      \code{\link[lattice:panel.xyplot]{panel.xyplot}} for details. For the basic
#'                      \code{flowFrame} method without a detailed formula, the addtional type
#'                          \code{discrete} is available, which plots a smoothed average of the flow
#'                          cytometry values against time.
#' @param binSize The size of a bin (i.e., the number of events within a bin)
#' used for the smoothed average timeline plots.

#' @export 
#' @rdname xyplot
panel.xyplot.flowframe.time <- 
    function(x, y, frame, time, type="discrete", nrpoints=0, binSize=100, ...)
{
    y <- as.character(y)
    expr <- exprs(frame)
    xx <- expr[, time]
    yy <- expr[, y]
    ## We record in the state environment which type of plot we produce
#	plotType("gtime", c("Time", as.character(y)))
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



#' Methods implementing Lattice xyplots for flow data.
#' 
#' 
#' These functions create Trellis scatter plots (a.k.a. dot plots in the Flow
#' Cytometry community) from flow cytometry data.
#' 
#' 
#' The implementation of \code{xyplot} in \code{flowViz} is very close to the
#' original \code{lattice} version. Concepts like conditioning and the use of
#' panels apply directly to the flow cytometry data. The single fundamental
#' difference is that conditioning variables are not evaluated in the context
#' of the raw data, but rather in the \code{phenoData} slot environment (only
#' for the \code{flowSet} methods. Thus, we can directly condition on pheotypic
#' variables like sample groups, patients or treatments.
#' 
#' In the formula interface, the primary and secondary variables (separated by
#' the tilde) have to be valid parameter names. Please note that frequently
#' used variants like \code{FSC-H} and \code{SSC-H} are not syntactically
#' correct R symbols, and need to be wrapped in \code{` `}. E.g.,
#' \code{`FSC-H`}. For \code{flowSets}, the use of a conditioning variable is
#' optional. We implicitely condition on \code{flowFrames} and the default is
#' to arrange panels by sample names.
#' 
#' @name xyplot
#' @param x A formula describing the structure of the plot and the variables to
#' be used in the display.  In the \code{prepanel} and \code{panel} functions,
#' also the names of \code{\link[flowCore:flowFrame-class]{flowFrames}} or any
#' of the annotation data columns in the \code{phenoData} slot.
#' @param data,y,frame a \code{flowSet} or
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} object that serves as the
#' source of data. For the workflow methods, this can also be various
#' \code{\link[flowCore:view-class]{view}} or
#' \code{\link[flowCore:actionItem-class]{actionItem}} objects.
#' @param prepanel The prepanel function. See
#' \code{\link[lattice:xyplot]{xyplot}}.
#' @param panel The panel function. See \code{\link[lattice:xyplot]{xyplot}}.
#' @param xlab,ylab Labels for data axes, with suitable defaults taken from the
#' formula.

#' @param \dots
#' 
#' marker.only \code{logical} specifies whether to show both channel and marker
#' names
#' 
#' More arguments, usually passed on to the underlying lattice methods.
#' @section Methods:
#' 
#' \describe{
#' 
#' \item{xyplot}{\code{signature(x = "flowFrame", data = "missing")}: Creates
#' diagnostic time series plots of flow parameter values against time. These
#' plots are useful to detect quality issues in the raw data. If not provided
#' explicitely via the \code{tine} argument, the time parameter will be
#' automatically detected. The additional arguments \code{xlab}, \code{ylab},
#' \code{nrpoints}, and \code{layout} are only listed because \code{flowViz}
#' provides different defaults. Internally, they are directly passed on to the
#' underlying \code{lattice} functions}. Argument \code{type} can be a
#' combination of any of the types allowed in \code{lattice xyplots}, or
#' \code{discrete}, in which case a smoothed average of the parameter against
#' time is plotted. \code{binSize} controls the binning that is used for the
#' smoothing procedure.
#' 
#' \item{xyplot}{\code{signature(x = "formula", data = "flowFrame")}: Creates
#' scatter plots (a.k.a. dot plots) of a pair of FCM channels. Depending on the
#' setting of the \code{smooth} argument, the data will be rendered as a
#' partially smoothed density estimate (\code{smooth=TRUE}, the default) or as
#' a regular scatter plot with separate points for individual events. The
#' formula interface allows for fairly general plotting, however there are
#' certain limitations on the use of expressions as part of the formulae.
#' Unless you are sure about what you are doing, you should transform the raw
#' data in a separate step using one of the tools in the
#' \code{\link[flowCore:flowCore-package]{flowCore}} package rather than inline
#' using the formula interface. The method allows to superimpose gating results
#' though the \code{filter} argument. If \code{smooth=TRUE}, we try to add
#' spherical 2D representations of the gates if applicable. For
#' \code{smooth=FALSE}, gates are indicated by a grouping mechanism using
#' different point shapes or colors (unless \code{outline} is also \code{TRUE},
#' in which case the gate outlines are superimposed in addition to the
#' grouping). Argument \code{margins} controls how events on the margins of the
#' measurement range are treated. The default (\code{TRUE}) is to discard them
#' from any density estimation and later add them as separate glyphs. Argument \code{par.settings} can be used
#' to supply lists of graphical parameters. See \code{\link{flowViz.par.set}} for details on controlling graphical
#' parameters in these plots. }
#' 
#' \item{xyplot}{\code{signature(x = "formula", data = "flowSet")}: Scatter
#' plots from a \code{flowSet} object. We allow for conditioning on variables
#' in the \code{phenoData} slot of the \code{flowSet}. All additional arguments
#' that apply to the \code{flowFrame} method are also valid for
#' \code{flowSets}. }
#' 
#' }
#' @author F. Hahne, D. Sarkar
#' @seealso
#' 
#' Not all standard lattice arguments will have the intended effect, but many
#' should.  For a fuller description of possible arguments and their effects,
#' consult documentation on lattice.
#' @keywords methods dplot
#' @examples
#' 
#' 
#' data(GvHD)
#' GvHD <- GvHD[pData(GvHD)$Patient %in% 5:6]
#' 
#' ## a bivariate scatterplot
#' ## by default ('smooth=TRUE') panel.smoothScatter is used
#' xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], nbin = 100,
#' main="A single flowFrame")
#' 
#' ## A non-smooth version of the same data
#' xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], nbin = 100,
#' main="A single flowFrame", smooth=FALSE)
#' 
#' ## A non-smooth version of the same data with customerized color scheme
#' require(IDPmisc)
#' colramp <- colorRampPalette(IDPcolorRamp(21))
#' xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], nbin = 100,
#'        main="A single flowFrame", smooth=FALSE,
#'        colramp=colramp, pch=20, cex=0.1)
#' 
#' ## A hexbin version of non-smooth scatter plot  
#' xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], xbin = 128
#'        ,main="A single flowFrame", smooth=FALSE)
#' 
#' 
#' ## Visual artifacts created by the pileup of margin events
#' xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], nbin = 100,
#'        main="A single flowFrame", margin=FALSE)
#' 
#' 
#' ## simple bivariate scatter plot (a.k.a. dot plot)
#' ## for the whole flowSet, conditioning on Patient and
#' ## Visit
#' xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD)
#' 
#' ## Same bivariate scatter plot with replacing default color
#' require(IDPmisc)
#' cols <- colorRampPalette(IDPcolorRamp(21))
#' xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD, colramp=cols)
#' 
#' ## several examples with time on the X axis
#' ## first for a flowFrame
#' xyplot(GvHD[[1]])
#' 
#' ## and for flowSets
#' xyplot(`FSC-H` ~ Time | Visit, GvHD, 
#'        smooth = FALSE, type = "l", 
#'        subset = (Patient == 5), xbin = 32)
#' 
#' xyplot(`FSC-H` ~ Time | Patient+Visit, GvHD, 
#'        smooth = FALSE, type = "a",
#'        strip = FALSE, strip.left = TRUE,
#'        aspect = "xy", xbin = 32)
#' 
#' 
#' ## combine plots for two channels
#' ssc.time <- 
#' 
#'     xyplot(`SSC-H` ~ Time | factor(Patient):factor(Visit), GvHD, 
#'            smooth = FALSE, type = "a",
#'            strip = FALSE,
#'            strip.left = strip.custom(horizontal = TRUE),
#'            par.strip.text = list(lines = 3),
#'            between = list(y = rep(c(0, 0.5), c(6, 1))),
#'            scales = list(x = list(axs = "i"), y = list(draw = FALSE)),
#'            layout = c(1, 14), xbin = 32)
#' 
#' fsc.time <- 
#'     xyplot(`FSC-H` ~ Time | factor(Patient):factor(Visit), GvHD, 
#'            smooth = FALSE, type = "a",
#'            strip = FALSE,
#'            strip.left = strip.custom(horizontal = TRUE),
#'            par.strip.text = list(lines = 3),
#'            between = list(y = rep(c(0, 0.5), c(6, 1))),
#'            scales = list(x = list(axs = "i"), y = list(draw = FALSE)),
#'            layout = c(1, 14), xbin = 32)
#' 
#' plot(fsc.time, split = c(1, 1, 2, 1))
#' plot(ssc.time, split = c(2, 1, 2, 1), newpage = FALSE)
#' 
#' 
#' ## saving plots as variables allows more manipulation
#' plot(update(fsc.time[8:14], layout = c(1, 7)),
#'      split = c(1, 1, 1, 2))
#' 
#' plot(update(ssc.time[8:14], layout = c(1, 7)),
#'      split = c(1, 2, 1, 2), newpage = FALSE)
#' 
#' 
#' ## displaying filters
#' n2gate <- norm2Filter("SSC-H", "FSC-H")
#' 
#' xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD,
#'        filter=n2gate, subset=Patient==5)
#' 
#' xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit,
#'        data=transform("SSC-H"=asinh,"FSC-H"=asinh) %on% GvHD,
#'        smooth=FALSE, filter=n2gate, subset=Patient == 5, xbin = 32)
#' 
#' 
#' ## displaying filters with stats
#' n2gate.results <- filter(GvHD, n2gate)
#' 
#' xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD,
#'        subset=Patient == "6",
#'        filter=n2gate.results, smooth=FALSE, xbin = 32
#'        ,stats=TRUE
#'        ,abs=TRUE
#'        ,digits=3 
#'        )
#' 
#'        
#' ## displaying multiple filters in one panel with stats
#' recGate1<-rectangleGate("FL3-H"=c(2.3,4.1),"FL2-H"=c(6.8,9))
#' recGate2<-rectangleGate("FL3-H"=c(1,3),"FL2-H"=c(4,6))
#' filters1<-filters(list(recGate1,recGate2))
#' trans<-transform("FL2-H"=asinh,"FL3-H"=asinh)
#' trans_data<-transform(GvHD[1:2],trans)
#' #replicate filters object across samples
#' flist <- list(filters1 , filters1)
#' names(flist) <- sampleNames(trans_data)
#' xyplot(`FL2-H` ~ `FL3-H`
#' 	   ,data=trans_data 
#'        ,filter= flist
#'        ,stats=TRUE
#'        ,margin=FALSE
#'        , xbin = 32
#'        , smooth = FALSE
#'        )
#' 
#' #display recGate2 as a overlay 
#' overlay <- Subset(trans_data,recGate1)
#' xyplot(`FL2-H` ~ `FL3-H`
#' 	   ,data=trans_data 
#'        ,filter=recGate2
#'        ,stats=TRUE
#'        ,margin=FALSE
#'        , smooth = FALSE
#'        , xbin = 32
#'        ,overlay= list(rect2 = overlay)
#'        ,par.settings = list(overlay.symbol = list(cex = 0.1))
#'        )
#' 
#' 
#' 
#' 
#' @export
#' @aliases xyplot,formula,flowFrame-method
setMethod("xyplot",
    signature=signature(x="formula",
        data="flowFrame"),
    definition=function(x,
                        data
                        , filter = NULL
                        , overlay= NULL#a list of flowFrames
                        , stats = FALSE
                        , strip.text = NULL#remove strip by setting cond as NULL
                        , ...)
{
  
  if(is.null(strip.text))   {
    defaultCond <- NULL
  }else
  {
    #try to use strip.text as the dummy sample name if it is given 
    #in order to display the customized text in strip
    defaultCond <- "name"
    data <- list(data)
    names(data) <- strip.text
  }
  
  data <- as(data, "flowSet")
  
  sn <- sampleNames(data)
  if(!is.null(filter)){
    filter <- list(filter)
    names(filter) <- sn
  }
  
  overlay <- .process_flowFrame_overlay(overlay, sn)
  
  if(!is.null(stats)){
    stats <- list(stats)
    names(stats) <- sn
  }
  
    
  xyplot(x, data, filter = filter
          , overlay = overlay, stats = stats
          , defaultCond = defaultCond , ...)  
})

#' convert a single flowFrame or a list of flowFrames to a list of flowSet
#' @param overlay flowFrame or a list of flowFrame objects
#' @param sn sample name
.process_flowFrame_overlay <- function(overlay, sn){
  if(!is.null(overlay)){
    # convert single flowFrame to a list
    if(!is.list(overlay))
      overlay <- list(dummy = overlay) #panel.xyplot.flowFrame expect a named list
    # convert a list of flowFrames to a list of flowSets
    overlay <- lapply(overlay, function(thisOverlay){
          thisOverlay <- as(list(thisOverlay), "flowSet")
          sampleNames(thisOverlay) <- sn
          thisOverlay
        })
    
  }
  overlay
}

## Prepanel function to set up dimensions. We want to use the instrument measurement
## range instead of the absolute range of the data. We also record the data ranges
## in the internal state environment for further use.
#' @param xlim,ylim limits for data axes. If not given, they are taken from the
#' ranges stored in flowFrame
#' @export 
#' @rdname xyplot
prepanel.xyplot.flowframe <- 
    function(frame, channel.x.name, channel.y.name, x, y, xlim, ylim,  ...)
{
    ranges <- range(frame)
#    browser()
    if(missing(xlim)){
      #use default calculation
      xlim <- if(!length(grep("`", channel.x.name, fixed=TRUE))){
          tmp <- ranges[, channel.x.name]
          xd <- diff(tmp)/15
          tmp + c(-1,1)*xd
        }else NULL
    }
    
    if(missing(ylim)){
      #use default calculation
        ylim <- if(!length(grep("`", channel.y.name, fixed=TRUE))){
          tmp <- ranges[, channel.y.name]
          yd <- diff(tmp)/15
          tmp + c(-1,1)*yd
        }else NULL
    }
    
    
    return(list(xlim=xlim, ylim=ylim))
}

#old version of panel function that still expects x,y vectors
panel.xyplot.flowframe.old <- function(x,y,frame,
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
    gp
    ,xbins=0
    ,binTrans=sqrt
    ,stats = FALSE
    ,pos=0.5
    ,digits=2
    ,abs=FALSE
    ,overlay.x=NULL
    ,overlay.y=NULL
    ,...)
{
  
  
  ## graphical parameter defaults
  limits<-prepanel.xyplot.flowframe(frame,channel.x.name,channel.y.name)
  xlim<-limits$xlim
  ylim<-limits$ylim
  argcolramp <- list(...)$colramp
  gpar <- flowViz.par.get()
  parameters<-c(channel.x.name, channel.y.name)
  if(!is.null(gp))
    gpar <- lattice:::updateList(gpar, gp)
  if(is.null(gpar$gate$cex))
    gpar$gate$cex <- cex
  if(is.null(gpar$gate$pch))
    gpar$gate$pch <- pch
#   browser()
  if(is.null(gpar$gate$plotType))
    gpar$gate$plotType<-"l"
  if(is.null(gpar$density))
    gpar$density<-TRUE
  ## Whenever we have a function call in the formula we might no longer be the
  ## original scale in which the gate was defined, and we have no clue how to
  # plot it
  validName <- !(length(grep("\\(", channel.x.name)) ||
        length(grep("\\(", channel.y.name)))
  
  ## in order to display overlay ,smooth needs to be set as TRUE
  if(!is.null(overlay.x))
    smooth<-TRUE
#browser()
  if(nrow(frame)==0)
    return (NULL)
  ## We remove margin events before passing on the data to panel.smoothScatter
  ## and after plotting indicate those events by grayscale lines on the plot
  ## margins
  l <- length(x)
  
  
  if (smooth){
    if(margin){
      r <- range(frame, c(channel.x.name, channel.y.name))
#               l <- length(x)
      inc <- apply(r, 2, diff)/1e5
      dots <- list(...)
      nb <- if("nbin" %in% names(dots)) rep(dots$nbin, 2) else rep(64, 2)
      selxL <- x > r[2,channel.x.name]-inc[1]
      selxS <- x < r[1,channel.x.name]+inc[1]
      selyL <- y > r[2,channel.y.name]-inc[2]
      selyS <- y < r[1,channel.y.name]+inc[2]
      allsel <- !(selxL | selxS | selyL | selyS)
      
      if(sum(allsel)>0)
      {
        panel.smoothScatter(x[allsel], y[allsel],
            range.x=list(r[,1], r[,2]), ...)
        addMargin(r[1,channel.x.name], y[selxS], r, l, nb)
        addMargin(r[2,channel.x.name], y[selxL], r, l, nb, b=TRUE)
        addMargin(x[selyS], r[1,channel.y.name], r, l, nb)
        addMargin(x[selyL], r[2,channel.y.name], r, l, nb, b=TRUE)
      }
      else
      {
        panel.smoothScatter(x, y, ...)
      }
    }else{
      panel.smoothScatter(x, y, ...)
    }
    ptList<-plotType("gsmooth", c(channel.x.name, channel.y.name))
  }else{
    #for non-smoothed plot:
    #1:always remove boundary events for hexbin version 
    #since they are going to affect the color encoding for density
    #2.for non-hexbin version,when gpar$density=FALSE,which means one-color non-densityscatter plot
    #we then keep the boundary events 
    if(margin){
      
      r <- range(frame, c(channel.x.name, channel.y.name))
#               l <- length(x)
      inc <- apply(r, 2, diff)/1e5
      dots <- list(...)
      nb <- if("nbin" %in% names(dots)) rep(dots$nbin, 2) else rep(64, 2)
      selxL <- x > r[2,channel.x.name]-inc[1]
      selxS <- x < r[1,channel.x.name]+inc[1]
      selyL <- y > r[2,channel.y.name]-inc[2]
      selyS <- y < r[1,channel.y.name]+inc[2]
      allsel <- !(selxL | selxS | selyL | selyS)
      #we may want to skip marginal events in non-smoothed version to save time           
      if(sum(allsel)>0)
      {
        addMargin(r[1,channel.x.name], y[selxS], r, l, nb)
        addMargin(r[2,channel.x.name], y[selxL], r, l, nb, b=TRUE)
        addMargin(x[selyS], r[1,channel.y.name], r, l, nb)
        addMargin(x[selyL], r[2,channel.y.name], r, l, nb, b=TRUE)
      }
      
      x<-x[allsel]
      y<-y[allsel]
      
    }
#           browser()
    if(xbins>0)
    {
      #using hexbin package to do the hexagon plot    
      bin<-hexbin(x,y,xbins=xbins)
      if (is.null(argcolramp))
        argcolramp<-flowViz.par.get("argcolramp")
      #           if (is.null(argcolramp))
      #               argcolramp<-colorRampPalette(c("blue","green","yellow","red"),bias=1)
      #           browser()
      grid.hexagons(bin,colramp = argcolramp,trans=binTrans)      
      
      
    }else
    {
      if (is.null(argcolramp))
        argcolramp<-flowViz.par.get("argcolramp")
      if(gpar$density)
        col <- densCols(x, y, colramp=argcolramp)
      panel.xyplot(x, y, col=col, cex=cex, pch=pch, alpha=alpha, ...)
      
    }
    
    ptList<-plotType("gpoints", c(channel.x.name, channel.y.name))
  }
#   }
#    browser()
  #plot gate
  if(!is.null(filter) && validName){
    if(is(filter,"filters"))
    {
      
      mapply(filter,stats,FUN=function(curFilter,curStats){
            
            
            if(gpar$gate$plotType=="p")##highlight the dots within gate,by default it is now disabled 
            {
              if(!is(curFilter, "filterResult"))
                curFilter <- filter(frame, curFilter)
              rest <- Subset(frame, !filter)
              x <- exprs(rest[,channel.x.name])
              y <- exprs(rest[,channel.y.name])
              
              
              glpoints(curFilter, frame,
                  channels=c(channel.x.name, channel.y.name),
                  verbose=FALSE, gpar=gpar, strict=FALSE,ptList=ptList, ...)
              if(outline)
                glpolygon(curFilter, frame,
                    channels=c(channel.x.name, channel.y.name),
                    verbose=FALSE, gpar=gpar, names=FALSE,
                    strict=FALSE,ptList=ptList,xlim=xlim
                    ,ylim=ylim)
              
            }else
            {
              
              names <- .getStats(curFilter,curStats, frame, digits, ...)
              
              glpolygon(curFilter, frame,
                  
                  verbose=FALSE, gpar=gpar
                  , names=names
                  ,strict=FALSE
                  ,pos=pos
                  ,abs=abs
                  ,ptList=ptList
                  ,xlim=xlim
                  ,ylim=ylim
              )
            }
          })
    }else
    {
      if(gpar$gate$plotType=="p")##highlight the dots within gate,by default it is now disabled 
      {
        if(!is(filter, "filterResult"))
          filter <- filter(frame, filter)
        rest <- Subset(frame, !filter)
        x <- exprs(rest[,channel.x.name])
        y <- exprs(rest[,channel.y.name])
        
        
        glpoints(filter, frame,
            channels=c(channel.x.name, channel.y.name),
            verbose=FALSE, gpar=gpar, strict=FALSE,ptList=ptList, ...)
        if(outline)
          glpolygon(filter, frame,
              channels=c(channel.x.name, channel.y.name),
              verbose=FALSE, gpar=gpar, names=FALSE,
              strict=FALSE,ptList=ptList,xlim=xlim
              ,ylim=ylim)
        
      }else
      {
        names <- .getStats(filter,stats, frame, digits, ...) 
        
#               browser()
        glpolygon(filter, frame,
#                       channels=c(channel.x.name, channel.y.name),
            verbose=FALSE, gpar=gpar
            , names=names
            ,strict=FALSE
            ,pos=pos
            ,abs=abs
            ,ptList=ptList
            ,xlim=xlim
            ,ylim=ylim
        )
      }
    }
  }
  
  if(!is.null(overlay.x)&&!is.null(overlay.y))
  {
    lpoints(overlay.x,overlay.y,col="red",cex=cex*3,pch=pch)
    #plot stats for bool gates
    if(stats&&is.null(filter))
    {
      p.stats<-length(overlay.x)/l
#           p.stats<-sprintf(paste("%.",prec,"f%%",sep=""),p.stats*100)
      p.stats<-paste(format(p.stats*100,digits=digits),"%",sep="")
#           browser()
      xx<-xlim
      yy<-ylim
      
      pos <- rep(pos, length=2)[1:2]
      xx<-xx[1]+diff(xx)*pos[1]
      yy<-yy[1]+diff(yy)*pos[2]
      
      gltext(xx, yy, labels=p.stats, adj=0.5, gp=gpar$gate.text)
    }
  }
  
  
}


## Panel function that allows us to add filters on the plot. The actual plotting
## is done by either the panel.smoothScatter or the default lattice panel.xyplot
## function
##when xbins>0, we do the hexagon plot provided by hexbin package to improve the speed
#overlay is a list of flowFrames, which are the extra points need to be plotted on top of the x,y



#' @param layout These arguments are passed unchanged to the corresponding
#' methods in lattice, and are listed here only because they provide different
#' defaults.  See documentation for the original methods for details.
#' @param overlay The extra cell events plotted on top of the current cell
#' population.  It is a \code{flowSet} for \code{panel.xyplot.flowset} function
#' and a \code{flowFrame} for \code{xyplot(c("formula","flowFrame"))} method.
#' @param channel.x.name,channel.y.name Character strings giving corresponding
#' names used to match filter parameters if applicable.
#' @param smooth Logical. If \code{TRUE}, \code{panel.smoothScatter} is used to
#' display a partially smoothed version of the data.  Otherwise, events are
#' plotted individually, as in a standard scatter plot. If \code{FALSE}, a
#' graphical parameter \code{colramp} can be used to obtain a coloring of
#' points that is indicative of their local density.
#' @param filter A \code{\link[flowCore:filter-class]{filter}},
#' \code{\link[flowCore:filterResult-class]{filterResult}} or
#' \code{\link[flowCore:filterResult-class]{filterResultList}} object or a list
#' of such objects of the same length as the \code{flowSet}.  Also a
#' \code{\link[flowCore:filters-class]{filters}} or A
#' \code{\link[flowCore:filtersList-class]{filtersList}} can be passed to
#' xyplot in order to plot multiple filters/gates(with the same x,y parameters)
#' on one panel to represent multiple sub-populations.  The appropriate
#' spherical 2D representation of this filter will be superimposed on the plot
#' if \code{smooth=TRUE}, or the result of the filtering operation will be
#' indicated by grouping if \code{smooth=FALSE}. The software will figure out
#' whether the \code{filter} needs to be evaluated in order to be plotted (in
#' which case providing a \code{filterResult} can speed things up
#' considerably).
#' @param margin Logical indicating whether to truncate the density estimation
#' on the margins of the measurement range and plot margin events as lines if
#' \code{smooth=TRUE}. To avoid visual artifacts it is highly recommended to
#' set this option to \code{TRUE}.
#' @param outline Logical, specifying whether to add the boundaries of a gate
#' to the plot when \code{smooth=FALSE} in addition to the grouping. Defaults
#' to \code{FALSE}.
#' @param pch,cex,col,alpha Graphical parameters used when \code{smooth=FALSE}.
#' These mostly exist for conveniance and much more control is available
#' throught the \code{lattice}-like \code{par.setting} and
#' \code{flowViz.par.set} customization. See \code{\link{flowViz.par.set}} for
#' details.
#' @param gp A list of graphical parameters that are passed down to the low
#' level panel functions. This is for internal use only. The public user
#' interface to set graphical parameters is either \code{par.settings} for
#' customization of a single call or \code{flowViz.par.set} for customization
#' of session-wide defaults.
#' @param channel.x,channel.y Expressions defining the x and y variables in
#' terms of columns in the data.  Can involve functions or multiple columns
#' from the data, however this usage is discouraged.
#' @param frames An environment containing frame-specific data.
#' @param xbins The argument passed to \code{\link[hexbin:hexbin]{hexbin}}
#' ,which is the number of bins partitioning the range of xbnds.  It is set as
#' 0 by default,which plots all the events without binning.  When it is larger
#' than 0,hexbin plot engine is used for the faster plotting.  Note that it is
#' only valid when \code{smooth} is set as FALSE .
#' @param binTrans The argument passed to
#' \code{\link[hexbin:grid.hexagons]{grid.hexagons}} ,which is a transformation
#' function (or NULL) for the count.  It is \code{\link[base:sqrt]{sqrt}} by
#' default.
#' @param stats,pos,digits,abs Arguments to control statistics that is
#' associated with \code{\link[flowCore:filter-class]{filter}} to be plotted
#' Currently only population proportion/percentage is supported.  \code{stats}
#' is a \code{logical} scalar indicating whether to display statistics. Default
#' is FALSE.  \code{pos} is the \code{numeric} scalar (range within c(0,1)) or
#' vector(length of 2,first is for x-axis,second for y-axis) to control the
#' position of the statistics label. It is set as 0.5,which is the center.
#' \code{digits} is an \code{integer} indicating the number of significant
#' digits to be used when displaying the percentage of population
#' statistics,Default is 2. see more details from \code{\link[base]{format}}
#' \code{abs} is a \code{logical} scalar indicating whether the \code{pos} is
#' relative to the gate boundary or the entire xy-axis(absolute position).  By
#' default it is set as FALSE,which indicates the position is relative to gate.
#' @param overlay.symbol list of the lattice graphic parameters to format the
#' overlay points.
#' @param sample.ratio \code{numeric} the ratio of sub-sampling of events to
#' speed up plotting.
#' @param strip.text A \code{character} that customizes the text in strip.
#' Default is NULL, which does not display the strip box at all. It is only
#' valid when plotting a \code{flowFrame}
#' @param checkName \code{logical} indicating whether to skip checking the
#' bracket '(' in channel name
#' @export panel.xyplot.flowframe
#' @importFrom hexbin hexbin grid.hexagons
#' @importFrom IDPmisc IDPcolorRamp col2hsv
#' @export 
#' @rdname xyplot
panel.xyplot.flowframe <- function(frame,
                                   filter=NULL,
                                   smooth=TRUE,
                                   margin=TRUE,
                                   outline=FALSE,
                                   channel.x.name,
                                   channel.y.name,
                                   pch=gp$flow.symbol$pch,
                                   alpha=gp$flow.symbol$alpha,
                                   cex=gp$flow.symbol$cex,
                                   col=gp$flow.symbol$col,
                                   gp
								   ,xbins=0
						   		   ,binTrans=sqrt
						   			,stats = FALSE
									,pos=0.5
									,digits=2
									,abs=FALSE
									,overlay = NULL
                                    ,checkName = TRUE
                                    ,sample.ratio = 1
                                    , overlay.symbol = NULL
						   			,...)
{
    
	
    
    
    ## graphical parameter defaults
	limits<-prepanel.xyplot.flowframe(frame,channel.x.name,channel.y.name)
	xlim<-limits$xlim
	ylim<-limits$ylim
    argcolramp <- list(...)$colramp
	
	parameters<-c(channel.x.name, channel.y.name)
    if(is.null(gp$gate$cex))
        gp$gate$cex <- cex
    if(is.null(gp$gate$pch))
        gp$gate$pch <- pch
#	browser()
	if(is.null(gp$gate$plotType))
		gp$gate$plotType<-"l"
	if(is.null(gp$density))
		gp$density<-TRUE
    ## Whenever we have a function call in the formula we might no longer be the
    ## original scale in which the gate was defined, and we have no clue how to
    # plot it
    if(checkName)
      validName <- !(length(grep("\\(", channel.x.name)) ||
                   length(grep("\\(", channel.y.name)))
    else
      validName <- TRUE
	
    if(!validName)
      warning("Gate will not be plotted because channel names contain '(' character! Try to set checkName to FALSE to skip this check.")

    ## in order to display overlay ,smooth needs to be set as TRUE
#    if(!is.null(overlay.x))
#         smooth<-TRUE
    
	if(nrow(frame)==0)
		return (NULL)
    ## We remove margin events before passing on the data to panel.smoothScatter
    ## and after plotting indicate those events by grayscale lines on the plot
    ## margins
    dat <- exprs(frame)
    l <- nrow(dat)
    if(sample.ratio < 1 && sample.ratio > 0){
      ind <- sample.int(l, size = l * sample.ratio)
      x <- dat[ind, channel.x.name]
      y <- dat[ind, channel.y.name]  
    }else
    {
      x <- dat[, channel.x.name]
      y <- dat[, channel.y.name]
    }
      
    
		
	    if (smooth){
			if(margin){
				r <- range(frame, c(channel.x.name, channel.y.name))
#				l <- length(x)
				inc <- apply(r, 2, diff)/1e5
				dots <- list(...)
				nb <- if("nbin" %in% names(dots)) rep(dots$nbin, 2) else rep(64, 2)
				selxL <- x > r[2,channel.x.name]-inc[1]
				selxS <- x < r[1,channel.x.name]+inc[1]
				selyL <- y > r[2,channel.y.name]-inc[2]
				selyS <- y < r[1,channel.y.name]+inc[2]
				allsel <- !(selxL | selxS | selyL | selyS)
				
	            if(sum(allsel)>0)
	            {
	                panel.smoothScatter(x[allsel], y[allsel],
	                                    range.x=list(r[,1], r[,2]), ...)
	                addMargin(r[1,channel.x.name], y[selxS], r, l, nb)
	                addMargin(r[2,channel.x.name], y[selxL], r, l, nb, b=TRUE)
	                addMargin(x[selyS], r[1,channel.y.name], r, l, nb)
	                addMargin(x[selyL], r[2,channel.y.name], r, l, nb, b=TRUE)
	            }
	            else
	            {
	                panel.smoothScatter(x, y, ...)
	            }
	        }else{
	            panel.smoothScatter(x, y, ...)
	        }
			ptList<-plotType("gsmooth", c(channel.x.name, channel.y.name))
	    }else{
			#for non-smoothed plot:
			#1:always remove boundary events for hexbin version 
			#since they are going to affect the color encoding for density
			#2.for non-hexbin version,when gp$density=FALSE,which means one-color non-densityscatter plot
			#we then keep the boundary events 
			if(margin){
				
				r <- range(frame, c(channel.x.name, channel.y.name))
#				l <- length(x)
				inc <- apply(r, 2, diff)/1e5
				dots <- list(...)
				nb <- if("nbin" %in% names(dots)) rep(dots$nbin, 2) else rep(64, 2)
				selxL <- x > r[2,channel.x.name]-inc[1]
				selxS <- x < r[1,channel.x.name]+inc[1]
				selyL <- y > r[2,channel.y.name]-inc[2]
				selyS <- y < r[1,channel.y.name]+inc[2]
				allsel <- !(selxL | selxS | selyL | selyS)
				#we may want to skip marginal events in non-smoothed version to save time			
				if(sum(allsel)>0)
				{
					addMargin(r[1,channel.x.name], y[selxS], r, l, nb)
					addMargin(r[2,channel.x.name], y[selxL], r, l, nb, b=TRUE)
					addMargin(x[selyS], r[1,channel.y.name], r, l, nb)
					addMargin(x[selyL], r[2,channel.y.name], r, l, nb, b=TRUE)
				}
					
				x<-x[allsel]
				y<-y[allsel]
				
			}
#			browser()
        
            if (is.null(argcolramp))
              argcolramp <- flowViz.par.get("argcolramp")
            if(!is.null(overlay)){
              argcolramp <- .colRmpPlt(alpha = gp$overlay.symbol$bg.alpha)
            }
              
			if(xbins > 0)
			{
              
				#using hexbin package to do the hexagon plot	
				bin <- hexbin(x,y,xbins=xbins)
                
				gridRes <- try(grid.hexagons(bin,colramp = argcolramp, trans=binTrans, border =0), silent = TRUE)		
                
                if(class(gridRes) == "try-error"){
                  #if error then try without trans
                  grid.hexagons(bin,colramp = argcolramp, border =0)
                }				
                
			}else
			{
				
				if(gp$density)
					col <- densCols(x, y, colramp=argcolramp)
#                browser()
                dots <- list(...)
                dots$darg <- NULL
                dots$type <- NULL
                do.call(panel.xyplot, args = c(list(x = x, y = y
                                                , col = col
                                                , cex = cex
                                                , pch = pch
                                                , alpha = alpha
                                                )
                                                , dots
                                                )
                                                
                                            )
                    	
			}
            
			ptList<-plotType("gpoints", c(channel.x.name, channel.y.name))
		}

#    browser()
	#plot gate
	if(!is.null(filter) && validName){
		if(is(filter,"filters"))
		{

			mapply(filter,stats,FUN=function(curFilter,curStats){
						
						
						if(gp$gate$plotType=="p")##highlight the dots within gate,by default it is now disabled 
						{
							if(!is(curFilter, "filterResult"))
								curFilter <- filter(frame, curFilter)
							rest <- Subset(frame, !filter)
							x <- exprs(rest[,channel.x.name])
							y <- exprs(rest[,channel.y.name])
							
							
							glpoints(curFilter, frame,
									channels=c(channel.x.name, channel.y.name),
									verbose=FALSE, gpar=gp, strict=FALSE,ptList=ptList, ...)
							if(outline)
								glpolygon(curFilter, frame,
										channels=c(channel.x.name, channel.y.name),
										verbose=FALSE, gpar=gpar, names=FALSE,
										strict=FALSE,ptList=ptList,xlim=xlim
										,ylim=ylim)
							
						}else
						{
                          
                            names <- .getStats(curFilter,curStats, frame, digits, ...)
							
							glpolygon(curFilter, frame,

									verbose=FALSE, gpar=gp
									, names=names
									,strict=FALSE
									,pos=pos
									,abs=abs
									,ptList=ptList
									,xlim=xlim
									,ylim=ylim
							)
						}
					})
		}else
		{
			if(gp$gate$plotType=="p")##highlight the dots within gate,by default it is now disabled 
			{
				if(!is(filter, "filterResult"))
					filter <- filter(frame, filter)
				rest <- Subset(frame, !filter)
				x <- exprs(rest[,channel.x.name])
				y <- exprs(rest[,channel.y.name])
				
				
				glpoints(filter, frame,
						channels=c(channel.x.name, channel.y.name),
						verbose=FALSE, gpar=gp, strict=FALSE,ptList=ptList, ...)
				if(outline)
					glpolygon(filter, frame,
							channels=c(channel.x.name, channel.y.name),
							verbose=FALSE, gpar=gp, names=FALSE,
							strict=FALSE,ptList=ptList,xlim=xlim
							,ylim=ylim)
				
			}else
			{
              names <- .getStats(filter,stats, frame, digits, ...) 
                  
#				browser()
				glpolygon(filter, frame,
#						channels=c(channel.x.name, channel.y.name),
						verbose=FALSE, gpar=gp
						, names=names
						,strict=FALSE
						,pos=pos
						,abs=abs
						,ptList=ptList
						,xlim=xlim
						,ylim=ylim
				)
			}
		}
	}
#	browser()
	if(!is.null(overlay))
	{
        overlayNames <- names(overlay)
        for(overLayName in overlayNames){
          thisOverlay <- overlay[[overLayName]]
          thisDat <- exprs(thisOverlay)
          overlay.x <- thisDat[, channel.x.name]
          overlay.y <- thisDat[, channel.y.name]
          
          this.overlay.symbol <- gp$overlay.symbol#default setting
          user.overlay.symbol <- overlay.symbol[[overLayName]]#customized settings
          #update the default with customized settings
          if(!is.null(user.overlay.symbol))
            this.overlay.symbol <- lattice:::updateList(this.overlay.symbol, overlay.symbol[[overLayName]])
          lpoints(overlay.x, overlay.y
                    , col = this.overlay.symbol[["fill"]]
                    , alpha = this.overlay.symbol[["alpha"]]
                    , cex= this.overlay.symbol[["cex"]] 
                    , pch = this.overlay.symbol[["pch"]]
                    )
        }
        #plot stats for bool gates
        if(is.null(filter))
        {
          if(is.logical(stats)&&stats)
          {
            p.stats<-length(overlay.x)/l
            
            p.stats<-paste(format(p.stats*100,digits=digits),"%",sep="")
            #			browser()
            xx<-xlim
            yy<-ylim
            
            pos <- rep(pos, length=2)[1:2]
            xx<-xx[1]+diff(xx)*pos[1]
            yy<-yy[1]+diff(yy)*pos[2]
            
            gltext(xx, yy, labels=p.stats, adj=0.5, gp=gp$gate.text)
          }
        }
	}
		
	
}




##############################################################################
##                            flowSet methods                               ##
##############################################################################
## xyplot method for flowSets with formula.
#' @export 
#' @rdname xyplot
#' @param par.settings A list of lists of graphical parameters.  See \code{\link{flowViz.par.set}} for details.
setMethod("xyplot",
          signature=signature(x="formula",
                              data="flowSet"),
          definition= function(x,data, ...){
            thisTrellisObj <- .xyplot.flowSet(x, data, plotType = "xyplot", ...)
            thisData <- thisTrellisObj[["panel.args.common"]][["frames"]]
            thisTrellisObj[["panel.args.common"]][["frames"]] <- thisData@frames
            thisTrellisObj
          })
      
## flowViz:::.xyplot.flowSet now passes data instead of data@frames 
## within flowViz::xyplot method that changes it back to data@frames
## however ncdfFlow::xyplot keeps data as it is
#'  \item{marker.only}{ \code{ligcal} specifies whether to show both channel and marker names }
.xyplot.flowSet <- function(x,
                              data,
                              smooth=TRUE,
                              filter=NULL,
                              as.table=TRUE,
                              prepanel=prepanel.xyplot.flowset,
                              panel=panel.xyplot.flowset
                              , xlab= NULL 
                              , ylab= NULL 
                              , par.settings=NULL
                              , axis= axis.grid
                              ,defaultCond = "name" #to override the default conditional variable 'name'
                                            #mainly used for plotting single flowFrame
                              , between = list(x=0.2,y=0.2)
                              , plotType = "xyplot"
                              , marker.only = FALSE
                              , ...)
      {
       
          ## no conditioning variable, we chose 'name' as default
          if (length(x[[3]]) == 1){
              if(!is.null(defaultCond)){
                tmp <- x[[3]]
                strFormula <- "~dummy"  
                strFormula <- paste(strFormula,defaultCond,sep = "|")
                strFormula <-as.formula(strFormula)
                x[[3]] <- (strFormula)[[2]]
                x[[3]][[2]] <- tmp
              }
          }
          if(! "name" %in% names(pData(data)))
              pData(data)$name <- sampleNames(data)
          ## par.settings will not be passed on to the panel functions, so
          ## we have to fetch it from ... and stick the gate relevant stuff
          ## back it in there manually
           
		  this.par.settings <- par.settings #copy the customized settings
		  gpar <- flowViz.par.get()# grab the default theme
		  
          #update the customized settings 
		  if(!is.null(this.par.settings))
		  {
            this.par.settings <- lattice:::updateList(gpar, this.par.settings)
              
		  }else
            this.par.settings <- gpar
            
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
          
          
          frm <- data[[1, use.exprs = FALSE]]
          xObj <- getChannelMarker(frm, channel.x.name)
          yObj <- getChannelMarker(frm, channel.y.name)
          channel.x.name <- xObj[["name"]]
          channel.y.name <- yObj[["name"]]
          if(marker.only){
            default_xlab <- as.character(ifelse(is.na(xObj[,"desc"]), channel.x.name, xObj[,"desc"]))
            default_ylab <- as.character(ifelse(is.na(yObj[,"desc"]), channel.y.name, yObj[,"desc"]))
            
          }else
          {
            default_xlab <- sub("NA","",paste(unlist(xObj),collapse=" "))
            default_ylab <- sub("NA","",paste(unlist(yObj),collapse=" "))
          }
         
#          browser()
          if(is.null(xlab)){
            xlab <- default_xlab
          }else
          {
            if(is.list(xlab))
              xlab <- lattice:::updateList(list(label = default_xlab), xlab) #update default scales if non-null scales are specified
          }
          
          if(is.null(ylab)){
            ylab <- default_ylab
          }else
          {
            if(is.list(ylab))
              ylab <- lattice:::updateList(list(label = default_ylab), ylab) #update default scales if non-null scales are specified
          }
          
          channel.x <- as.expression(channel.x)
          channel.y <- as.expression(channel.y)
          ## use densityplot or bwplot method with dedicated panel and prepanel
          ## functions to do the actual plotting
          
          if(plotType == "xyplot"){
            densityplot(x, data=pd, prepanel=prepanel, panel=panel,
                frames=data, channel.x=channel.x,
                channel.y=channel.y, channel.x.name=channel.x.name,
                channel.y.name=channel.y.name, xlab=xlab, ylab=ylab,
                smooth=smooth, gp=this.par.settings, as.table=as.table, filter=filter,
                par.settings=this.par.settings, axis = axis, between = between, ...)
          }else if(plotType %in% c("densityplot", "histogram")){     #use bwplot instead due to the conflicts with the default darg argument of lattice::densityplot      
            bwplot(x, data=pd, prepanel=prepanel, panel=panel,
                      frames=data, channel.x=channel.x,
                      channel.y=channel.y, channel.x.name=channel.x.name,
                      channel.y.name=channel.y.name, xlab=xlab, ylab=ylab,
                      smooth=smooth, gp=this.par.settings, as.table=as.table, filter=filter,
                      par.settings=this.par.settings, axis = axis, between = between, plotType = plotType, ...)
        }else
          stop("unknown plot type: ", plotType)
      }



## Prepanel function to set up dimensions. We want to use the instrument measurement
## range instead of the absolute range of the data. We also record the data ranges
## in the internal state environment for further use.
#' @export 
#' @rdname xyplot
prepanel.xyplot.flowset <- 
    function(x, frames, channel.x.name, channel.y.name
      , xlim , ylim, ...)
{
  if (length(nm <- as.character(x)) > 1)
    stop("must have only one flow frame per panel!Please ensure the conditional variable is not set properly so that each group just has one sample. ")

    if (length(nm) == 1){
        ranges <- range(frames[[nm, use.exprs = FALSE]])
        if(missing(xlim)){
          #use default calculation
          xlim <- if(!length(grep("`", channel.x.name, fixed=TRUE))){
              tmp <- ranges[, channel.x.name]
              xd <- diff(tmp)/15
              tmp + c(-1,1)*xd
          }else NULL
      }
      if(missing(ylim)){
        #use default calculation
        ylim <- if(!length(grep("`", channel.y.name, fixed=TRUE))){
            tmp <- ranges[, channel.y.name]
            yd <- diff(tmp)/15
            tmp + c(-1,1)*yd
        }else NULL
      }
        return(list(xlim=xlim, ylim=ylim))
    }else
		return(list())
}

## 'filter' either has to be a single filter, or a list of filters matching
## the flowSet's sample names, or a filterResultList.
#' @param nm \code{character} sample names
#' @param lc \code{numeric} number of channels, only used for stacked densityplot validity check, ignored for xyplot and non-stacked densityplot
#' @param which.channel \code{numeric} index of channel, only used for stacked densityplot validity check, ignored for xyplot and non-stacked densityplot
.processFilter <- function(filter, nm, lc = 0, which.channel = 0){
#  browser()
  if(!is.null(filter)){

    if(lc > 1){
      if(!is.list(filter) && is(filter, "filterResultList") ||
          length(filter) != lc)
        stop("You are plotting several channels. Argument 'filter'/n",
            "must be a list of the same length as number of channels.",
            call.=FALSE)
      filter <- filter[[which.channel]]
    }
    
    if(!is.list(filter)){
      if(is(filter, "filter")){
        filter <- lapply(seq_along(nm), function(x) filter)
        names(filter) <- nm
      }
    }else if(!is(filter, "filterResultList")&&!is(filter, "filtersList"))
      filter <- as(filter, "filterResultList")
    
    
    if(any(!nm %in% names(filter)) || any(sapply(nm,function(thisNM)
                                                  !(is(filter[[thisNM]] ,"filter")||is(filter[[thisNM]] ,"filters"))
                                            )))
    {
          warning("'filter' must either be a filtersList,filterResultList, a single\n",
          "filter object or a named list of filter objects.",
          call.=FALSE)
      filter <- NULL
    }
  }
  filter
  
  
    
   
  
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
#' @export 
#' @rdname xyplot
panel.xyplot.flowset <- function(x,
                                 frames,
                                 filter = NULL,
                                 channel.x,
                                 channel.y,
								 overlay= NULL#a list of flowSet
                                ,stats = FALSE
						 ,...)
{
    nm <- as.character(x)
    if (length(nm) < 1) return()
   
    filter <- .processFilter(filter, nm)

    if(!is.list(stats)){
        stats <- lapply(seq_along(nm), function(x) stats)
        names(stats) <- nm
      
    }
    overlay <- .process_overlay_flowSet(overlay, nm)
    
    
    panel.xyplot.flowframe(frame=frames[[nm]], filter=filter[[nm]]
                              , overlay = overlay
                              , stats = stats[[nm]]
                               , ...)
}

#' extract the respective flowFrame from each flowSet based on the given sampleName
#' @param overlay a list of flowSet
#' @param nm sample name
.process_overlay_flowSet <- function(overlay, nm){
  if(!is.null(overlay)){
    if(!is.list(overlay))
      stop("overlay must be a list of flowSet!")
    overlay <- sapply(overlay, function(thisOverlay){
          thisClass <- class(thisOverlay)
          if(!(extends(thisClass, "flowSet")||extends(thisClass, "ncdfFlowList")))
            stop("overlay must be a list of 'flowSet'")
          thisOverlay[[nm]]
          
        })
  }
  overlay
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
#' @rdname xyplot
#' @aliases xyplot,formula,view-method
setMethod("xyplot",
          signature(x="formula", data="view"),
          function(x, data, ...) xyplot(x, Data(data), ...))

#' @rdname xyplot
#' @aliases xyplot,view,missing-method
setMethod("xyplot",
          signature(x="view", data="missing"),
          function(x, data, ...) xyplot(Data(x), ...))



## ==========================================================================
## Plot a gateView object. Essentially, this is calling the plot
## method on the parent data of the view, adding gates if appropriate.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @rdname xyplot
#' @aliases xyplot,formula,gateView-method
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



