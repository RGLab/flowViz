#' One-dimensional density/histogram plots for flow data
#' 
#' 
#' For \code{\link[flowCore:flowSet-class]{flowSets}} the idea is to
#' horizontally stack plots of density estimates for all frames in the
#' \code{flowSet} for one or several flow parameters. In the latter case, each
#' parameter will be plotted in a separate panel, i.e., we implicitly condition
#' on parameters.
#' 
#' 
#' Not all standard lattice arguments will have the intended effect, but many
#' should.  For a fuller description of possible arguments and their effects,
#' consult documentation on lattice (Trellis docs would also work for the
#' fundamentals).
#' 
#' @name densityplot
#' 
#' @param x A formula describing the structure of the plot and the variables to
#' be used in the display. The structure of the formula is \code{factor ~
#' parameter}, where \code{factor} can be any of the phenotypic factors in the
#' \code{phenoData} slot or an appropriate factor object and \code{parameter}
#' is a flow parameter. Panels for multiple parameters are drawn if the formula
#' structure is similar to \code{factor ~ parameter1 + parameter2}, and
#' \code{factor} can be missing, in which case the sample names are used as
#' y-variable. To facilitate programatic access, the formula can be of special
#' structure \code{factor ~ .}, in which case the optional \code{channel}
#' argument is considered for parameter selection. For the workflow methods,
#' \code{x} can also be one of the several workflow objects.
#' 
#' @param data A flow data object that serves as a source of data, either a
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} or
#' \code{\link[flowCore:flowSet-class]{flowSet}}
#' 
#' @param \dots More arguments, usually passed on to the underlying lattice
#' methods.  
#' \itemize{
#'
#'  \item channels A character vector of parameters that are supposed to be
#'                  plotted when the formula in \code{x} is of structure \code{factor ~ .}.
#' 
#'  \item xlab: Label for data x axis, with suitable defaults taken from the
#'  formula
#' 
#'  \item prepanel: The prepanel function.  See \code{\link[lattice]{xyplot}}
#' 
#'  \item panel: the panel function.  See \code{\link[lattice]{xyplot}}
#' 
#'  \item axis: axis function passed to lattice, default is \code{axis.grid}
#' 
#'  \item ... : other arguments passed to panel.densityplot.flowset.stack or
#'              panel.histogram.flowframe
#' 
#' }
#' @examples
#' 
#' library(flowStats)
#' data(GvHD)
#' GvHD <- GvHD[pData(GvHD)$Patient %in% 6:7]
#' 
#' densityplot(~ `FSC-H`, GvHD)
#' 
#' densityplot(~ `FSC-H` + `SSC-H`, GvHD)
#' 
#' densityplot(~ ., GvHD[1:3])
#' 
#' ## include a filter
#' densityplot(~ `FSC-H`, GvHD, filter=curv1Filter("FSC-H"))
#' 
#' #display the gate by its boundaries with statistics 
#' densityplot(~ `FSC-H`, GvHD[1:2], filter=curv1Filter("FSC-H"),fitGate=FALSE,stats=TRUE)
#' 
#' ## plot a single flowFrame
#' densityplot(~ `SSC-H`, GvHD[[1]], margin=FALSE)
#' 
#' ## plot histogram
#' histogram(~ `SSC-H`, GvHD[[1]]) #default type is 'density'
#' #change the type to 'count' and adjust breaks
#' histogram(~ `SSC-H`, GvHD[[1]], margin=FALSE, type = "count", breaks = 50)
#' 
#' @aliases densityplot,formula,flowSet-method
#' @export 
#' @rdname densityplot
setMethod("densityplot",
    signature(x = "formula", data = "flowSet"),
    function(x, data, ...){
      #construct lattice object
      thisTrellisObj <- .densityplot.adapor(x, data, ...) 
      
      #pass frames slot
      thisData <- thisTrellisObj[["panel.args.common"]][["frames"]]
      thisTrellisObj[["panel.args.common"]][["frames"]] <- thisData@frames
      thisTrellisObj
    })

## Dedicated prepanel function to set up dimensions
#' @export 
#' @rdname densityplot
prepanel.densityplot.flowset.stack <- 
    function(x, y, frames, 
             overlap=0.3, subscripts, ...,
             which.channel)
{
    channel.name <- unique(which.channel[subscripts])
    stopifnot(length(channel.name) == 1)
    if(extends(class(frames),"flowSet"))
      ranglist <- eapply(frames@frames, function(fr)range(fr)[, channel.name])
    else
      ranglist <- lapply(sampleNames(frames), function(sn)range(frames[[sn, use.exprs = FALSE]])[, channel.name])
    xl <- range(ranglist, finite=TRUE)
    list(xlim=xl + c(-1,1)*0.07*diff(xl))   
}


## add a barchart indicating the margin events
mbar <- function(dat, p, r, i, col, m)
{
    rl <- r + c(-1,1)*min(100, 0.06*diff(r))
    off <- min(5, 0.005*diff(r))
    lx <- length(dat)
    des <- -20
    fac <- 0.7
    for(j in seq_along(p)){
        sp <- sum(p[[j]], na.rm=TRUE)
        if(sp > lx*m)
        {
            panel.segments(rl[j], i, rl[j], i+fac, col="gray")
            panel.lines(rep(rl[j],2), c(i, i+(sp/lx*fac)),
                        col=desat(col, des), lwd=4)
            panel.segments(rl[j]-off, i, rl[j]+off, i, col="gray")
            panel.segments(rl[j]-off, i+fac, rl[j]+off, i+fac, col="gray")

        }
    }
}


    
## Dedicated panel function to do the plotting and add gate boundaries
#' @param overlay see help(xyplot).
#' 
#' @param overlap The amount of overlap between stacked density plots. This
#' argument is ignored for the \code{flowFrame} method.
#' 
#' @param filter A \code{\link[flowCore:filter-class]{filter}},
#' \code{\link[flowCore:filterResult-class]{filterResult}} or
#' \code{\link[flowCore:filterResult-class]{filterResultList}} object or a list
#' of such objects of the same length as the \code{flowSet}. If applicable, the
#' gate region will be superiposed on the density curves using color shading.
#' The software will figure out whether the \code{filter} needs to be evaluated
#' in order to be plotted (in which case providing a \code{filterResult} can
#' speed things up considerably).
#' 
#' @param groups Use identical colors for grouping. The value of the argument
#' is expected to be a phenotypic variable in the \code{flowSet}, or a factor.
#' 
#' 
#' @param subscripts,which.channel,channel.name,y Internal indices necessary to
#' map panels to parameters.
#' 
#' @param frames An environment containing frame-specific data.
#' @param channel The name of the currently plotted flow parameter.
#' @param darg These arguments are passed unchanged to the corresponding
#' methods in lattice, and are listed here only because they provide different
#' defaults.  See documentation for the original methods for details.
#' \code{darg} gets passed on to \code{\link[stats]{density}}.
#' @param col,fill,lty,lwd,alpha Graphical parameters. These mostly exist for
#' conveniance and much more control is available throught the
#' \code{lattice}-like \code{par.setting} and \code{flowViz.par.set}
#' customization. The relevant parameter category for density plots is
#' \code{gate.density} with available parameters \code{col}, \code{fill},
#' \code{lwd}, \code{alpha} and \code{lty}. See
#' \code{\link[flowViz:flowViz.par.get]{flowViz.par.set}} for details.
#' @param refline Logical. Add one ore more vertical reference lines to the
#' plot.  This argument is directly passed to
#' \code{\link[lattice:panel.functions]{panel.abline}}.
#' @param margin Either Logical value 'FALSE' or Numeric valuein \code{[0,1]}.
#' When 'FALSE', it doesn't do anything to the margin events.  When Numeric
#' value, it indicates margin events by horizontal bars. The value of
#' \code{margin} is interpreted as the proportion of events on the margin over
#' which the bars are added. E.g., a value of \code{0,5} means to indicate
#' margin events if there are more than \code{0.5} times the total number of
#' events. \code{1} means to ignore margin events completetly. For \code{0}
#' bars are added even if there is only a single margin event.
#' @param stats,pos,digits,abs Arguments to control statistics that is
#' associated with \code{\link[flowCore:filter-class]{filter}} to be plotted.
#' see \link[flowViz]{xyplot} for details.
#' @param fitGate A \code{logical} scalar indicating whether to display the
#' gate as fitted 1d density gate region or simply display the gate boundaries
#' using vertical lines. The latter would be helpful to display the gate when
#' the gated density region is too small to see.
#' @param checkName A \code{logical} scalar indicating whether to validity
#' check the channel name. Default is TRUE, which consider '(' as invalid
#' character in channel names
#' @param gp A list of graphical parameters that are passed down to the low
#' level panel functions. This is for internal use only. The public user
#' interface to set graphical parameters is either \code{par.settings} for
#' customization of a single call or \code{flowViz.par.set} for customization
#' of session-wide defaults.
#' @param plotType either 'densityplot' or 'histogram'
#' @param hist.type see 'type' argument in 'help(panel.histogram)'
#' @param breaks see 'help(hist)' 
#' @importFrom stats approxfun cor density median ppoints qnorm quantile qunif var
#' @export
#' @rdname densityplot
panel.densityplot.flowset.stack <-
    function(x, y, darg=list(n=50, na.rm=TRUE), frames, channel,
             overlap = 0.3, channel.name, filter=NULL,
             fill=superpose.polygon$col,
             lty=superpose.polygon$lty,
             lwd=superpose.polygon$lwd,
             alpha=superpose.polygon$alpha,
             col=superpose.polygon$border,
             groups=NULL, refline=NULL,
             margin=0.005
			 ,stats=FALSE
			 ,pos=0.5
			 ,digits=2
			 ,abs=FALSE
			 ,fitGate=TRUE
             ,checkName = TRUE
             , plotType = "densityplot"
             , hist.type = "density"
             , breaks = "Sturges"
             ,gp, ...)
{
    which.channel <- tail(which.packet(), 1)
    
    lc <- length(channel)
    channel <- channel[[which.channel]]
    channel.name <- channel.name[which.channel]
    ycode <- as.numeric(y)
    validName <- !length(grep("\\(", channel.name))
    if(checkName)
      validName <- !(length(grep("\\(", channel.name)) ||
            length(grep("\\(", channel.name)))
    else
      validName <- TRUE
    
    if(!validName)
      warning("Gate will not be plotted because channel names contain '(' character! Try to set checkName to FALSE to skip this check.")

    if (any(duplicated(ycode)))
        warning("Some combinations seem to have multiple samples.  \n  ",
                "Only one will be used.")
    nnm <- as.character(x)
    filter <- .processFilter(filter, nnm, lc = lc, which.channel = which.channel)  
    ny <- nlevels(y)
	superpose.polygon <- flowViz.par.get("superpose.polygon")
    border <- rep(col, length = ny)
    col <- rep(fill, length = ny)
    if(!is.null(groups))
        col <- col[groups]
    lty <- rep(lty, length = ny)
    lwd <- rep(lwd, length = ny)
    alpha <- rep(alpha, length = ny)
    x <- as.character(x)
    height <- (1 + overlap)
    parm <- gsub("`", "", as.character(channel))
	ptList<-plotType("gdensity", parm)
    for (i in rev(seq_len(ny))){
        if (i %in% ycode)
        {
            nm <- x[match(i, ycode)]
            xx <- evalInFlowFrame(channel, frames[[nm]])
            r <- unlist(range(frames[[nm]])[, channel.name])
            
            if(!is.logical(margin)||isTRUE(margin))#when margin is logical FALSE we skip marginal events filtering
            {
              margin <- min(1, max(0, margin))
              pl <- xx<=r[1]
              pr <- xx>=r[2]
              xxt <- xx[!(pl | pr)]
              ## we indicate piled up data by vertical lines (if > 1%) unless
              ## margin=FALSE
              if(margin<1)
                  mbar(xx, list(pl, pr), r, i, col[i], margin)
              }else
                xxt <- xx
            ## we need a smaller bandwidth than the default and keep it constant
            if(length(xxt)){
                
                if(plotType == "densityplot"){
                  #densityplot
                  if(!("adjust" %in% names(darg)))
                    darg[["adjust"]] <- 2
                  h <- do.call(density, c(list(x=xxt), darg))
                  x.val <- h$x
                  y.val <- h$y
                  n <- length(x.val)
                }else{
                  fitGate <- FALSE
                  #histogram
                  h <- lattice:::hist.constructor(xxt, breaks = breaks, ...)
                  y.val <- switch(hist.type, count = h$counts, percent = 100 * h$counts/length(x), density = h$density)
                  x.val <- h$breaks
                  n <- length(x.val)
                  if (length(y.val) != n - 1) 
                    warning("problem with 'hist' computations")
                }

                max.d <- max(y.val)
                xl <- x.val[c(1, 1:n, n)]
                yl <- i + height * c(0, y.val, 0) / max.d
                
                if(plotType == "densityplot")
                  panel.polygon(x=xl,y=yl, col=col[i], border=NA, alpha=alpha[i])
                else
                  panel.rect(x = xl[-n], y = 0, height = yl, width = diff(xl), 
                      col = col[i], alpha = alpha[i], border = border[i], just = c("left", "bottom"), identifier = "histogram")
                
                ## add the filterResult if possible, we get them from the output of
                ## glpolygon (with plot=FALSE)
				
                if(!is.null(filter[[nm]]) && validName){
                  curFilter <- filter[[nm]]
                  
                  if(is.list(stats))
                    curStats <- stats[[nm]]
                  else
                    curStats <- stats
#                browser()
                if(is(curFilter,"filters"))
                {
                  mapply(curFilter,curStats,FUN=function(thisFilter,thisStats){
                         names <- .getStats(thisFilter,thisStats, frames[[nm]], digits, ...)
                        

                        #this plot routine is only for 2-d scatter plot
                        #thus set plot as FALSE,just use it to get bounds
                        #add plot stats/names
                        bounds <- glpolygon(thisFilter, frames[[nm]],
                            channels=parm
                            ,ptList=ptList
                            ,verbose=FALSE
                            , plot=FALSE 
                            ,names=names
                            ,xlim=r
                            ,ylim=c(i,i+1)
                            ,strict=FALSE)
                        oo <- options(warn=-1)
                        on.exit(options(oo))
                        if(any(!is.na(bounds))){
                          ## iterate over gate regions
                          for(j in seq_along(bounds)){
                            if(any(!is.na(bounds[[j]]))){
                              tb <- bounds[[j]]
                              if(fitGate)
                              {
                                if(ncol(tb) == 1 && colnames(tb) == parm){
                                  sel <- xl >= min(tb) & xl <= max(tb)
                                  if(any(sel)){
                                    afun <- approxfun(xl, yl)
                                    xr <- c(min(tb), seq(min(tb), max(tb), len=100),
                                            max(tb))
                                    yr <- c(i, afun(xr[-c(1, length(xr))]), i)
                                    gpd<-gp$gate.density
                                    panel.polygon(xr, yr
                                                  , border=gpd$col, col=gpd$fill,
                                                  alpha=gpd$alpha, lwd=gpd$lwd,
                                                  lty=gpd$lty
                                    )
                                  }
                                }	
                              }else
                              {
                                gpg<-gp$gate
                                panel.lines(x=c(tb[1],tb[1]),y=c(i,i+height)
                                            ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                                panel.lines(x=c(tb[2],tb[2]),y=c(i,i+height)
                                            ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                              }
                            }
                          }
                        }
                        
                        
                        
                        options(oo)
                      })
                }else
                {
					curFilter<-filter[[nm]]
#					browser()
                    if(is.list(stats))
                      thisStats <- stats[[nm]]
                    else
                      thisStats <- stats
                    names <- .getStats(curFilter,thisStats, frames[[nm]], digits, ...)
					
#					browser()
					#this plot routine is only for 2-d scatter plot
					#thus set plot as FALSE,just use it to get bounds
					#add plot stats/names
                    bounds <- glpolygon(curFilter, frames[[nm]],
                                        channels=parm
										,ptList=ptList
                                        ,verbose=FALSE
										, plot=FALSE 
										,names=names
										,xlim=r
										,ylim=c(i,i+1)
										,strict=FALSE)
                    oo <- options(warn=-1)
                    on.exit(options(oo))
                    if(any(!is.na(bounds))){
                        ## iterate over gate regions
                        for(j in seq_along(bounds)){
                          if(any(!is.na(bounds[[j]]))){
                            tb <- bounds[[j]]
                            if(fitGate)
                            {
                              if(ncol(tb) == 1 && colnames(tb) == parm){
                                sel <- xl >= min(tb) & xl <= max(tb)
                                if(any(sel)){
                                  afun <- approxfun(xl, yl)
                                  xr <- c(min(tb), seq(min(tb), max(tb), len=100),
                                          max(tb))
                                  yr <- c(i, afun(xr[-c(1, length(xr))]), i)
                                  gpd<-gp$gate.density
                                  panel.polygon(xr, yr
                                                , border=gpd$col, col=gpd$fill,
                                                alpha=gpd$alpha, lwd=gpd$lwd,
                                                lty=gpd$lty
                                  )
                                }
                              }	
                            }else
                            {
                              gpg<-gp$gate
                              panel.lines(x=c(tb[1],tb[1]),y=c(i,i+height)
                                          ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                              panel.lines(x=c(tb[2],tb[2]),y=c(i,i+height)
                                          ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                            }
                          }
                        }
                    }
					
					
					
                    options(oo)
                }
               }
#				browser()
                panel.lines(x=xl,y=yl, col=border[i], lty=lty[i],lwd=lwd[i])

            }else{
                panel.lines(rl, rep(i,2), col="black")
            }
			
        }
    }
    if(!is.null(refline))
        panel.abline(v=refline)
}


prepanel.densityplot.flowset <- 
    function(x, frames, channel.x.name, xlim ,margin = TRUE, hist.type = "density", breaks = "sturges", ...)
{
  if (length(nm <- as.character(x)) > 1)
    stop("must have only one flow frame per panel")
  
  if (length(nm) == 1){
    ranges <- range(frames[[nm]])
    if(missing(xlim)){
      #use default calculation
      xlim <- if(!length(grep("`", channel.x.name, fixed=TRUE))){
            tmp <- ranges[, channel.x.name]
            xd <- diff(tmp)/15
            tmp + c(-1,1)*xd
          }else NULL
    }
    
    if(hist.type %in% c("count", "percent")){
      #run hist to estimate the max count
      data <- as.vector(exprs(frames[[nm, channel.x.name]]))
      
      if(margin){
        r <- ranges[, channel.x.name]
        pl <- data<=r[1]
        pr <- data>=r[2]
        data <- data[!(pl | pr)]  
      }
      h <- lattice:::hist.constructor(data, breaks = breaks, ...)
      y.val <- switch(hist.type, count = h$counts, percent = 100 * h$counts/length(data), density = h$density)
      ylim = c(0, max(y.val))
    }else
      ylim = c(0,1)#when hist.type = density, y will be normalized to c(0,1), thus no need to compute ylim  
    return(list(xlim=xlim, ylim = ylim))
  }else
    return(list())
}

panel.densityplot.flowset <- function(x,
    frames,
    filter = NULL,
    channel.x,
    channel.y
    ,stats = FALSE
    , overlay= NULL#a list of flowSet
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
  
  panel.densityplot.flowFrame(frame=frames[[nm]], filter=filter[[nm]], stats = stats[[nm]], overlay = overlay, ...)
}
panel.densityplot.flowFrame <-
    function(frame,
        darg=list(n=50, na.rm=TRUE),
        filter=NULL,
        margin=0.005,
        outline=FALSE,
        channel.x.name,
        channel.y.name,
        fill= "#619CFF",
        lty= 1,
        lwd= 1,
        alpha= 1,
        col= "transparent",
        gp
        ,stats = FALSE
        ,pos=0.5
        ,digits=2
        ,abs=FALSE
        ,fitGate=TRUE
        ,refline=NULL
        ,checkName = TRUE
        , overlay = NULL
        , overlay.symbol = NULL
        , plotType = "densityplot"
        , hist.type = "density"
        , breaks = "Sturges"
        ,...
        )
{
          
  
          
          validName <- !length(grep("\\(", channel.x.name))
          if(checkName)
            validName <- !(length(grep("\\(", channel.x.name)) ||
                  length(grep("\\(", channel.y.name)))
          else
            validName <- TRUE
          
          if(!validName)
            warning("Gate will not be plotted because channel names contain '(' character! Try to set checkName to FALSE to skip this check.")
          border <- col
          parm <- channel.x.name
          ptList<-plotType("gdensity", parm)
         
          xx <- exprs(frame)[,channel.x.name]
          r <- unlist(range(frame)[, channel.x.name])
          
          ## we ignore data that has piled up on the margins
          rl <- r + c(-1,1)*min(100, 0.06*diff(r))
          if(!is.logical(margin)||isTRUE(margin))#when margin is logical FALSE we skip marginal events filtering
          {
            margin <- min(1, max(0, margin))
            pl <- xx<=r[1]
            pr <- xx>=r[2]
            xxt <- xx[!(pl | pr)]
            ## we indicate piled up data by vertical lines (if > 1%) unless
            ## margin=FALSE
            if(margin<1)
              mbar(xx, list(pl, pr), r, 1, col, margin)
          }else
            xxt <- xx
          if(length(xxt)){
            if(plotType == "densityplot"){
              #densityplot
              if(!("adjust" %in% names(darg)))
                darg[["adjust"]] <- 2
              
              h <- do.call(density, c(list(x=xxt), darg))
              x.val <- h$x
              y.val <- h$y
              n <- length(x.val)
            }else{
              fitGate <- FALSE
              #histogram
              h <- lattice:::hist.constructor(xxt, breaks = breaks, ...)
              y.val <- switch(hist.type, count = h$counts, percent = 100 * h$counts/length(xxt), density = h$density)
              x.val <- h$breaks
              n <- length(x.val)
              if (length(y.val) != n - 1) 
                warning("problem with 'hist' computations")
            }
            
            max.d <- max(y.val)
            xl <- x.val[c(1, 1:n, n)]
            yl <- c(0, y.val, 0)
            if(hist.type == "density")
              yl <- yl / max.d
            
            if(plotType == "densityplot")
              panel.polygon(x=xl,y=yl, col=fill, border=NA, alpha=alpha)
            else
              panel.rect(x = xl[-n], y = 0, height = yl, width = diff(xl), 
                  col = fill, alpha = alpha, border = col, just = c("left", "bottom"), identifier = "histogram")
                        
#             browser()           
            if(!is.null(overlay))
            {
              overlayNames <- names(overlay)
              for(overLayName in overlayNames){
                thisOverlay <- overlay[[overLayName]]
                if(!is.null(thisOverlay)){
                  thisDat <- exprs(thisOverlay)
                  overlay.x <- thisDat[, channel.x.name]
                  
                  
                  this.overlay.symbol <- gp$overlay.symbol#default setting
                  user.overlay.symbol <- overlay.symbol[[overLayName]]#customized settings
                  #update the default with customized settings
                  if(!is.null(user.overlay.symbol))
                    this.overlay.symbol <- lattice:::updateList(this.overlay.symbol, overlay.symbol[[overLayName]])
#                browser()                
                  h <- do.call(density, c(list(x=overlay.x), darg))
                  n <- length(h$x)
                  
                  xl_overlay <- h$x[c(1, 1:n, n)]
                  yl_overlay <- c(0, h$y, 0) / max(h$y) 
                  
                  
                  panel.polygon(x=xl_overlay,y=yl_overlay 
                      ,  border=NA
                      , col = this.overlay.symbol[["fill"]]
                      , alpha = this.overlay.symbol[["alpha"]]
                  
                  )  
                }
                
              }
              
            }                        
            ## add the filterResult if possible, we get them from the output of
            ## glpolygon (with plot=FALSE)
            if(!is.null(filter) && validName){    
              if(is(filter,"filters"))
              {
                fitGate <- FALSE #it is hard to plot multiple fitted gates
                
                mapply(filter,stats,FUN=function(curFilter,curStats){
                      
#                      browser()
                      names <- .getStats(curFilter,curStats, frame, digits, ...)
                      
                      
                      #this plot routine is only for 2-d scatter plot
                      #thus set plot as FALSE,just use it to get bounds
                      #add plot stats/names
                      bounds <- glpolygon(curFilter, frame,
                          channels=parm
                          ,ptList=ptList
                          ,verbose=FALSE
                          , plot=FALSE 
                          ,names=names
                          ,xlim=r
                          ,ylim=c(0,1)
                          ,strict=FALSE)
                      oo <- options(warn=-1)
                      on.exit(options(oo))
                      if(any(!is.na(bounds))){
                        ## iterate over gate regions
                        for(j in seq_along(bounds)){
                          if(any(!is.na(bounds[[j]]))){
                            tb <- bounds[[j]]
                            gpg<-gp$gate
                            panel.lines(x=c(tb[1],tb[1]),y=c(0,1)
                                        ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                            panel.lines(x=c(tb[2],tb[2]),y=c(0,1)
                                        ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                          }
                        }
                      }
                      
                      options(oo)
                    })
              }else
              {
                  names <- .getStats(filter,stats, frame, digits, ...)
                  
        
                  #this plot routine is only for 2-d scatter plot
                  #thus set plot as FALSE,just use it to get bounds
                  #add plot stats/names
                  bounds <- glpolygon(filter, frame,
                      channels=parm
                      ,ptList=ptList
                      ,verbose=FALSE
                      , plot=FALSE 
                      ,names=names
                      ,xlim=r
                      ,ylim=c(0,1)
                      ,strict=FALSE)
                  oo <- options(warn=-1)
                  on.exit(options(oo))
                  if(any(!is.na(bounds))){
                    ## iterate over gate regions
                    for(j in seq_along(bounds)){
                      if(any(!is.na(bounds[[j]]))){
                        tb <- bounds[[j]]
                        if(fitGate)
                        {
                          tb_min <- min(tb)
                          tb_max <- max(tb)
                          
                          tb_min <- max(min(xl),tb_min)
                          tb_max <- min(max(xl),tb_max)
                          if(ncol(tb) == 1 && colnames(tb) == parm){
                            sel <- xl >= tb_min & xl <= tb_max
                            if(any(sel)){
                              afun <- approxfun(xl, yl)
                              xr <- c(tb_min, seq(tb_min, tb_max, len=100), tb_max)
                              #                        browser()
                              yr <- c(0,afun(xr[-c(1, length(xr))]),0)
                              gpd<-gp$gate.density
                              panel.polygon(xr, yr
                                            , border=gpd$col
                                            , col=gpd$fill,
                                            alpha=gpd$alpha, lwd=gpd$lwd,
                                            lty=gpd$lty
                              )
                            }
                          }	
                        }else
                        {
                          gpg<-gp$gate
                          panel.lines(x=c(tb[1],tb[1]),y=c(0,1)
                                      ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                          panel.lines(x=c(tb[2],tb[2]),y=c(0,1)
                                      ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                        }
                      }
                    }
                  }
                  
                  
                  
                  options(oo)
              }
            }

            panel.lines(x=xl,y=yl, col=border, lty=lty,lwd=lwd)
    
          }else{
            panel.lines(rl, rep(0,2), col="black")
          }
      if(!is.null(refline))
        panel.abline(v=refline)
}

analyzeDensityFormula <- function(x, dot.names)
{
    ans <- list()
    if (length(x) == 2) {
        ans$left <- FALSE
        x <- x[[2]]
    }
    else if (length(x) == 3) {
        ans$left <- TRUE
        ans$left.symbol <- x[[2]]
        x <- x[[3]]
    }
    else stop("unexpected formula structure")
    if (length(x) == 1) {
        ans$conditioned <- FALSE
        ans$multiple.right <- FALSE
        ans$right.symbol <- x
    }
    else if (length(x) == 3) {
        switch(as.character(x[[1]]),
               "+" = {
                   ans$conditioned <- FALSE
                   ans$multiple.right <- TRUE
                   ans$right.symbol <- x
               },
               "|" = {
                   ans$conditioned <- TRUE
                   ans$cond.symbol <- x[[3]]
                   ans$right.symbol <- x[[2]]
                   ans$multiple.right <-
                       (length(x[[2]]) > 1 &&
                        as.character(x[[2]][[1]]) == "+")
               })
    }
    else stop("unexpected formula structure")
    if (ans$multiple.right) {
        x <- ans$right.symbol
        right.comps <- list()
        while ((length(x) > 1) && (as.character(x[[1]]) == "+")) {
            right.comps <- c(list(x[[3]]), right.comps)
            x <- x[[2]]
        }
        ans$right.comps <- c(list(x), right.comps)
    }
    else if (length(ans$right.symbol) == 1 && as.character(ans$right.symbol) == ".") {
        ans$multiple.right <- TRUE
        ans$right.comps <- lapply(dot.names, as.symbol)
    }
    else ans$right.comps <- list(ans$right.symbol)
    ans
}




#' dispatch to different trellis object contructing function based on stack argument
#' @param stack \code{logical} indicating whether to stack `name` on y axis
.densityplot.adapor <- function(x, data, stack = TRUE, plotType = "densityplot", ...){
  
      plotType <- match.arg(plotType, c("densityplot", "histogram"))
      if(stack)
        thisObj <- .densityplot.flowSet.stack(x, data, plotType = plotType, ...)
      else
      {
    
        #add dummy y term in order to use .xyplot.flowSet to construct lattice object
        if(length(x[[2]]) == 1)
          xTerm <- x[[2]]
        else
          xTerm <- x[[2]][[2]]
        thisFormula <- eval(substitute(thisX~y, list(thisX = xTerm)))
        thisFormula[[3]] <- x[[2]]
        thisObj <- .xyplot.flowSet(thisFormula, data
            , panel = panel.densityplot.flowset
            , prepanel = prepanel.densityplot.flowset
            , ylab = ""
            , plotType = plotType
            ,...)
        #append channel.name to work ncdfFlow::densityplot
        thisObj$panel.args.common$channel.name <- thisObj$panel.args.common$channel.x.name
      }
      thisObj
}

.densityplot.flowSet.stack <- function(x, data, channels, xlab,
                   as.table = TRUE, overlap = 0.3,
                   prepanel = prepanel.densityplot.flowset.stack,
                   panel = panel.densityplot.flowset.stack,
                   filter=NULL, scales=list(y=list(draw=F)),
                   groups, axis= axis.grid
#                    , marker.only = FALSE
                    , ...)
      {
          ocall <- sys.call(sys.parent())
          ccall <- match.call(expand.dots = TRUE)
          
          
          if(! "name" %in% names(pData(data)))
              pData(data)$name <- sampleNames(data)
          pd <- pData(phenoData(data))
          uniq.name <- createUniqueColumnName(pd)
          ## ugly hack to suppress warnings about coercion introducing
          ## NAs (needs to be `undone' inside prepanel and panel
          ## functions):
          pd[[uniq.name]] <-
              factor(sampleNames(data),
                     levels=unique(sampleNames(data))) 
          if (missing(channels))
              channels <-
                  setdiff(colnames(data),
                          flowCore:::findTimeChannel(data))

          formula.struct <- analyzeDensityFormula(x, dot.names = channels)
          ## we want to add a column to pd for each channel, repeating
          ## pd as necessary.  We might want to skip this if there is
          ## only one channel, but for now we'll use it for
          ## conditioning even then.

          channel.name <-
              sapply(formula.struct$right.comps, expr2char)
          
          frm <- data[[1, use.exprs = FALSE]]
          xObjs <- sapply(channel.name,  getChannelMarker, frm = frm, simplify = FALSE)
          channel.name <- sapply(xObjs, "[[", "name")
          formula.struct$right.comps <- lapply(channel.name, function(i)as.symbol(i))
          
          
          pd <- rep(list(pd), length(channel.name))
#          browser()
#          if(marker.only){
#            names(pd) <-  sapply(xObjs, function(xObj)as.character(ifelse(is.na(xObj[,"desc"]), xObj[,"name"], xObj[,"desc"])))
#               
#          }else{
#            names(pd) <-  sapply(xObjs, function(xObj)sub("NA","",paste(unlist(xObj),collapse=" "))) 
#          }
          names(pd) <- channel.name
          pd <- do.call(lattice::make.groups, pd)
          ## FIXME: this won't work if pd already has a column named
          ## 'which'.  Should deal with that case somehow.

          ## Next task is to manipulate the formula.  The details of
          ## the transformation depends on whether there is a
          ## conditioning variable alread.
          ## y ~ channel ==> y ~ sample | which
          ## y ~ channel | var ==> y ~ sample | which + var

         
          new.x <- d1 ~ d2 | d3
          new.x[[2]] <- ## d1
              if (formula.struct$left) formula.struct$left.symbol
              else as.name("name")
          new.x[[3]][[2]] <- ## d2
              as.name(uniq.name)
          new.x[[3]][[3]] <- ## d3
              if (formula.struct$conditioned) {
                  ans <- (~.+.)[[2]]
                  ans[[3]] <- as.name("which")
                  ans[[2]] <- formula.struct$cond.symbol
                  ## probably not the ideal order, but I don't see how
                  ## to easily ake 'which' the first conditioning
                  ## variable (in case there is more than one
                  ## conditioning variable to begin with)
                  ans
              }
              else as.name("which")
          if (missing(xlab))
              xlab <- ""
            
          gp <- list(...)[["par.settings"]]
          gpar <- flowViz.par.get()
          if(!is.null(gp))
              gpar <- lattice:::updateList(gpar, gp)
          if(!missing(groups))
              ccall$groups <- as.factor(eval(substitute(groups), pData(data)))
          ccall$gp <- gpar
          ccall$x <- new.x
          ccall$data <- pd
          ccall$prepanel <- prepanel
          ccall$panel <- panel
          ccall$par.settings <- gpar
          ccall$channel <- formula.struct$right.comps ## channel
          ## That is super ugly!!! How do we get to the channel name
          ## from the formula???
          ccall$channel.name <- gsub("^.*\\(`|`\\).*$", "", channel.name)
          ccall$frames <- data
          ccall$as.table <- as.table
          overlap <- max(-0.5, overlap)
          ccall$overlap <- overlap
          ccall$axis <- axis
          ccall$xlab <- xlab
          ccall$horizontal <- TRUE
          ccall$subscripts <- TRUE
          ccall$default.scales <- list(x = list(relation = "free"))
          ccall$which.channel <-
              gsub("^.*\\(`|`\\).*$", "", as.character(pd$which))
          ypad <- lattice.getOption("axis.padding")$numeric * (length(data) + overlap)
          ccall$lattice.options <-
              list(axis.padding = list(factor = c(ypad, ypad + 1 + overlap)))
          ccall[[1]] <- quote(lattice::bwplot)
          ans <- eval.parent(ccall)
          ans$call <- ocall
          ans
      }



#' @export
#' @rdname  densityplot
setMethod("densityplot",
          signature(x="formula", data="flowFrame"),
          function(x, data, overlay = NULL, ...){
            sn <- identifier(data)
            data <- as(data, "flowSet")
            sampleNames(data) <- sn
            
            overlay <- .process_flowFrame_overlay(overlay, sn)
                        
            densityplot(x, data, overlap = 0, scales = list(y=list(draw=FALSE)), overlay = overlay, ...)              
          })


#' @export
#' @rdname  densityplot
setMethod("densityplot",
          signature(x="formula", data="view"),
          function(x, data, ...) densityplot(x, Data(data), ...))





