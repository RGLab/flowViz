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



## xyplot for flowFrames with a formula.  We'll make this very simple;
## the drawback being that the expression matrix will be copied,
## the upshot being that all the fancy xyplot formula stuff will be valid.
#setMethod("xyplot",
#          signature=signature(x="formula",
#                              data="flowFrame"),
#          definition=function(x,
#                              data,
#                              smooth=TRUE,
#                              prepanel=prepanel.xyplot.flowframe,
#                              panel=panel.xyplot.flowframe,
#							  overlay=NULL,#a flowframe
#                              ...)
#      {
#          ## par.settings will not be passed on to the panel functions, so
#          ## we have to fetch it from ... and stick the gate relevant stuff
#          ## back it in there manually
#          gp <- list(...)[["par.settings"]]
#          
#          ## deparse the formula structure
#          ## FIXME: shouldn't all.vars be helpful here?!?
#          channel.y <- x[[2]]
#          channel.x <- x[[3]]
#          if (length(channel.x) == 3)
#              channel.x <- channel.x[[2]]
#          channel.x.name <- expr2char(channel.x)
#          channel.y.name <- expr2char(channel.y)
##		  browser()
#		if(!is.null(overlay))
#		{
#			overlay.x <- flowViz:::evalInFlowFrame(channel.x, overlay)
#			overlay.y <- flowViz:::evalInFlowFrame(channel.y, overlay)
#		}else
#		{
#			overlay.x <- NULL
#			overlay.y <-NULL
#		}
#        data <- data[,c(channel.x.name,channel.y.name)]
#        
#          ## call data.frame xyplot method with our panel function
#          xyplot(x, data=as.data.frame(exprs(data)), smooth=smooth,
#                 prepanel=prepanel, panel=panel, frame=data,
#                 channel.x.name=channel.x.name,
#                 channel.y.name=channel.y.name,
#				 overlay.x=overlay.x,overlay.y=overlay.y,
#                 gp=gp, ...)
#      })

setMethod("xyplot",
    signature=signature(x="formula",
        data="flowFrame"),
    definition=function(x,
                        data
                        , filter = NULL
                        , overlay= NULL#a flowset
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
  
  if(!is.null(overlay)){
    overlay <- list(overlay)
    names(overlay) <- sn
  }
  
  if(!is.null(stats)){
    stats <- list(stats)
    names(stats) <- sn
  }
  
    
  xyplot(x, data, filter = filter
          , overlay = overlay, stats = stats
          , defaultCond = defaultCond , ...)  
})
## Prepanel function to set up dimensions. We want to use the instrument measurement
## range instead of the absolute range of the data. We also record the data ranges
## in the internal state environment for further use.
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
#overlay is a list(x=,y=), which is the extra points need to be plotted on top of the x,y
#' @param checkName \code{logical} indicating whether to skip checking the bracket '(' in channel name
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
									,overlay.x=NULL
									,overlay.y=NULL
                                    ,checkName = TRUE
						   			,...)
{
  
	x <-exprs(frame)[,channel.x.name]
    y <-exprs(frame)[,channel.y.name]
    

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
	l <- length(x)
		
		
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
            if(!is.null(overlay.x)){
              argcolramp <- .colRmpPlt(alpha = gp$overlay.symbol$bg.alpha)
            }
              
			if(xbins > 0)
			{
				#using hexbin package to do the hexagon plot	
				bin <- hexbin(x,y,xbins=xbins)
    
				grid.hexagons(bin,colramp = argcolramp, trans=binTrans, border =0)		
				
				
			}else
			{
				
				if(gp$density)
					col <- densCols(x, y, colramp=argcolramp)
            
				panel.xyplot(x, y, col=col, cex=cex, pch=pch, alpha=alpha, ...)
	
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
	if(!is.null(overlay.x)&&!is.null(overlay.y))
	{
        lpoints(overlay.x, overlay.y
                  , col = gp$overlay.symbol$fill
                  , alpha = gp$overlay.symbol$alpha
                  , cex= gp$overlay.symbol$cex 
                  , pch = gp$overlay.symbol$pch
                  )
		
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
setMethod("xyplot",
          signature=signature(x="formula",
                              data="flowSet"),
          definition= function(x,data, ...){
            thisTrellisObj <- .xyplot.flowSet(x, data, ...)
            thisData <- thisTrellisObj[["panel.args.common"]][["frames"]]
            thisTrellisObj[["panel.args.common"]][["frames"]] <- thisData@frames
            thisTrellisObj
          })
      
## flowViz:::.xyplot.flowSet now passes data instead of data@frames 
## within flowViz::xyplot method that changes it back to data@frames
## however ncdfFlow::xyplot keeps data as it is
.xyplot.flowSet <- function(x,
                              data,
                              smooth=TRUE,
                              filter=NULL,
                              as.table=TRUE,
                              prepanel=prepanel.xyplot.flowset,
                              panel=panel.xyplot.flowset,
                              xlab=channel.x.name,
                              ylab=channel.y.name,
                              par.settings=NULL
                              , axis= axis.grid
                              ,defaultCond = "name" #to override the default conditional variable 'name'
                                            #mainly used for plotting single flowFrame
                              , between = list(x=0.2,y=0.2)                              
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
          channel.x <- as.expression(channel.x)
          channel.y <- as.expression(channel.y)
          ## use densityplot method with dedicated panel and prepanel
          ## functions to do the actual plotting
          densityplot(x, data=pd, prepanel=prepanel, panel=panel,
                      frames=data, channel.x=channel.x,
                      channel.y=channel.y, channel.x.name=channel.x.name,
                      channel.y.name=channel.y.name, xlab=xlab, ylab=ylab,
                      smooth=smooth, gp=this.par.settings, as.table=as.table, filter=filter,
                      par.settings=this.par.settings, axis = axis, between = between, ...)
          
      }



## Prepanel function to set up dimensions. We want to use the instrument measurement
## range instead of the absolute range of the data. We also record the data ranges
## in the internal state environment for further use.
prepanel.xyplot.flowset <- 
    function(x, frames, channel.x.name, channel.y.name
      , xlim , ylim, ...)
{
  if (length(nm <- as.character(x)) > 1)
    stop("must have only one flow frame per panel")

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
panel.xyplot.flowset <- function(x,
                                 frames,
                                 filter = NULL,
                                 channel.x,
                                 channel.y,
								 overlay= NULL#a flowset
                                ,stats = FALSE
						 ,...)
{
    nm <- as.character(x)
    if (length(nm) < 1) return()
   
    filter <- .processFilter(filter, nm)

	
	if(!is.null(overlay))
	{
		overlay.x <- flowViz:::evalInFlowFrame(channel.x, overlay[[nm]])
		overlay.y <- flowViz:::evalInFlowFrame(channel.y, overlay[[nm]])
	}else
	{
		overlay.x <- NULL
		overlay.y <-NULL
	}
	
    if(!is.list(stats)){
        stats <- lapply(seq_along(nm), function(x) stats)
        names(stats) <- nm
      
    }
#    browser()
    panel.xyplot.flowframe(frame=frames[[nm]], filter=filter[[nm]]
                              , overlay.x=overlay.x,overlay.y=overlay.y
                              , stats = stats[[nm]]
                               , ...)
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



