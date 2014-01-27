## Helper functions that will be needed for all of the gate plotting functions
## These are mostly functions that do some sanity checking.
#require(grDevices)  # RColorBrewer assumes 'rgb' on the search path

## Store state info in this internal environment
flowViz.state <- new.env(hash = FALSE)


# based on colorRampPalette
# add the optional alpha control
#' @param alpha \code{numeric} range from 0 to 1
.colRmpPlt <- function (alpha = 1, ...) 
{
  
  colors <- rev(brewer.pal(11, "Spectral"))
  ramp <- colorRamp(colors, ...)
  function(n) {
    x <- ramp(seq.int(0, 1, length.out = n))
    rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255, alpha = alpha * 255)
  }
}

## FIXME: We hardcode a lattice theme for the X11xcairo device for now which will be
## used for all other devices as well. Later this should follow the example from the
## lattice package, i.e., real device-specific themes that are set up once the device
## is first called.
.flowViz.par.init <- function(lattice.par){
  gate.par <- list(gate=list(alpha=1,
                              cex=NULL,
                              pch=NULL,
                              col="#9E0142"#"red",
                              ,fill="transparent",
                              lwd=1,
                              lty="solid"),
                        gate.text=list(font=1,
                                        col="#000000",
                                        alpha=1,#0.4,
                                        cex=0.8,#1.2,
                                        lineheight=0.8#1.2
                                        ,background=list(fill="white"
                                                        ,col="transparent"
                                                        ,alpha=1
                                                    )
                                    ),
                        overlay.symbol = list(alpha = 0.5
                                              ,bg.alpha = 0.3
                                              ,col = "transparent"
                                              ,fill = "red"
                                              ,cex = 0.5
                                              , pch = 19
                                          ),                 
                        flow.symbol=list(alpha=1,
                                            cex=0.8,
                                            pch=".",
                                            col="black",
                                            fill="transparent"),
                        gate.density=list(alpha=1,
                                          fill="#FFFFFFB3",
                                          col="black",
                                          lwd=1,
                                          lty="dotted")
                                      
                        ,argcolramp = .colRmpPlt()
                    )
  flowViz.state[["lattice.theme"]] <- list(X11cairo= c(gate.par, lattice.par))
}
# can't do it within Namespace events hooks (e.g. .onLoad) 
# since it needs to load namespace (e.g. 'rgb' through 'brewer.pal' call
# other than base namespace
.flowViz.par.init(lattice.par = ggplot2like())

## Get and set graphical defaults for plots in flowViz. To a great extend, this uses
## the trellis.par.get and trellis.par.set infrastructure defined in the lattice package.
## Currently it is not possible to directly register additional setting there, instead
## we provide our own wrappers which look for available defaults in lattice first and
## fetch or set additional defaults within the internal flowViz environment if they can't
## be found in the lattice defaults.
flowViz.par.get <- function (name = NULL) 
{
  ## FIXME: We want this to be device-specific, needs to be set up in the same
  ## way as it is done in lattice.
  lPars <- flowViz.state[["lattice.theme"]]$X11cairo
  if (is.null(name)) 
    lPars
  else if (name %in% names(lPars)) 
    lPars[[name]]
  else NULL
}


flowViz.par.set <- function (name, value, ..., theme, warn = TRUE, strict = FALSE, reset = FALSE) 
{
  if (missing(theme)) {
    if (!missing(value)) {
      theme <- list(value)
      names(theme) <- name
    }
    else if (!missing(name) && is.list(name)) {
      theme <- name
    }
    else theme <- list(...)
  }
  else {
    if (is.character(theme)) 
      theme <- get(theme)
    if (is.function(theme)) 
      theme <- theme()
    if (!is.list(theme)) {
      warning("Invalid 'theme' specified")
      theme <- NULL
    }
  }
  
  #reset by dropping all the existing lattic par
  if(reset){
    .flowViz.par.init(lattice.par = theme)  
  }else{
    fvPars <- names(theme) %in% names(flowViz.state[["lattice.theme"]]$X11cairo)
    
    if (strict)
      flowViz.state[["lattice.theme"]]$X11cairo[names(theme[fvPars])] <- theme[fvPars]  
    else flowViz.state[["lattice.theme"]]$X11cairo <-
          lattice:::updateList(flowViz.state[["lattice.theme"]]$X11cairo, theme[fvPars])
    
  }
  invisible()
}




## return the state information from the internal environment
state <- function(x) flowViz.state[[x]]

plotType <- function(type, parms){
	return (list(plotted= TRUE
				,type=type
				,parameters= parms)
			)
}




## We only know how to add the gate boundaries if the definiton of that gate
## contains exactly two dimensions, and even then we need to guess that
## they match the plotted data, hence we warn
fmatchWarn <- function(parms, verbose=TRUE)
{
    if(verbose)
        warning("The filter is defined for parameters '",
                paste(parms, collapse="' and '"), "'.\nPlease make sure ",
                "that they match the plotting parameters.", call.=FALSE)
}

checkParameterMatch <- function(parms,channels, verbose=TRUE, strict=TRUE,ptList, ...)
{
	parameters<-ptList$parameters
	plotted<-ptList$plotted
    if(missing(channels)){
        if(plotted && length(parameters==2)){
            sparms <- parameters
            err <-   function(strict=TRUE){
                if(strict)
                    stop("The flow parameters used in the last plot don't match ",
                         "the\nparameters provided by the filter or via the ",
                         "'channels' argument.\n   Plotted: ",
                         paste(sparms, collapse=" vs. "), "\n   Provided: ",
                         paste(parms, collapse=", "),
                         "\nPlease check or redraw the plot.", call.=FALSE)
            }
            mt <- sparms %in% parms
            if(length(parms)<2){
                if(!parms %in% sparms){
                    err(strict)
                    return(NA)
                }
            }else if(!all(mt)){
                err(strict)
                return(NA)
            }
            return(sparms)
        }else{
            fmatchWarn(parms, verbose=verbose)
            return(parms)
        }
    }else{
        parms <- channels
        if(length(parms)!=2 && strict)
            stop("The filter definition contains the following parameters:\n",
                 paste(parms, collapse=", "), "\nDon't know how to match to",
                 " the plotted data.\nPlease specify plotting parameters as ",
                 "an additional argument.", call.=FALSE)
        return(parms)
    }
}


## We check that the IDs of a flowFrame and a filter are matching
checkIdMatch <- function(x, f)
{
    if(!identifier(f) %in% x@frameId)
        stop("The filter was evaluated on flowFrame(s)\n'",
             paste(x@frameId, collapse="', '", sep=),
             "'\n  but the frame provided is '",
             identifier(f), "'.", call.=FALSE)
}


## Warning when we ignore an argument that is not needed
dropWarn <- function(type, gate, verbose=FALSE)
{
    if(verbose)
        warning("No '", type, "' needed to plot '", gate,
                "'.\nArgument is ignored.", call.=FALSE)
}


## Cast error that we need the data to evaluate the filter
evalError <- function(type)
{
    stop("'", type, "' need to be evaluated for plotting.\n",
         "Either provide a 'flowFrame' or an appropriate ",
         "'filterResult' \nas second argument.", call.=FALSE)
}

## check that a filterResult and a filter actually match
checkFres <- function(filter, fres, verbose=TRUE)
{
    fd <- filterDetails(fres, identifier(filter))
    if(!identical(identifier(filter), identifier(fres)))
        stop("The 'filterResult' and the '", class(filter),
             "' don't match.", call.=FALSE)
}

## We don't know how to draw 'type' filters, hence we warn
nnWarn <- function(type, verbose=TRUE)
{
    if(verbose)
        warning("Don't know how to plot outline for a '", type,
                "'.", call.=FALSE)
    return(invisible(NULL))
}

## replace infinite values by something useful
fixInf <- function(x, replacement)
{
#	browser()
    for(i in seq_along(x)){
        y <- x[i]
        if(is.infinite(y) && y<0)
            x[i] <- replacement[1]-diff(replacement)*10
        else if(is.infinite(y) && y>0)
            x[i] <- replacement[2]+diff(replacement)*10
    }
    return(x)
}
##fixInf returned the exagerated bounaries which does not give good esitmation of label position for gates
##so we add another version here specifically for addName methods
fixBound_addName <- function(x, range)
{
#	browser()
	for(i in seq_along(x)){
		y <- x[i]
		if(y<range[1])
			x[i] <- range[1]
		else if(y>range[2])
			x[i] <- range[2]
	}
	return(x)
}

## convert a norm2Filter into a polygonGate
norm2Polygon <- function(fd, parms=fd$parameters)
{
    ## get the ellipse lines
    norm.center <- fd$center[parms]
    norm.cov <- fd$cov[parms, parms]
    norm.radius <- fd$radius
    chol.cov <- t(chol(norm.cov))
    t <- seq(0, 2 * base::pi, length = 50)
    ans <- norm.center +
        (chol.cov %*% rbind(x = norm.radius * cos(t),
                            y = norm.radius * sin(t)))
    ans <- as.data.frame(t(ans))
    names(ans) <- parms
    ## create a polygonGate
    polygonGate(.gate=ans)
}

## convert a norm2Gate into an ellipsoidGate
norm2Ell <- function(fd, parms=fd$parameters)
{
    ellipsoidGate(mean=fd$center[parms],
                  cov=fd$cov[parms, parms],
                  distance=fd$radius)
}

## convert an ellipseoidalFilter into a polygonGate
ell2Polygon <- function(fd, parms=parameters(fd))
{
    ## get the ellipse lines
    center <- fd@mean[parms]
    if(is.null(rownames(fd@cov)))
      rownames(fd@cov) <- colnames(fd@cov)
    cov <- fd@cov[parms, parms]
    radius <- fd@distance
    chol.cov <- t(chol(cov))
    t <- seq(0, 2 * base::pi, length = 50)
    ans <- center +
        (chol.cov %*% rbind(x = radius * cos(t),
                            y = radius * sin(t)))
    ans <- as.data.frame(t(ans))
    names(ans) <- parms
    ## create a polygonGate
    polygonGate(.gate=ans)
}


## Helper function that calls panel.x in a more user-fiendly way
glrect <- function (xleft, ybottom, xright, ytop, ..., gp) 
{
    panel.rect(xleft, ybottom,  xright, ytop, col=gp$fill,
               border=gp$col, lty=gp$lty, lwd=gp$lwd,
               alpha=gp$alpha, ...)
}

gltext <- function (x, y, labels, ..., gp) 
{   
#    browser()
#    gp_txt <- gp
#    gp_txt$background <- NULL
#    txtObj <- textGrob(label = labels, gp = gpar(gp_txt), default.units = "native")
    cex <- gp$cex    
	#add rectange as white background for the better visual effect when label is plotted against color 
	grid.rect(x=unit(x,"native")
			,y=unit(y,"native")
#			,width= unit(cex, "grobwidth", data = txtObj)#unit(1,'strwidth',labels)
#			,height=unit(cex, "grobheight", data = txtObj) #unit(1,'strheight',labels)
            ,width= unit(cex,'strwidth',labels)
            ,height=unit(cex,'strheight',labels)
			, gp=gpar(fill=gp$background$fill
					,col=gp$background$col
					,alpha=gp$background$alpha
					)
	)
#	browser()
		
    panel.text(x, y, labels=labels, col=gp$col, cex=gp$cex,
               linheight=gp$lineheigt, alpha=gp$alpha, font=gp$font)
}

glpoly <- function (x, y, ..., gp) 
{
    panel.polygon(x, y, labels=labels, border=gp$col, col=gp$fill,
                  lty=gp$lty, lwd=gp$lwd, alpha=gp$alpha, ...)
}





## Default function to plot points for a single gate region. This uses the
## existing Subset architecture, hence it only works for filters that
## produce logicalFilterResults so far
addLpoints <- function(x, data, channels, verbose=TRUE,
                       filterResult=NULL, gpar, ...)
{
    ## We check if the filterResult matches the filter and subset with that
    if(!is.null(filterResult)){
        if(!identical(identifier(x), identifier(filterResult)) ||
           class(x) != class(filterDetails(filterResult)[[1]]$filter))
            stop("The 'filterResult' and the filter object ",
                 "don't match.", call.=FALSE)
        x <- filterResult
    }
    exp <- exprs(Subset(data, x))
    
    panel.points(exp[,channels[1]], exp[,channels[2]],
                 pch=gpar$pch, cex=gpar$cex, col=gpar$col,
                 fill=gpar$fill, alpha=gpar$alpha, ...)
    return(invisible(NULL))
}


## Helper function to add points for multipleFilterResults.
## We check that the filterResult matches the filter and split by that
## The first component is everything outside not in the filter and we
## drop that
multFiltPoints <-  function(x, data, channels, verbose=TRUE,
                            filterResult=NULL, gpar, ...)
{
    if(is.null(filterResult))
        filterResult <- filter(data, x)
    checkIdMatch(x=filterResult, f=data)
    x <- filterResult 
    datsplit <- split(data, x)[-1]
    ld <- length(datsplit)
    ## we want to be able to use different colors for each population
    col <- rep(gpar$col, ld)[1:ld]
    for(i in 1:ld){
        gpar$col <- col[i]
        panel.points(exprs(datsplit[[i]])[,channels[1]],
                     exprs(datsplit[[i]])[,channels[2]],
                     pch=gpar$pch, cex=gpar$cex,
                     col=gpar$col, alpha=gpar$alpha,
                     fill=gpar$fill, ...)
    }
    return(invisible(NULL))
}





