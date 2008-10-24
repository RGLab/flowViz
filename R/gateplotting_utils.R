## Helper functions that will be needed for all of the gate plotting functions
## These are mostly functions that do some sanity checking.

## Store state info in this internal environment
flowViz.state <- new.env(hash = FALSE)
flowViz.state[["plotted"]] <- FALSE
flowViz.state[["par"]] <- get.gpar()

## FIXME: We hardcode a lattice theme for the X11xcairo device for now which will be
## used for all other devices as well. Later this should follow the example from the
## lattice package, i.e., real device-specific themes that are set up once the device
## is first called.
flowViz.state[["lattice.theme"]] <-
    list(X11cairo=list(gate=list(alpha=1,
                                 cex=NULL,
                                 pch=NULL,
                                 col="red",
                                 fill="transparent",
                                 lwd=1,
                                 lty="solid"),
                       gate.text=list(font=1,
                                      col="#000000",
                                      alpha=0.4,
                                      cex=1.2,
                                      lineheight=1.2),
                       flow.symbol=list(alpha=1,
                                        cex=0.8,
                                        pch=".",
                                        col="black",
                                        fill="transparent"),
                       gate.density=list(alpha=1,
                                         fill="#FFFFFFB3",
                                         col="black",
                                         lwd=1,
                                         lty="dotted")))
                                      


## return the state information from the internal environment
state <- function(x) flowViz.state[[x]]

## set the state information in the internal environment from a named
## vector or list
setState <- function(x)
{
    if(is.null(names(x)))
        stop("'x' must be a named vector or list")
    x <- as.list(x)
    for(i in seq_along(x))
        flowViz.state[[names(x)[i]]] <- x[[i]]
}


## set or return graphical parameters from the internal environment
flowViz.par <- function(x){
    if(missing(x))
        return(flowViz.state[["par"]])
    if(!is.null(names(flowViz.state[["par"]])))
        flowViz.state[["par"]][x]
    else
        NULL
}
set.flowViz.par <- function(x)
   flowViz.state[["par"]] <- modifyList(flowViz.state[["par"]], x)


## record the type of plot or the plotting dimensions in the internal
## environment
plotType <- function(type, parms){
    flowViz.state[["plotted"]] <- TRUE
    flowViz.state[["type"]] <- type
    flowViz.state[["parameters"]] <- parms
    return(invisible(NULL))
}
plotLims <- function(xlim, ylim){
    if(!missing(xlim))
       flowViz.state[["xlim"]] <- xlim
    if(!missing(ylim))
        flowViz.state[["ylim"]] <- ylim
    return(invisible(NULL))
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

checkParameterMatch <- function(parms, channels, verbose=TRUE, strict=TRUE, ...)
{
    if(missing(channels)){
        if(state("plotted") && length(state("parameters")==2)){
            sparms <- state("parameters")
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
    for(i in seq_along(x)){
        y <- x[i]
        if(is.infinite(y) && y<0)
            x[i] <- replacement[1]-diff(replacement)*10
        else if(is.infinite(y) && y>0)
            x[i] <- replacement[2]+diff(replacement)*10
    }
    return(x)
}

## convert a norm2Filter into a polygonGate
norm2Polygon <- function(fd, parms)
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


## convert an ellipseoidalFilter into a polygonGate
ell2Polygon <- function(fd, parms)
{
    ## get the ellipse lines
    center <- fd@mean[parms]
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
    opar <- flowViz.par()
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
    lPars <- c(flowViz.state[["lattice.theme"]]$X11cairo, trellis.par.get(name))
    if (is.null(name)) 
        lPars
    else if (name %in% names(lPars)) 
        lPars[[name]]
    else NULL
}


flowViz.par.set <- function (name, value, ..., theme, warn = TRUE, strict = FALSE) 
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
    fvPars <- names(theme) %in% names(flowViz.state[["lattice.theme"]]$X11cairo)
    trellis.par.set(theme=theme[!fvPars], warn=warn, strict=strict)
    if (strict)
        flowViz.state[["lattice.theme"]]$X11cairo[names(theme[fvPars])] <- theme[fvPars]  
    else flowViz.state[["lattice.theme"]]$X11cairo <-
        lattice:::updateList(flowViz.state[["lattice.theme"]]$X11cairo, theme[fvPars])
    invisible()
}
