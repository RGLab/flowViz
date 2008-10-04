## Helper functions that will be needed for all of the gate plotting functions
## These are mostly functions that do some sanity checking.

## Store state info in this internal environment
flowViz.state <- new.env(hash = FALSE)
flowViz.state[["plotted"]] <- FALSE
flowViz.state[["par"]] <- get.gpar()


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


## record the type of plot or the limits in the internal environment
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

checkParameterMatch <- function(parms, channels, verbose=TRUE, ...)
{
    if(missing(channels)){
        if(state("plotted") && length(state("parameters")==2)){
            sparms <- state("parameters")
            err <-   function()
                stop("The flow parameters used in the last plot don't match ",
                     "the\nparameters provided by the filter or via the ",
                     "'channels' argument.\n   Plotted: ",
                     paste(sparms, collapse=" vs. "), "\n   Provided: ",
                     paste(parms, collapse=", "),
                     "\nPlease check or redraw the plot.", call.=FALSE)
            mt <- sparms %in% parms
            if(length(parms)<2){
                if(!parms %in% sparms)
                    err()
            }else{
                if(!all(mt))
                    err()
            }
            return(sparms)
        }else{
            fmatchWarn(parms, verbose=verbose)
            return(parms)
        }
    }else{
        parms <- channels
        if(length(parms)!=2)
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
    if(!identical(identifier(filter), identifier(fres)) ||
       class(filter) != class(fd$filter))
        stop("The 'filterResult' and the '", class(filter),
             "' don't match.", call.=FALSE)
    if("transformation" %in% slotNames(fd$filter) &&
       length(fd$filter@transformation) > 0 && verbose)
        warning("'result' appears to have been applied on ",
                "transformed data.\nThese are not supported yet.")
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
            x[i] <- replacement[1]
        else if(is.infinite(y) && y>0)
            x[i] <- replacement[2]
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
    polygonGate(boundaries=ans)
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
    polygonGate(boundaries=ans)
}


## A helper function that calls grid.rect in a more user-fiendly way
glrect <- function (xleft, ybottom, xright, ytop, x=(xleft + xright)/2, 
                    y=(ybottom + ytop)/2, width=xright - xleft,
                    height=ytop - ybottom, just="center",
                    hjust=NULL, vjust=NULL, gp, ...) 
{
    class(gp) <- "gpar"
    grid.rect(x=x, y=y, width=width, height=height, default.units="native", 
        just=just, hjust=hjust, vjust=vjust, gp=gp)
}


## Default function to plot points for a single gate region. This uses the
## existing Subset architecture, hence it only works for filters that
## produce logicalFilterResults so far
addLpoints <- function(x, data, channels, verbose=TRUE,
                      filterResult=NULL, gpar, pch=".", ...)
{
    parms <- parameters(x)
    channels <- checkParameterMatch(channels, verbose=verbose)
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
    if(!is.null(gpar))
        opar <- modifyList(opar, gpar)
    class(opar) <- "gpar"
    grid.points(exp[,channels[1]], exp[,channels[2]],
               default.units = "native", gp=opar, pch=pch)
}


## Helper function to add points for multipleFilterResults.
## We check that the filterResult matches the filter and split by that
## The first component is everything outside not in the filter and we
## drop that
multFiltPoints <-  function(x, data, channels, verbose=TRUE,
                            filterResult=NULL, pch=".", gpar, ...)
{

    channels <- checkParameterMatch(channels, verbose=verbose)
    if(!is.null(filterResult)){
        checkIdMatch(x=x, f=data)
        x <- filterResult
    }
    datsplit <- split(data, x)[-1]
    ld <- length(datsplit)
    ## we want to be able to use different colors for each population
    opar <- flowViz.par()
    col <- 2:(ld+1)
    if(!missing(gpar) && !is.null(gpar)){
        if("col" %in% names(gpar))
            col <- rep(gpar$col, ld)
        opar <- modifyList(opar, gpar)
    }
    class(opar) <- "gpar"
    pch <- rep(pch, ld)
    for(i in 1:ld){
        opar$col <- col[i]
        grid.points(exprs(datsplit[[i]])[,channels[1]],
                    exprs(datsplit[[i]])[,channels[2]],
                    default.units = "native", gp=opar,
                    pch=pch[i])
    }
    return(invisible(NULL))
}
