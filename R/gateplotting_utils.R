## Helper functions that will be needed for all of the gate plotting functions
## This is mostly function that do some sanity checking.


## return the state information from the internal environment
state <- function(x) flowViz.state[[x]]

## return the state information in the internal environment from a named
## vector or list
setState <- function(x)
{
    if(is.null(names(x)))
        stop("'x' must be a named vector or list")
    x <- as.list(x)
    for(i in seq_along(x))
        flowViz.state[[names(x)[i]]] <- x[[i]]
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

checkParameterMatch <- function(parms, verbose=TRUE)
{
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
    }else if(length(parms)!=2)
        stop("The filter definition contains the following parameters:\n",
             paste(parms, collapse=", "), "\nDon't know how to match to",
             " the plotted data.\nPlease specify plotting parameters as ",
             "an additional argument.", call.=FALSE)
    else{
        fmatchWarn(parms=parms, verbose=verbose)
        return(parms)
    }
}


## We check that the IDs of a flowFrame and a filter are matching
checkIdMatch <- function(x, f)
{
    if(x@frameId != identifier(f))
        stop("The filter was evaluated on flowFrame '",
             x@frameId, "'\n  but the frame provided is '",
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
