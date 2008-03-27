## ==========================================================================
## Draw points within a gate  on an existing plot. All graphical
## parameters will be passed on to lpoints in the end, so
## these are basically extended lpoints methods for gates and filters.
## We need to evaluate a filter to plot points (unless we pass a
## filterResult) and we always need the data to plot.
## FIXME: Need to figure out how to make points a generic without masking
##        the default function

## Default function to plot points for a single gate region. This uses the
## existing Subset architecture, hence it only works for filters that
## produce logicalFilterResults
addLpoints <- function(x, data, channels, verbose=TRUE,
                      filterResult=NULL, col, ...)
{
    parms <- parameters(x)
    if(length(channels) != 2)
        stop("Plotting parameters need to be character vector of ",
             "length 2", call.=FALSE)
    mt <- match(parms, channels)
    if(any(is.na(mt)))
        stop("The filter is not defined for the following ",
             "parameter(s):\n", paste(channels[is.na(mt)], collapse=", "),
             call.=FALSE)
    mt2 <- match(channels, colnames(data)) 
    if(any(is.na(mt2)))
        stop("The gate definition includes parameters\n",
             paste(parms, collapse=", "), "\nbut the data only ", 
             "provides\n", paste(colnames(data), collapse=", "),
             call.=FALSE)
    ## We check if the filterResult matches the filter and subset with that
    if(!is.null(filterResult)){
        if(!identical(identifier(x), identifier(filterResult)) ||
           class(x) != class(filterDetails(filterResult)[[1]]$filter))
            stop("The 'filterResult' and the filter object ",
                 "don't match.", call.=FALSE)
        x <- filterResult
    }
    if(missing(col))
        col="red"
    lpoints(exprs(Subset(data, x))[,channels], col=col[1], ...)
}



## ==========================================================================
## for all filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We only know how to add points within the gate if the definiton of that gate
## fits the parameters of the data provided. We also need to know about the
## plotted parameters, or guess if there are only two in the gate.
## If the first argument is a filterResult, we extract the filter definiton
## and pass that on to the next method as a separate argument
setMethod("glpoints",
          signature(x="filter", data="flowFrame", channels="missing"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          if(is(x, "filterResult")){
              filterResult <- x
              x <- filterDetails(x)[[1]]$filter
          }
          parms <- parameters(x)
          mt <- match(parms, colnames(data))
          if(length(parms)!=2)
              stop("The filter definition contains the following parameters:\n",
                   paste(parms, collapse=", "), "\nDon't know how to match to",
                   " the plotted data.\nPlease specify plotting parameters as ",
                   "an additional argument.", call.=FALSE)
          if(any(is.na(mt)))
              stop("The gate definition includes parameters\n",
                   paste(parms, collapse=", "), "\nbut the data only ", 
                   "provides\n", paste(colnames(data), collapse=", "),
                   call.=FALSE)
          if(verbose)
              warning("The filter is defined for parameters '",
                      paste(parms, collapse="' and '"), "'.\nPlease make sure ",
                      "that they match the plotting parameters.", call.=FALSE)
          glpoints(x=x, data=data, channels=parms, verbose=verbose,
                  filterResult=filterResult, ...)
      })

## Need this to deal with filterResults when the plotting channels are
## specified
setMethod("glpoints",
          signature(x="filterResult", data="flowFrame", channels="character"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          filterResult <- x
          x <- filterDetails(x)[[1]]$filter
          glpoints(x=x, data=data, channels=channels, verbose=verbose,
                   filterResult=filterResult, ...)
      })    


## ==========================================================================
## for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Only warning that we don't use the filterResult, the rest is pretty
## much the default
setMethod("glpoints",
          signature(x="rectangleGate", data="flowFrame", channels="character"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          if(verbose & !is.null(filterResult))
               warning("No 'filterResult' needed to plot 'rectangleGates'.\n",
                       "Argument is ignored.", call.=FALSE)
          addLpoints(x=x, data=data, channels=channels, verbose=verbose,
                      filterResult=NULL, ...)
      })



## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Only warning that we don't use the filterResult, the rest is pretty
## much the default
setMethod("glpoints",
          signature(x="polygonGate", data="flowFrame", channels="character"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          if(verbose & !is.null(filterResult))
               warning("No 'filterResult' needed to plot 'polygonGates'.\n",
                       "Argument is ignored.", call.=FALSE)
          addLpoints(x=x, data=data, channels=channels, verbose=verbose,
                      filterResult=NULL, ...)
      })


## ==========================================================================
## for norm2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We don't need to reevaluate the filter if the filterResult is used,
## so no need for a warning here
setMethod("glpoints",
          signature(x="norm2Filter", data="flowFrame", channels="character"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          addLpoints(x=x, data=data, channels=channels, verbose=verbose,
                      filterResult=filterResult, ...)
      })




## ==========================================================================
## for curv2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## This filter produces a multipleFilterResult, so we can't subset directly.
## Instead, we split the original frame and plot each component separately
setMethod("glpoints",
          signature(x="curv2Filter", data="flowFrame", channels="character"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, col, ...)
      {
          ## We check that the filterResult matches the filter and split by that
          if(!is.null(filterResult)){
              if(!identical(identifier(x), identifier(filterResult)) ||
                 class(x) != class(filterDetails(filterResult)[[1]]$filter))
                  stop("The 'filterResult' and the filter object ",
                       "don't match.", call.=FALSE)
              x <- filterResult
          }
          datsplit <- split(data, x)[-1]
          ld <- length(datsplit)
          if(missing(col))
              col <-  colorRampPalette(brewer.pal(9, "Set1"))(ld)
          else
              col <- rep(col, ld)[1:ld]
          mapply(function(z, co, ...) lpoints(exprs(z)[,channels], col=co, ...),
                 z=datsplit, co=col, MoreArgs=list(...))
      })




## ==========================================================================
## for curv1Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## This filter produces a multipleFilterResult, so we can't subset directly.
## Instead, we split the original frame and plot each component separately
setMethod("glpoints",
          signature(x="curv1Filter", data="flowFrame", channels="character"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, col, ...)
      {
          ## We check that the filterResult matches the filter and split by that
          if(!is.null(filterResult)){
              if(!identical(identifier(x), identifier(filterResult)) ||
                 class(x) != class(filterDetails(filterResult)[[1]]$filter))
                  stop("The 'filterResult' and the filter object ",
                       "don't match.", call.=FALSE)
              x <- filterResult
          }
          datsplit <- split(data, x)[-1]
          ld <- length(datsplit)
          if(missing(col))
              col <-  colorRampPalette(brewer.pal(9, "Set1"))(ld)
          else
              col <- rep(col, ld)[1:ld]
          mapply(function(z, co, ...) lpoints(exprs(z)[,channels], col=co, ...),
                 z=datsplit, co=col, MoreArgs=list(...))
      })




## ==========================================================================
## for kmeansFilters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## This filter produces a multipleFilterResult, so we can't subset directly.
## Instead, we split the original frame and plot each component separately
setMethod("glpoints",
          signature(x="kmeansFilter", data="flowFrame", channels="character"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, col, ...)
      {
          ## We check that the filterResult matches the filter and split by that
          if(!is.null(filterResult)){
              if(!identical(identifier(x), identifier(filterResult)) ||
                 class(x) != class(filterDetails(filterResult)[[1]]$filter))
                  stop("The 'filterResult' and the filter object ",
                       "don't match.", call.=FALSE)
              x <- filterResult
          }
          datsplit <- split(data, x)
          ld <- length(datsplit)
          if(missing(col))
              col <-  colorRampPalette(brewer.pal(9, "Set1"))(ld)
          else
              col <- rep(col, ld)[1:ld]
          mapply(function(z, co, ...) lpoints(exprs(z)[,channels], col=co, ...),
                 z=datsplit, co=col, MoreArgs=list(...))
      })







          
