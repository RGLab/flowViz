## ==========================================================================
## Draw points within a gate  on an existing plot. All graphical
## parameters will be passed on to graphics::points in the end, so
## these are basically extended points methods for gates and filters.
## We need to evaluate a filter to plot points (unless we pass a
## filterResult) and we always need the data to plot.
## FIXME: Need to figure out how to make points a generic without masking
##        the default function

## Default function to plot points for a single gate region. This uses the
## existing Subset architecture, hence it only works for filters that
## produce logicalFilterResults
addPoints <- function(x, data, channels, verbose=TRUE,
                      filterResult=NULL, ...)
{
    parms <- parameters(x)
    channels <- checkParameterMatch(channels, verbose=verbose,...)
    ## We check if the filterResult matches the filter and subset with that
    if(!is.null(filterResult)){
        fd <- filterDetails(filterResult, identifier(x))
        if(!identical(identifier(x), identifier(filterResult)) ||
           class(x) != class(fd$filter))
            stop("The 'filterResult' and the filter object ",
                 "don't match.", call.=FALSE)
        x <- filterResult
    }
    points(exprs(Subset(data, x))[,channels], ...)
}




## ==========================================================================
## for all filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We only know how to add points within the gate if the definiton of that gate
## fits the parameters of the data provided. We also need to know about the
## plotted parameters, or guess if there are only two in the gate.
## If the first argument is a filterResult, we extract the filter definiton
## and pass that on to the next method as a separate argument
#' Adding points within a gate to a plot
#' 
#' These methods extend the basic graphics \code{\link{points}} methods for
#' drawing of points contained within a
#' \code{\link[flowCore:filter-class]{filter}}. They allow for multiple
#' dispatch, since not all \code{\link[flowCore:filter-class]{filter}} types
#' need to be evaluated for plotting, but this decision should be made
#' internally. In any case, we need the raw data in the form of a
#' \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' 
#' When plotting \code{\link[flowCore:flowFrame-class]{flowFrame}}s using the
#' \code{plot} method provided by \code{flowViz}, the plotted parameters are
#' recorded, which makes it possible to correctly overlay the points within
#' \code{\link[flowCore:filter-class]{filter}}s assuming that they are defined
#' for the respective parameters. Warnings and error will be cast for the cases
#' where the parameters are non-distinct or ambigious.
#' 
#' @name gpoints-methods
#' @aliases gpoints gpoints gpoints-methods
#' gpoints,curv1Filter,flowFrame,character-method
#' gpoints,curv2Filter,flowFrame,character-method
#' gpoints,filter,flowFrame,missing-method gpoints,filter,missing,ANY-method
#' gpoints,filterResult,flowFrame,character-method
#' gpoints,kmeansFilter,flowFrame,character-method
#' gpoints,norm2Filter,flowFrame,character-method
#' gpoints,polygonGate,flowFrame,character-method
#' gpoints,quadGate,flowFrame,character-method
#' gpoints,rectangleGate,flowFrame,character-method
#' @docType methods
#' @section Methods:
#' 
#' \describe{
#' 
#' \item{x = "filter", data = "flowFrame", channels = "missing"}{ General
#' method for all objects inheriting from
#' \code{\link[flowCore:filter-class]{filter}}. This is used as the default
#' when no more explicit method is found. It tries to find the plotted
#' parameters from the internal \code{flowViz.state} environment. This only
#' works if the flow data has been plotted using the \code{plot} methods
#' provided by this \code{flowViz} package. }
#' 
#' \item{x = "filter", data = "missing", channels = "ANY"}{ This gives a useful
#' error message when we don't get what we need. }
#' 
#' \item{x = "filterResult", data = "flowFrame", channels = }{ We can get all
#' the information about a \code{\link[flowCore:filter-class]{filter}} from its
#' \code{\link[flowCore:filterResult-class]{filterResult}} without the need to
#' re-evaluate.}\item{ "character"}{ We can get all the information about a
#' \code{\link[flowCore:filter-class]{filter}} from its
#' \code{\link[flowCore:filterResult-class]{filterResult}} without the need to
#' re-evaluate.}
#' 
#' \item{x = "curv1Filter", data = "ANY"}{ We either need a
#' \code{\link[flowCore:filterResult-class]{filterResult}} or the raw data as a
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} for
#' \code{\link[flowStats:curv1Filter-class]{curv1Filter}}s. }
#' 
#' \item{x = "curv1Filter", data = "flowFrame"}{ see above }
#' 
#' \item{x = "curv1Filter", data = "missing"}{ see above }
#' 
#' \item{x = "curv1Filter", data = "multipleFilterResult"}{ see above }
#' 
#' \item{x = "curv2Filter", data = "ANY"}{ We either need a
#' \code{\link[flowCore:filterResult-class]{filterResult}} or the raw data as a
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} for
#' \code{\link[flowStats:curv2Filter-class]{curv2Filter}}s.}
#' 
#' \item{x = "curv1Filter", data = "flowFrame", channels = }{ We evaluate the
#' \code{\link[flowCore:filter-class]{filter}} on the
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} and plot the subset of
#' selected points. By default, every subpopulation (if there are any) is
#' colored differently.}\item{ "character"}{ We evaluate the
#' \code{\link[flowCore:filter-class]{filter}} on the
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} and plot the subset of
#' selected points. By default, every subpopulation (if there are any) is
#' colored differently.}
#' 
#' \item{x = "curv2Filter", data = "flowFrame", channels = "character"}{ see
#' above }
#' 
#' \item{x = "kmeansFilter", data = "flowFrame", channels = }{ see above
#' }\item{ "character"}{ see above }
#' 
#' \item{x = "norm2Filter", data = "flowFrame", channels = }{ see above }\item{
#' "character"}{ see above }
#' 
#' \item{x = "polygonGate", data = "flowFrame", channels = }{ see above }\item{
#' "character"}{ see above }
#' 
#' \item{x = "quadGate", data = "flowFrame", channels = "character"}{ see above
#' }
#' 
#' \item{x = "rectangleGate", data = "flowFrame", channels = }{ see above
#' }\item{ "character"}{ see above }
#' 
#' }
#' @author F. Hahne
#' @seealso
#' 
#' \code{\link[flowCore:filter-class]{filter}},
#' \code{\link[flowCore:flowFrame-class]{flowFrame}}, \code{\link{glines}},
#' \code{\link{gpolygon}}
#' @keywords methods
#' @importFrom graphics abline axis barplot contour hist layout lines pairs par points polygon rect smoothScatter
#' @export
#' @param x filter, filterResult or any derived filter class
#' @param data flowFrame or missing
#' @param channels character or missing
#' @param verbose logical
#' @param filterResult filterResult class
#' @param ... other arguments
setMethod("gpoints",
          signature(x="filter", data="flowFrame", channels="missing"), 
          function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          if(is(x, "filterResult")){
              filterResult <- x
              fd <- filterDetails(x)
              x <- fd[[length(fd)]]$filter
          }
          channels <- checkParameterMatch(parameters(x), verbose=verbose,...)
          gpoints(x=x, data=data, channels=channels, verbose=verbose,
                  filterResult=filterResult, ...)
      })

## Need this to deal with filterResults when the plotting channels are
## specified
setMethod("gpoints",
          signature(x="filterResult", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          if(x@frameId != identifier(data))
              stop("The filter was evaluated on flowFrame '",
                   x@frameId, "'\n  but the frame provided is '",
                   identifier(data), "'.", call.=FALSE)
          filterResult <- x
          fd <- filterDetails(x)
          x <- fd[[length(fd)]]$filter
          channels <- checkParameterMatch(channels, verbose=verbose,...)
          gpoints(x=x, data=data, channels=channels, verbose=verbose,
                   filterResult=filterResult, ...)
      })    

## A useful error message when we don't get what we need
setMethod("gpoints",
          signature(x="filter", data="missing", channels="ANY"), 
          function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          stop("Need the 'flowFrame' in order to add points.", call.=FALSE)
      })



          
## ==========================================================================
## for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Only warning that we don't use the filterResult, the rest is pretty
## much the default
setMethod("gpoints",
          signature(x="rectangleGate", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      { 
          if(!is.null(filterResult))
              dropWarn("filterResult", "rectangleGates", verbose=verbose)
          addPoints(x=x, data=data, channels=channels, verbose=verbose,
                      filterResult=NULL, ...)
      })




## ==========================================================================
## for quadGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We plot this as four individual rectangle gates
setMethod("gpoints",
          signature(x="quadGate", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE, col, ...)
      {
          if(!is.null(filterResult))
              dropWarn("filterResult", "quadGates", verbose=verbose)
          parms <- parameters(x)
          channels <- checkParameterMatch(channels, verbose=verbose,...)
          if(missing(col))
              col <-  colorRampPalette(brewer.pal(9, "Set1"))(4)
          else
              col <- rep(col,4)
          v <- x@boundary[channels[1]]
          h <- x@boundary[channels[2]]
          mat <- matrix(c(-Inf, v, h, Inf, v, Inf, h, Inf, -Inf, v, -Inf,
                          h, v, Inf, -Inf, h), byrow=TRUE, ncol=4)              
          for(i in 1:4){
              rg <- rectangleGate(.gate=matrix(mat[i,], ncol=2,
                                  dimnames=list(c("min", "max"), parms)))
              gpoints(rg, data=data, channels=channels, verbose=FALSE,
                      col=col[i], ...)
          }
      })




## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Only warning that we don't use the filterResult, the rest is pretty
## much the default
setMethod("gpoints",
          signature(x="polygonGate", data="flowFrame", channels="character"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          if(!is.null(filterResult))
              dropWarn("filterResult", "polygonGates", verbose=verbose)
          addPoints(x=x, data=data, channels=channels, verbose=verbose,
                      filterResult=NULL, ...)
      })


## ==========================================================================
## for norm2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We don't need to reevaluate the filter if the filterResult is used,
## so no need for a warning here
setMethod("gpoints",
          signature(x="norm2Filter", data="flowFrame", channels="character"), 
          definition=function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          addPoints(x=x, data=data, channels=channels, verbose=verbose,
                      filterResult=filterResult, ...)
      })



## ==========================================================================
## for kmeansFilters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## This filter produces a multipleFilterResult, so we can't subset directly.
## Instead, we split the original frame and plot each component separately
setMethod("gpoints",
          signature(x="kmeansFilter", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE,
                   filterResult=NULL, col, ...)
      {
          ## We check that the filterResult matches the filter and split by that
          channels <- checkParameterMatch(channels, verbose=verbose,...)
          if(!is.null(filterResult)){
              if(!identical(identifier(x), identifier(filterResult)) ||
                 class(x) != class(filterDetails(filterResult,
                                                 identifier(x))$filter))
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
          mapply(function(z, co, ...) points(exprs(z)[,channels], col=co, ...),
                 z=datsplit, co=col, MoreArgs=list(...))
          return(invisible(NULL))
      })







          
