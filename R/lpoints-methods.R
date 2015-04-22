## ==========================================================================
## Draw points within a gate  on an existing plot. All graphical
## parameters will be passed on to lpoints in the end, so
## these are basically extended lpoints methods for gates and filters.
## We need to evaluate a filter to plot points (unless we pass a
## filterResult) and we always need the data to plot, hence the methods
## always dispatch on filters (or filteResults since
## they inherits from filter) and on a flowFrame.
## The optional channels argument can be used to pass along
## the plotting parameters, and it takes precendence over everything else.
## FIXME: Need to figure out how to make points a generic without masking
##        the default function



#' Adding points within a gate to a plot
#' 
#' These methods extend the lattice \code{\link[lattice:llines]{lpoints}}
#' methods for drawing of points contained within a
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
#' where the parameters are non-distinct or ambigious. These methods are meant
#' to be used within lattice panel functions and are probably not of much use
#' outside of those.
#' 
#' @name glpoints-methods
#' @aliases glpoints glpoints-methods
#' glpoints,curv1Filter,flowFrame,character-method
#' glpoints,curv2Filter,flowFrame,character-method
#' glpoints,filter,flowFrame,missing-method glpoints,filter,missing,ANY-method
#' glpoints,filterResult,flowFrame,character-method
#' glpoints,kmeansFilter,flowFrame,character-method
#' glpoints,norm2Filter,flowFrame,character-method
#' glpoints,polygonGate,flowFrame,character-method
#' glpoints,quadGate,flowFrame,character-method
#' glpoints,rectangleGate,flowFrame,character-method
#' glpoints,ellipsoidGate,flowFrame,character-method
#' glpoints,filterResult,flowFrame,ANY-method
#' glpoints,subsetFilter,flowFrame,ANY-method
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
#' \item{x = "ellipsoidGate", data = "flowFrame", channels = }{ see above
#' }\item{ "character"}{ see above }
#' 
#' }
#' @author F. Hahne
#' @seealso
#' 
#' \code{\link[flowCore:filter-class]{filter}},
#' \code{\link[flowCore:flowFrame-class]{flowFrame}}, \code{\link{glpolygon}}
#' @keywords methods
#' @export
#' @param x filter, filterResult or any derived filter class
#' @param data flowFrame or missing
#' @param channels character or missing
#' @param verbose logical
#' @param filterResult filterResult class
#' @param ... other arguments
setMethod("glpoints",
          signature(x="filter", data="flowFrame", channels="missing"), 
          function(x, data, channels, verbose=TRUE,
          filterResult=NULL, ...)
      {
          glpoints(x=x, data=data, channels=parms, verbose=verbose,
                   filterResult=filterResult, strict=strict, ...)
      })

## Extract the filter definiton from a filterResult and pass that on
## along with it. We decide later if we really need it or not.
setMethod("glpoints",
          signature(x="filterResult", data="flowFrame", channels="ANY"), 
          function(x, data, channels, verbose=TRUE,
                   filterResult=NULL, ...)
      {
          checkIdMatch(x=x, f=data)
          filterResult <- x
          fd <- filterDetails(x)
          x <- fd[[length(fd)]]$filter
          if(missing(channels))
              glpoints(x=x, data=data, verbose=verbose,
                       filterResult=filterResult, ...)
          else
              glpoints(x=x, data=data, channels=channels,
                       verbose=verbose, filterResult=filterResult, ...)
      })    

## A useful error message when we don't get what we need.
setMethod("glpoints",
          signature(x="filter", data="missing", channels="ANY"), 
          function(x, data, channels, ...)
          stop("Need the 'flowFrame' in order to add points.", call.=FALSE))




## ==========================================================================
## for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Only warning that we don't use the filterResult, the rest is pretty
## much the default.
setMethod("glpoints",
          signature(x="rectangleGate", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE, filterResult=NULL,
                   gpar=flowViz.par.get(), names=FALSE, ...)
      {
          if(!is.null(filterResult))
              dropWarn("filterResult", "rectangleGates", verbose=verbose)
          addLpoints(x=x, data=data, channels=channels, verbose=verbose,
                     filterResult=NULL, gpar=gpar$gate, ...)
          ## add names if necessary
          addName(x, names, channels, gpar$gate.text,...)
      })




## ==========================================================================
## for ellipsoidGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We convert to a polygon gate and pass that on
setMethod("glpoints",
          signature(x="ellipsoidGate", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE, filterResult=NULL,
                   gpar=flowViz.par.get(), names=FALSE, ...)
      {
          xe <- ell2Polygon(x, parameters(x))
          if(!is.null(filterResult))
              dropWarn("filterResult", "ellipsoidGates", verbose=verbose)
          addLpoints(x=xe, data=data, channels=channels, verbose=verbose,
                     filterResult=NULL, gpar=gpar$gate, ...)
          addName(x, names, channels, gpar$gate.text,...)
      })



## ==========================================================================
## for quadGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We plot this as four individual rectangle gates.
setMethod("glpoints",
          signature(x="quadGate", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE, filterResult=NULL,
                   gpar=flowViz.par.get(), names=FALSE, strict=TRUE, ...)
      {
          if(!is.null(filterResult))
              dropWarn("filterResult", "quadGates", verbose=verbose)
          v <- x@boundary[channels[1]]
          h <- x@boundary[channels[2]]
          mat <- matrix(c(-Inf, v, h, Inf, v, Inf, h, Inf, -Inf, v, -Inf,
                          h, v, Inf, -Inf, h), byrow=TRUE, ncol=4)
          ## we want to be able to use different colors for each population
          col <- rep(gpar$gate$col, 4)
          for(i in 1:4){
              gpar$gate$col <- col[i]
              rg <- rectangleGate(.gate=matrix(mat[i,], ncol=2,
                                               dimnames=list(c("min", "max"),
                                                             channels)))
              glpoints(x=rg, data=data, channels=channels, verbose=FALSE,
                       gpar=gpar, ...)
              addName(x, names, data, gp=gpar$gate.text,...)
          }
      })




## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Only warning that we don't use the filterResult, the rest is pretty
## much the default.
setMethod("glpoints",
          signature(x="polygonGate", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE,
                   filterResult=NULL, gpar=flowViz.par.get(), names=FALSE,
                   ...)
      {
          if(!is.null(filterResult))
              dropWarn("filterResult", "polygonGates", verbose=verbose)
          addLpoints(x=x, data=data, channels=channels, verbose=verbose,
                     filterResult=NULL, gpar=gpar$gate, ...)
          addName(x, names, channels, gp=gpar$gate.text,...)
      })




## ==========================================================================
## for norm2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We don't need to re-evaluate the filter if the filterResult is used,
## so no need for a warning here.
setMethod("glpoints",
          signature(x="norm2Filter", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE, filterResult=NULL,
                   names=FALSE, ...)
      {
          if(is.null(filterResult))
              filterResult <- filter(data, x)
          checkFres(filter=x, fres=filterResult, verbose=verbose)
          fd <- filterDetails(filterResult, identifier(x))
          np <- norm2Polygon(fd, parameters(x))
          identifier(np) <- identifier(x)
          glpoints(x=np, data=data, channels=channels, verbose=FALSE,
                   names=names, ...)
      })







## ==========================================================================
## for kmeansFilters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## This filter produces a multipleFilterResult, so we can't subset directly.
## Instead, we split the original frame and plot each component separately
setMethod("glpoints",
          signature(x="kmeansFilter", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE,
                   filterResult=NULL, gpar=flowViz.par.get(), names=FALSE,
                   ...)
      {
          if(is.null(filterResult))
              filterResult <- filter(data, x)
          multFiltPoints(x=x, data=data, channels=channels, verbose=verbose,
                            filterResult=filterResult, gpar=gpar$gate, ...)
          addName(x, names, split(data[,channels], filterResult),
                  gp=gpar$gate.text,...)
      })




# ==========================================================================
## for subsetFilter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## For now we just plot the top-level filter.
## FIXME: We may want to be able to plot all filters
setMethod("glpoints",
          signature(x="subsetFilter", data="flowFrame", channels="ANY"), 
          function(x, data, channels, ...)
      {
          glpoints(x@filters[[1]], data, channels, ...)
      })






          
