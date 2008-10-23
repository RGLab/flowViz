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



## ==========================================================================
## for all filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We only know how to add points within the gate if the definiton of that
## gate fits the parameters of the data provided.
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
          addName(x, names, channels, gpar$gate.text)
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
          addName(x, names, channels, gpar$gate.text)
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
              addName(x, names, data, gp=gpar$gate.text)
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
          addName(x, names, channels, gp=gpar$gate.text)
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
## for curv1Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## This filter produces a multipleFilterResult, so we can't subset directly.
## Instead, we split the original frame and plot each component separately
setMethod("glpoints",
          signature(x="curv1Filter", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE,
                   filterResult=NULL, gpar=flowViz.par.get(), names=FALSE,
                   ...)  
      {
          if(is.null(filterResult))
              filterResult <- filter(data, x)
          multFiltPoints(x=x, data=data, channels=channels, verbose=verbose,
                         filterResult=filterResult, gpar=gpar$gate, ...)
          addName(x, name=names, data=filterResult, gp=gpar$gate.text)
      })




## ==========================================================================
## for curv2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## This filter produces a multipleFilterResult, so we can't subset directly.
## Instead, we split the original frame and plot each component separately.
setMethod("glpoints",
          signature(x="curv2Filter", data="flowFrame", channels="character"), 
          function(x, data, channels, verbose=TRUE,
                   filterResult=NULL, gpar=flowViz.par.get(), names=FALSE,
                   ...)
      {
          if(is.null(filterResult))
              filterResult <- filter(data, x)
          multFiltPoints(x=x, data=data, channels=channels, verbose=verbose,
                         filterResult=filterResult, gpar=gpar$gate, ...)
          addName(x, name=names, data=filterResult, gp=gpar$gate.text)
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
                  gp=gpar$gate.text)
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






          
