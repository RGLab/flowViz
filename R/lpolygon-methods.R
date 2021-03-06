## ==========================================================================
## Draw gate regions as polygons on an existing lattice plot. All graphical
## parameters will be passed on to grid.polygon in the end, so
## these are basically extended grid.polygon methods for gates and filters.
## Filters are only evaluated if that is needed for plotting of the
## regions. The methods always dispatch on filters (or filteResults since
## they inherits from filter) and on an optional additional argument
## data, which can be either a character vector, a filterResult or a
## flowFrame. The optional channels argument can be used to pass along
## the plotting parameters, and it takes precendence over everything else.
##
## Most of the methods will just call the next applicable method, the real
## workhorses are signature(x="rectangleGate", data="character") and
## signature(x="polygonGate", data="character") since in the end everything
## boils down to rectangles or polygons.
##
## Population names can be added to the plot if argument 'names' is true.
## Eventually we may want to fix the size and location of the names tags,
## but for now it is always using the gate centroids.



## ==========================================================================
## for all filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We only know how to add the gate boundaries if the definiton of that gate
## contains exactly two dimensions, and even then we need to guess that
## they match the plotted data, hence we warn unless the plotted parameters
## are explicitely provided by the channels argument.
#' Drawing filter regions
#' 
#' These methods extend the lattice \code{\link[lattice:llines]{lpolygon}}
#' methods for drawing of \code{\link[flowCore:filter-class]{filter}} regions.
#' They allow for multiple dispatch, since not all
#' \code{\link[flowCore:filter-class]{filter}} types need to be evaluated for
#' plotting, but this decision should be made internally.
#' 
#' When plotting \code{\link[flowCore:flowFrame-class]{flowFrames}} using the
#' any of the lattice-type \code{plot} method provided by \code{flowViz}, the
#' plotted parameters are recorded, which makes it possible to correctly
#' overlay the outlines of \code{\link[flowCore:filter-class]{filter}} assuming
#' that they are defined for the respective parameters. Warnings and error will
#' be cast for the cases where the parameters are non-distinct or ambigious.
#' These methods are meant to be used within lattice panel functions and are
#' probably not of much use outside of those.
#' 
#' @name glpolygon-methods
#' @aliases glpolygon glpolygon-methods glpolygon,curv1Filter,ANY-method
#' glpolygon,complementFilter,ANY-method glpolygon,curv1Filter,flowFrame-method
#' glpolygon,curv1Filter,missing-method
#' glpolygon,curv1Filter,multipleFilterResult-method
#' glpolygon,curv2Filter,ANY-method glpolygon,curv2Filter,flowFrame-method
#' glpolygon,curv2Filter,multipleFilterResult-method
#' glpolygon,filter,missing-method glpolygon,filterResult,flowFrame-method
#' glpolygon,filterResult,missing-method glpolygon,filterResult,ANY-method
#' glpolygon,kmeansFilter,ANY-method glpolygon,norm2Filter,ANY-method
#' glpolygon,norm2Filter,flowFrame-method
#' glpolygon,norm2Filter,logicalFilterResult-method
#' glpolygon,polygonGate,character-method
#' glpolygon,polygonGate,filterResult-method
#' glpolygon,polygonGate,flowFrame-method glpolygon,quadGate,character-method
#' glpolygon,quadGate,filterResult-method glpolygon,quadGate,flowFrame-method
#' glpolygon,rectangleGate,character-method
#' glpolygon,rectangleGate,filterResult-method
#' glpolygon,rectangleGate,flowFrame-method
#' glpolygon,ellipsoidGate,character-method
#' glpolygon,ellipsoidGate,filterResult-method
#' glpolygon,ellipsoidGate,flowFrame-method glpolygon,subsetFilter,ANY-method
#' @docType methods
#' @return The methods will return the outlines of the gate region as polygon
#' vertices.
#' @section Methods:
#' 
#' \describe{
#' 
#' \item{x = "filter", data = "missing"}{ General method for all objects
#' inheriting from \code{\link[flowCore:filter-class]{filter}}. This is used as
#' the default when no more explicit method is found. It tries to find the
#' plotted parameters from the internal \code{flowViz.state} environment. This
#' only works if the flow data has been plotted using the \code{plot} methods
#' provided by this \code{flowViz} package. }
#' 
#' \item{x = "filterResult", data = "missing"}{ General method for all
#' \code{\link[flowCore:filterResult-class]{filterResult}} object. This
#' basically extracts the \code{\link[flowCore:filter-class]{filter}} from the
#' \code{\link[flowCore:filterResult-class]{filterResult}} and dispatches on
#' that. }
#' 
#' \item{x = "filterResult", data = "flowFrame"}{ For some
#' \code{\link[flowCore:filter-class]{filter}} types we need the raw data to
#' re-evaluate the filter. }
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
#' \code{\link[flowStats:curv2Filter-class]{curv2Filter}}.}
#' 
#' \item{x = "curv2Filter", data = "flowFrame"}{ see above }
#' 
#' \item{x = "curv2Filter", data = "multipleFilterResult"}{ see above }
#' 
#' \item{x = "kmeansFilter", data = "ANY"}{ We don't know how to plot regions
#' of a \code{\link[flowCore:kmeansFilter]{kmeansFilter}}, hence we warn. }
#' 
#' \item{x = "norm2Filter", data = "ANY"}{ We either need a
#' \code{\link[flowCore:filterResult-class]{filterResult}} or the raw data as a
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} for
#' \code{\link[flowStats:norm2Filter-class]{norm2Filter}}.}
#' 
#' \item{x = "norm2Filter", data = "flowFrame"}{ see above }
#' 
#' \item{x = "norm2Filter", data = "logicalFilterResult"}{ see above }
#' 
#' \item{x = "polygonGate", data = "character"}{ We can plot a
#' \code{\link[flowCore]{polygonGate}} directly from the gate definition. }
#' 
#' \item{x = "polygonGate", data = "filterResult"}{ see above }
#' 
#' \item{x = "polygonGate", data = "flowFrame"}{ see above }
#' 
#' \item{x = "quadGate", data = "character"}{ We can plot a
#' \code{\link[flowCore]{quadGate}} directly from the gate definition. }
#' 
#' \item{x = "quadGate", data = "filterResult"}{ see above }
#' 
#' \item{x = "quadGate", data = "flowFrame"}{ see above }
#' 
#' \item{x = "rectangleGate", data = "character"}{ We can plot a
#' \code{\link[flowCore]{rectangleGate}} directly from the gate definition. }
#' 
#' \item{x = "rectangleGate", data = "filterResult"}{ see above }
#' 
#' \item{x = "rectangleGate", data = "flowFrame"}{ see above }
#' 
#' \item{x = "ellipsoidGate", data = "character"}{ We can plot a
#' \code{\link[flowCore]{rectangleGate}} directly from the gate definition. }
#' 
#' \item{x = "ellipsoidGate", data = "filterResult"}{ see above }
#' 
#' \item{x = "ellipsoidGate", data = "flowFrame"}{ see above }
#' 
#' }
#' @author F. Hahne
#' @seealso
#' 
#' \code{\link[flowCore:filter-class]{filter}},
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowViz:glpoints-methods]{glpoints}}
#' @keywords methods
#' @export 
#' @param x filter or filterResult or any derived filter class
#' @param data flowFrame or filterResult or character or missing or ANY
#' @param verbose logical
#' @param gpar a list of graphical parameters. see 'help(flowViz.par.get)' for details.
#' @param strict logical
#' @param ... other arguments
setMethod("glpolygon",
          signature(x="filter", data="missing"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(),
                   strict=TRUE, ...)
      {
          parms <- checkParameterMatch(parameters(x), verbose=verbose,
                                       strict=strict,...)
          if(!all(is.na(parms)))
              glpolygon(x=x, data=parms, gpar=gpar, verbose=FALSE,strict=strict, ...)
      })


## Extract the filter definiton from a filterResult and pass that on
## along with it. We decide later if we really need it or not.
setMethod("glpolygon",
          signature(x="filterResult", data="ANY"), 
          definition=function(x, data, verbose=TRUE,
          gpar=flowViz.par.get(), ...)
      {
          fd <- filterDetails(x)
          filt <- fd[[length(fd)]]$filter
          if(!missing(data) && is.character(data) &&
             ! ("channels" %in% names(list(...))))
              glpolygon(filt, x, verbose=FALSE, channels=data,
                        gpar=gpar, ...)
          else
              glpolygon(filt, x, verbose=FALSE, gpar=gpar, ...) 
      })


## We don't need the flowFrame if the filter is already evaluated, but we
## can check that the IDs are matching before dropping it.
setMethod("glpolygon",
          signature(x="filterResult", data="flowFrame"), 
          definition=function(x, data, verbose=TRUE,
          gpar=flowViz.par.get(), ...){
              checkIdMatch(x=x, f=data)
              dropWarn("flowFrame", "filterResults", verbose=verbose)
              glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
          })




## ==========================================================================
## for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector. If channels is
## provided this takes precedence.
setMethod("glpolygon",
          signature(x="rectangleGate", data="character"), 
          function(x, data, channels, verbose=TRUE,
                   gpar=flowViz.par.get(), plot=TRUE, names=FALSE,xlim,ylim, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
#          browser()
          if(plot){
              ## 1D rectangular gate (region gate). We cheat those by drawing a
              ## rectangle that is slightly bigger than the drawing limits in the
              ## missing channel and let clipping work for us.
              if(length(parms)==1){
                  mt <- match(parms, data)
                  if(mt==1){
                      glrect(fixInf(x@min, xlim), ylim[1]-diff(ylim)*10,
                             fixInf(x@max, xlim), ylim[2]+diff(ylim)*10,
                             gp=gpar$gate, ...)
                  }   
                  else if(mt==2){
                      glrect(xlim[1]-diff(xlim)*10, fixInf(x@min, ylim),
                             xlim[2]+diff(xlim)*10, fixInf(x@max, ylim),
                             gp=gpar$gate, ...)
                  }else
                  stop("How did you end up here????")
              }else{## 2D rectangular gate   
                  bl <- x@min[data]
                  tr <- x@max[data]
                  glrect(fixInf(bl[1], xlim), fixInf(bl[2], ylim),
                         fixInf(tr[1], xlim),  fixInf(tr[2], ylim),
                         gp=gpar$gate, ...)
              }
              
          }
		  ## add names if necessary
		  addName(x, names, data, gp = gpar$gate.text,xlim=xlim,ylim=ylim,...)
          res <- rbind(x@min[data], x@max[data])
          res <- res[,!apply(res, 2, function(z) all(is.na(z))), drop=FALSE]
          return(invisible(list(res)))
      })

## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="rectangleGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(), ...)
      {
          dropWarn("filterResult", "rectangleGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })

## We can drop the flowFrame, don't need it for rectangleGates.
setMethod("glpolygon",
          signature(x="rectangleGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(), ...)
      {
          dropWarn("flowFrame", "rectangleGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })




## ==========================================================================
## for ellipsoidGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("glpolygon",
          signature(x="ellipsoidGate", data="character"), 
          function(x, data, channels, verbose=TRUE, gpar=flowViz.par.get(),
                   plot=TRUE, names=FALSE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
	  ## We coerce to a polygon gate and plot that
          res <- as(x, "polygonGate")
          glpolygon(res, verbose=verbose, gpar=gpar, plot=plot, ...)
          addName(x, names, data, gp = gpar$gate.text,...)
          return(invisible(res@boundaries))
      })


## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="ellipsoidGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(), ...)
      {
          dropWarn("filterResult", "ellipsoidGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })

## We can drop the flowFrame, don't need it for rectangleGates.
setMethod("glpolygon",
          signature(x="ellipsoidGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(), ...)
      {
          dropWarn("flowFrame", "ellipsoidGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })






## ==========================================================================
## for quadGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We plot this as four individual rectangle gates. If channels is
## provided this takes precedence.
setMethod("glpolygon",
          signature(x="quadGate", data="character"), 
          function(x, data, channels, verbose=TRUE,
                   gpar=flowViz.par.get(), names=NULL,abs=FALSE,pos=NULL,...)
      {
#		  browser()
          if(!missing(channels))
              data <- channels
          v <- x@boundary[data[1]]
          h <- x@boundary[data[2]]
          mat <- matrix(c(-Inf, v, h, Inf, v, Inf, h, Inf, -Inf, v, -Inf,
                          h, v, Inf, -Inf, h), byrow=TRUE, ncol=4)
          ## we want to be able to use different colors for each population
          fill <- rep(gpar$gate$fill, 4)
          col <- rep(gpar$gate$col, 4)
          col <- col[c(2,1,3,4)]
          fill <- fill[c(2,1,3,4)]
          res <- vector(4, mode="list")
		  #quadrant iteration in this loop goes by :x-y+,x+y+,x-y-,x+y-
		  #	yet populations generated from quadGate is:++,-+,+-,--
		  #so we need to re-order the names character
		  if(length(names)==1)
			  names<-rep(names,4)
		  names<-names[c(2,1,4,3)]
		  # fix the location of gate labels at four corners for quadGate
		  
		  poslist<-vector(mode="list",4)
		  poslist[[1]]<-c(0.15,0.9)
		  poslist[[2]]<-c(0.9,0.9)
		  poslist[[3]]<-c(0.15,0.1)
		  poslist[[4]]<-c(0.9,0.1)
          for(i in 1:4){
              gpar$gate$col <- col[i]
              gpar$gate$fill <- fill[i]
              rg <- rectangleGate(.gate=matrix(mat[i,], ncol=2,
                                  dimnames=list(c("min", "max"),
                                  data)))
		  
              res[[i]] <- glpolygon(x=rg, data=data, verbose=FALSE,
                                    gpar=gpar, channels=data,names=names[i],abs=abs,pos=poslist[[i]],...)[[1]]
										
          }
          return(invisible(res))
      })

## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="quadGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(),
                   names=FALSE, ...)
      {
#		  browser()
          dropWarn("filterResult", "quadGates", verbose=verbose)
#          glpolygon(x=x, verbose=verbose, gpar=gpar,
#                    names=ifelse(names, names(x), FALSE), ...)  

			glpolygon(x=x, verbose=verbose, gpar=gpar,
                    names=names, ...)
      })

## We can drop the dataFrame, don't need it for rectangleGates.
setMethod("glpolygon",
          signature(x="quadGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(),
                   names=FALSE, ...)
      {
          dropWarn("flowFrame", "quadGates", verbose=verbose)
          res <- glpolygon(x=x, verbose=verbose, gpar=gpar,
                    names=FALSE, ...)
          addName(x, names, data=data, gp = gpar$gate.text,...)
          return(invisible(res))
      })




## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector. If channels is
## provided this takes precedence.
setMethod("glpolygon",
          signature(x="polygonGate", data="character"), 
          function(x, data, channels, verbose=TRUE,
                   gpar=flowViz.par.get(), plot=TRUE, names=FALSE,
                   inverse=FALSE,xlim,ylim, ...)
      {
		  
          if(!missing(channels))
              data <- channels
          xp <- x@boundaries[,data[1]]
          yp <- x@boundaries[,data[2]]
#          class(gpar) <- "gpar"
          if(plot){
              if(inverse){
                  mx <- which.min(yp)
                  xpi <- rep(xp, 2)
                  ypi <- rep(yp, 2)
                  nextDiff <- min(which(xpi[-(1:mx)] != xpi[mx]))+mx
                  if(xpi[nextDiff] < xp[mx]){
                      xp <- rev(xp)
                      yp <- rev(yp)
                      mx <- which.min(yp)
                  }
                  
                  xoff <- diff(xlim)
                  yoff <- diff(ylim)
                  boundingBox <- cbind(c(rep(xlim-xoff, 2),
                                         rep(xlim+xoff, 2)),
                                       c(ylim-yoff, ylim+yoff,
                                         ylim+yoff, ylim-yoff))
                  xpi <- c(rep(boundingBox[1,1], 2), xp, xp[1],
                          boundingBox[,1], boundingBox[1,1])
                  ypi <- c(boundingBox[1,2], yp[mx], yp, yp[1],
                          yp[mx], boundingBox[-1,2],
                             boundingBox[1,2])
                  gpari <- gpar$gate
                  gpari$col <- "transparent"
                  glpoly(xpi, ypi, gp=gpari)
                  gpari$col <- gpar$gate$col
                  gpari$fill <- "transparent"
                  glpoly(xp, yp, gp=gpari)
              }else{
                  glpoly(xp, yp, gp=gpar$gate)
              }
          }
          addName(x, names, data, gp = gpar$gate.text,xlim=xlim,ylim=ylim,...)
          res <- cbind(xp,yp)
          res <- res[,!apply(res, 2, function(z) all(is.na(z))), drop=FALSE]
          return(invisible(list(res)))
      })

## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="polygonGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(), ...)
      {
          dropWarn("filterResult", "polygonGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })

## we can drop the flowFrame, don't need it for polygonGates.
setMethod("glpolygon",
          signature(x="polygonGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(), ...)
      {
          dropWarn("flowFrame", "polygonGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })




# ==========================================================================
## for norm2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter.
setMethod("glpolygon",
          signature(x="norm2Filter", data="ANY"), 
          function(x, data, ...) evalError("norm2Filters"))

## Filter has been evaluated and the filterResult is provided. We convert
## this into a polygon gate and plot that.
setMethod("glpolygon",
          signature(x="norm2Filter", data="logicalFilterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(),
                   names=FALSE, ...){
              checkFres(filter=x, fres=data, verbose=verbose)
              fd <- filterDetails(data, identifier(x))
              np <- norm2Polygon(fd, parameters(x))
              identifier(np) <- identifier(x)
              glpolygon(x=np, verbose=FALSE, gpar=gpar, names=names, ...)
          })

## Evaluate the filter and pass on to the filterResult method.
setMethod("glpolygon",
          signature(x="norm2Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(), ...){
              fres <- filter(data, x)
               glpolygon(x=x, data=fres, verbose=verbose, gpar=gpar, ...)
          })





## ==========================================================================
## for kmeansFilters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We don't know how to plot these, hence we warn
setMethod("glpolygon",
          signature(x="kmeansFilter", data="ANY"), 
          function(x, data, verbose=TRUE, ...)
          nnWarn("kmeansFilter", verbose))



## ==========================================================================
## for complementFilter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We have to do this on the level of polygons or rectangles, so here we
## only set a flag and pass on the basic filter.
## FIXME: This doesn't work for multi-polygon filters so far.
setMethod("glpolygon",
          signature(x="complementFilter", data="ANY"), 
          function(x, data, verbose=TRUE, inverse=FALSE, ...)
      {
          glpolygon(x@filters[[1]], data, verbose=verbose,
                    inverse=!inverse, ...)
      })




# ==========================================================================
## for subsetFilter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## For now we just plot the top-level filter.
## FIXME: We may want to be able to plot all filters
setMethod("glpolygon",
          signature(x="subsetFilter", data="ANY"), 
          function(x, data, verbose=TRUE, inverse=FALSE, ...)
      {
          glpolygon(x@filters[[1]], data, verbose=verbose, ...)
      })
