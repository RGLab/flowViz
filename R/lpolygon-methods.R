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
setMethod("glpolygon",
          signature(x="filter", data="missing"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(),
                   strict=TRUE, ...)
      {
          parms <- checkParameterMatch(parameters(x), verbose=verbose,
                                       strict=strict)
          if(!all(is.na(parms)))
              glpolygon(x=x, data=parms, gpar=gpar, verbose=FALSE, ...)
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
                   gpar=flowViz.par.get(), plot=TRUE, names=FALSE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
          xlim <- state("xlim")
          ylim <- state("ylim")
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
              ## add names if necessary
              addName(x, names, data, gpar$gate.text)
          }
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
          res <- ell2Polygon(x, parameters(x))
          glpolygon(res, verbose=verbose, gpar=gpar, plot=plot, ...)
          addName(x, names, data, gpar$gate.text)
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
                   gpar=flowViz.par.get(), names=NULL, ...)
      {
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
          for(i in 1:4){
              gpar$gate$col <- col[i]
              gpar$gate$fill <- fill[i]
              rg <- rectangleGate(.gate=matrix(mat[i,], ncol=2,
                                  dimnames=list(c("min", "max"),
                                  data)))
              res[[i]] <- glpolygon(x=rg, data=data, verbose=FALSE,
                                    gpar=gpar, channels=data, ...)[[1]]
          }
          return(invisible(res))
      })

## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="quadGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(),
                   names=FALSE, ...)
      {
          dropWarn("filterResult", "quadGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar,
                    names=ifelse(names, names(x), FALSE), ...)  
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
          addName(x, names, data=data, gp=gpar$gate.text)
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
                   inverse=FALSE, ...)
      {
          if(!missing(channels))
              data <- channels
          xp <- x@boundaries[,data[1]]
          yp <- x@boundaries[,data[2]]
          class(gpar) <- "gpar"
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
                  xlim <- state("xlim")
                  ylim <- state("ylim")
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
          addName(x, names, data, gpar$gate.text)
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
## for curv1Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter.
setMethod("glpolygon",
          signature(x="curv1Filter", data="ANY"), 
          function(x, data, ...) evalError("curv1Filters"))

## Filter has been evaluated and the filterResult is provided. We plot this
## as a series of rectangleGates.
setMethod("glpolygon",
          signature(x="curv1Filter", data="multipleFilterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(),
                   names=FALSE, ...)
      {
          checkFres(filter=x, fres=data, verbose=verbose)
          fd <- filterDetails(data, identifier(x))
          bounds <- fd$boundaries
          if(all(is.na(bounds[[1]])))
              return(NA)
          lb <- length(bounds)
          ## we want to be able to use different colors for each population
          fill <- rep(gpar$gate$fill, lb)
          col <- rep(gpar$gate$col, lb)
          res <- vector(lb, mode="list")
          n <- if(is.logical(names) && names) names(data)[-1]
          else if(is.character(names)) rep(names, lb) else rep(FALSE, lb)
          for(i in 1:lb){
              tmp <- matrix(bounds[[i]], nrow=2)
              colnames(tmp) <- parameters(x)
              gpar$gate$col <- col[i]
              gpar$gate$fill <- fill[i]
              res[[i]] <- glpolygon(x=rectangleGate(.gate=tmp,
                                    filterId=as.character(n[i])),
                                    verbose=FALSE, gpar=gpar,
                                    names=n[i],
                                    ...)[[1]]
          }
          return(invisible(res))
      })
              
## Evaluate the filter and pass on to the filterResult method.
setMethod("glpolygon",
          signature(x="curv1Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(), ...)
      {
          fres <- filter(data, x)
          glpolygon(x=x, data=fres, verbose=verbose, gpar=gpar, ...)
      })




## ==========================================================================
## for curv2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter.
setMethod("glpolygon",
          signature(x="curv2Filter", data="ANY"), 
          function(x, data, ...)  evalError("curv2Filters")) 
         
## Filter has been evaluated and the filterResult is provided. We plot this
## as a series of polygonGates.
setMethod("glpolygon",
          signature(x="curv2Filter", data="multipleFilterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(),
                   names=FALSE, ...)
      {
          checkFres(filter=x, fres=data, verbose=verbose)
          fd <- filterDetails(data, identifier(x))
          polygons <- fd$polygons
          lf <- length(polygons)
          ## we want to be able to use different colors for each population
          fill <- rep(gpar$gate$fill, lf)
          col <- rep(gpar$gate$col, lf)
          res <- vector(lf, mode="list")
          n <- if(is.logical(names) && names) names(data)[-1]
          else if(is.character(names)) rep(names, lf) else rep(FALSE, lf)
          for(i in 1:lf){
              tmp <- cbind(polygons[[i]]$x, polygons[[i]]$y)
              colnames(tmp) <- parameters(x)
              gpar$gate$col <- col[i]
              gpar$gate$fill <- fill[i]
              res[[i]] <- glpolygon(x=polygonGate(.gate=tmp,
                                                  filterId=as.character(n[i])),
                                    verbose=FALSE, gpar=gpar,
                                    names=n[i], ...)[[1]]
          }
          return(invisible(res))
      })

## Evaluate the filter and pass on to the filterResult method.
setMethod("glpolygon",
          signature(x="curv2Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par.get(), ...)
      {
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
