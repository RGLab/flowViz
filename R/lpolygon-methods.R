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



## ==========================================================================
## for all filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We only know how to add the gate boundaries if the definiton of that gate
## contains exactly two dimensions, and even then we need to guess that
## they match the plotted data, hence we warn unless the plotted parameters
## are explicitely provided by the channels argument.
setMethod("glpolygon",
          signature(x="filter", data="missing"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
      {
          parms <- checkParameterMatch(parameters(x), verbose=verbose)
          glpolygon(x=x, data=parms, gpar=gpar, verbose=FALSE, ...)
      })


## Extract the filter definiton from a filterResult and pass that on
## along with it. We decide later if we really need it or not.
setMethod("glpolygon",
          signature(x="filterResult", data="ANY"), 
          definition=function(x, data, verbose=TRUE,
          gpar=flowViz.par(), ...)
      {
          filt <- filterDetails(x)$filter
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
          gpar=flowViz.par(), ...){
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
                   gpar=flowViz.par(), plot=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
          usr <- rep(c(-1e4, 1e4), 2)
          if(plot){
              if(length(parms)==1){
                  mt <- match(parms, data)## 1D rectangular gate (region gate)
                  if(mt==1)
                      glrect(fixInf(x@min, usr[1:2]), usr[3],
                             fixInf(x@max, usr[1:2]), usr[4], gp=gpar, ...)
                  else if(mt==2)
                      glrect(usr[1], fixInf(x@min, usr[3:4]), usr[2],
                             fixInf(x@max, usr[3:4]), gp=gpar, ...)
                  else
                      stop("How did you end up here????")
              }else{## 2D rectangular gate   
                  bl <- x@min[data]
                  tr <- x@max[data]
                  glrect(fixInf(bl[1], usr[1:2]), fixInf(bl[2], usr[3:4]),
                         fixInf(tr[1], usr[1:2]),  fixInf(tr[2], usr[3:4]),
                         gp=gpar, ...)
              }
          }
          res <- rbind(x@min[data], x@max[data])
          res <- res[,!apply(res, 2, function(z) all(is.na(z))), drop=FALSE]
          return(invisible(list(res)))
      })

## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="rectangleGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
      {
          dropWarn("filterResult", "rectangleGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })

## We can drop the flowFrame, don't need it for rectangleGates.
setMethod("glpolygon",
          signature(x="rectangleGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
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
          function(x, data, channels, verbose=TRUE, gpar=flowViz.par(),
                   plot=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
	  ## We coerce to a polygon gate and plot that
          res <- ell2Polygon(fd, parameters(x))
          glpolygon(res, verbose=verbose, gpar=gpar, plot=plot, ...)
          return(invisible(res@boundaries))
      })


## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="ellipsoidGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
      {
          dropWarn("filterResult", "ellipsoidGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })

## We can drop the flowFrame, don't need it for rectangleGates.
setMethod("glpolygon",
          signature(x="ellipsoidGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
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
                   gpar=flowViz.par(), ...)
      {
          if(!missing(channels))
              data <- channels
          v <- x@boundary[data[1]]
          h <- x@boundary[data[2]]
          mat <- matrix(c(-Inf, h, v, Inf, h, Inf, v, Inf, -Inf, h, -Inf,
                          v, h, Inf, -Inf, v), byrow=TRUE, ncol=4)
          ## we want to be able to use different colors for each population
          fill <- rep(gpar$fill, 4)
          col <- rep(gpar$col, 4)
          res <- vector(4, mode="list")
          for(i in 1:4){
              gpar$col <- col[i]
              gpar$fill <- fill[i]
              rg <- rectangleGate(.gate=matrix(mat[i,], ncol=2,
                                  dimnames=list(c("min", "max"),
                                  parameters(x))))
              res[[i]] <- glpolygon(x=rg, data=data, verbose=FALSE,
                                    gpar=gpar, channels=data, ...)[[1]]
          }
          return(invisible(res))
      })

## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="quadGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
      {
          dropWarn("filterResult", "quadGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })

## We can drop the dataFrame, don't need it for rectangleGates.
setMethod("glpolygon",
          signature(x="quadGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
      {
          dropWarn("flowFrame", "quadGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })




## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector. If channels is
## provided this takes precedence.
setMethod("glpolygon",
          signature(x="polygonGate", data="character"), 
          function(x, data, channels, verbose=TRUE,
                   gpar=flowViz.par(), plot=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          xp <- x@boundaries[,data[1]]
          yp <- x@boundaries[,data[2]]
          class(gpar) <- "gpar"
          if(plot)
              grid.polygon(xp, yp, default.units = "native", gp=gpar)
          res <- cbind(xp,yp)
          res <- res[,!apply(res, 2, function(z) all(is.na(z))), drop=FALSE]
          return(invisible(list(res)))
      })

## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="polygonGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
      {
          dropWarn("filterResult", "polygonGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar, ...)
      })

## we can drop the flowFrame, don't need it for polygonGates.
setMethod("glpolygon",
          signature(x="polygonGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
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
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...){
              checkFres(filter=x, fres=data, verbose=verbose)
              fd <- filterDetails(data, identifier(x))
              glpolygon(x=norm2Polygon(fd, parameters(x)),
                       verbose=FALSE, gpar=gpar, ...)
          })

## Evaluate the filter and pass on to the filterResult method.
setMethod("glpolygon",
          signature(x="norm2Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...){
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
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
      {
          checkFres(filter=x, fres=data, verbose=verbose)
          fd <- filterDetails(data, identifier(x))
          bounds <- fd$boundaries
          if(all(is.na(bounds[[1]])))
              return(NA)
          lb <- length(bounds)
          ## we want to be able to use different colors for each population
          fill <- rep(gpar$fill, lb)
          col <- rep(gpar$col, lb)
          res <- vector(lb, mode="list")
          for(i in 1:lb){
              tmp <- matrix(bounds[[i]], nrow=2)
              colnames(tmp) <- parameters(x)
              gpar$col <- col[i]
              gpar$fill <- fill[i]
              res[[i]] <- glpolygon(x=rectangleGate(.gate=tmp),
                                    verbose=FALSE, gpar=gpar, ...)[[1]]
          }
          return(invisible(res))
      })
              
## Evaluate the filter and pass on to the filterResult method.
setMethod("glpolygon",
          signature(x="curv1Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
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
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
      {
          checkFres(filter=x, fres=data, verbose=verbose)
          fd <- filterDetails(data, identifier(x))
          polygons <- fd$polygons
          lf <- length(polygons)
          ## we want to be able to use different colors for each population
          fill <- rep(gpar$fill, lf)
          col <- rep(gpar$col, lf)
          res <- vector(lf, mode="list")
          for(i in 1:lf){
              tmp <- cbind(polygons[[i]]$x, polygons[[i]]$y)
              colnames(tmp) <- parameters(x)
              gpar$col <- col[i]
              gpar$fill <- fill[i]
              res[[i]] <- glpolygon(x=polygonGate(boundaries=tmp),
                                     verbose=FALSE, gpar=gpar, ...)[[1]]
          }
          return(invisible(res))
      })

## Evaluate the filter and pass on to the filterResult method.
setMethod("glpolygon",
          signature(x="curv2Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(), ...)
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
