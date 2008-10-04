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
                   gpar=flowViz.par(), plot=TRUE, names=FALSE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
          xlim <- state("xlim")
          ylim <- state("ylim")
          if(plot){
              if(!is.character(names))
                  n <- paste(identifier(x), "+", sep="")
              else{
                  n <- names
                  names <- TRUE
              }
              if(length(parms)==1){
                  mt <- match(parms, data)## 1D rectangular gate (region gate)
                  if(mt==1){
                      glrect(fixInf(x@min, xlim*2), ylim[1]*2,
                             fixInf(x@max, xlim*2), ylim[2]*2, gp=gpar, ...)
                      if(names){
                          xx <- c(fixInf(x@min, xlim), fixInf(x@max, xlim))
                          panel.text(mean(xx), ylim[2], n, pos=1,
                                     col="#00000040", cex=1.2)
                      }   
                  }   
                  else if(mt==2){
                      glrect(xlim[1]*2, fixInf(x@min, ylim*2), xlim[2]*2,
                             fixInf(x@max, ylim*2), gp=gpar, ...)
                      if(names){
                          yy <- c(fixInf(x@min, ylim), fixInf(x@max, ylim))
                          panel.text(xlim[1], mean(yy), n, pos=4,
                                     col="#00000040", cex=1.2)
                      }   
                  }else
                      stop("How did you end up here????")
              }else{## 2D rectangular gate   
                  bl <- x@min[data]
                  tr <- x@max[data]
                  glrect(fixInf(bl[1], xlim*2), fixInf(bl[2], ylim*2),
                         fixInf(tr[1], xlim*2),  fixInf(tr[2], ylim*2),
                         gp=gpar, ...)
                  if(names){
                       xx <- c(fixInf(bl[1], xlim), fixInf(tr[1], xlim))
                       yy <- c(fixInf(bl[2], ylim), fixInf(tr[2], ylim))
                       panel.text(mean(xx), mean(yy), n, col="#00000040",
                                  adj=0.5, cex=1.2)
                   }
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
                   plot=TRUE, names=FALSE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
	  ## We coerce to a polygon gate and plot that
          res <- ell2Polygon(x, parameters(x))
          glpolygon(res, verbose=verbose, gpar=gpar, plot=plot, ...)
          n <- paste(identifier(x), "+", sep="")
          if(names)
              panel.text(x@mean[1], x@mean[2], n, adj=0.5,
                         col="#00000040", cex=1.2)
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
                   gpar=flowViz.par(), names=NULL, ...)
      {
          if(!missing(channels))
              data <- channels
          if(!is.null(names) && names==TRUE)
              stop("Either the data or a filterResult have to be provided ",
                   "to print names", call.=FALSE)
          v <- x@boundary[data[1]]
          h <- x@boundary[data[2]]
          mat <- matrix(c(-Inf, v, h, Inf, v, Inf, h, Inf, -Inf, v, -Inf,
                          h, v, Inf, -Inf, h), byrow=TRUE, ncol=4)
          ## we want to be able to use different colors for each population
          fill <- rep(gpar$fill, 4)
          col <- rep(gpar$col, 4)
          res <- vector(4, mode="list")
          for(i in 1:4){
              gpar$col <- col[i]
              gpar$fill <- fill[i]
              rg <- rectangleGate(.gate=matrix(mat[i,], ncol=2,
                                  dimnames=list(c("min", "max"),
                                  data)))
              res[[i]] <- glpolygon(x=rg, data=data, verbose=FALSE,
                                    gpar=gpar, channels=data, ...)[[1]]
          }
          if(is.character(names)){
              xlim <- state("xlim")
              ylim <- state("ylim")
              yoff <- diff(ylim)/50
              trans <- match(colnames(names)[1], data)-1
              if(trans)
                  names <- matrix(c(names[2,2], names[2,1], names[1,2],
                                    names[1,1]), ncol=2)
              panel.text(c(mean(c(xlim[1], v)), mean(c(v, xlim[2]))),
                         rep(ylim[2], 2)-yoff, names[1,], adj=c(0.5, 1),
                         col="#00000040", cex=1.2)
              panel.text(c(mean(c(xlim[1], v)), mean(c(v, xlim[2]))),
                         rep(ylim[1], 2)+yoff, names[2,], adj=c(0.5, 0),
                         col="#00000040", cex=1.2)
          }
              
          return(invisible(res))
      })

## We can ignore the filterResult, don't need it to plot the gate.
setMethod("glpolygon",
          signature(x="quadGate", data="filterResult"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(),
                   names=FALSE, ...)
      {
          dropWarn("filterResult", "quadGates", verbose=verbose)
          glpolygon(x=x, verbose=verbose, gpar=gpar,
                    names=ifelse(names, names(x), NULL), ...)  
      })

## We can drop the dataFrame, don't need it for rectangleGates.
setMethod("glpolygon",
          signature(x="quadGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, gpar=flowViz.par(),
                   names=FALSE, ...)
      {
          dropWarn("flowFrame", "quadGates", verbose=verbose)
          names <- if(names){
              desc <- flowCore:::popNames(data, x)
              n <- matrix(c(sprintf("%s-%s+", desc[1], desc[2]),
                            sprintf("%s-%s-", desc[1], desc[2]),
                            sprintf("%s+%s+", desc[1], desc[2]),
                            sprintf("%s+%s-", desc[1], desc[2])),
                          ncol=2)
              colnames(n) <- names(desc)
              n
          } else NULL
              glpolygon(x=x, verbose=verbose, gpar=gpar,
                        names=names, ...)
      })




## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector. If channels is
## provided this takes precedence.
setMethod("glpolygon",
          signature(x="polygonGate", data="character"), 
          function(x, data, channels, verbose=TRUE,
                   gpar=flowViz.par(), plot=TRUE, names=FALSE,
                   inverse=FALSE, ...)
      {
          if(!missing(channels))
              data <- channels
          xp <- x@boundaries[,data[1]]
          yp <- x@boundaries[,data[2]]
          class(gpar) <- "gpar"
          if(plot){
              if(!is.character(names))
                  n <- paste(identifier(x), "+", sep="")
              else{
                  n <- names
                  names <- TRUE
              }
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
                  #browser()
                  gpari <- gpar
                  gpari$col <- "transparent"
                  grid.polygon(xpi, ypi, default.units = "native", gp=gpari)
                  gpari$col <- gpar$col
                  gpari$fill <- "transparent"
                  grid.polygon(xp, yp, default.units = "native", gp=gpari)
              }else{
                  grid.polygon(xp, yp, default.units = "native", gp=gpar)
              }
              if(names)
                  panel.text(mean(xp), mean(yp), n, adj=c(0.5, 0.5),
                             col="#00000040", cex=1.2)
          }
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
          function(x, data, verbose=TRUE, gpar=flowViz.par(),
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
          function(x, data, verbose=TRUE, gpar=flowViz.par(),
                   names=FALSE, ...)
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
          n <- names(data)[-1]
          for(i in 1:lb){
              tmp <- matrix(bounds[[i]], nrow=2)
              colnames(tmp) <- parameters(x)
              gpar$col <- col[i]
              gpar$fill <- fill[i]
              res[[i]] <- glpolygon(x=rectangleGate(.gate=tmp,
                                    filterId=n[i]),
                                    verbose=FALSE, gpar=gpar,
                                    names=ifelse(names, n[i], FALSE),
                                    ...)[[1]]
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
          function(x, data, verbose=TRUE, gpar=flowViz.par(),
                   names=FALSE, ...)
      {
          checkFres(filter=x, fres=data, verbose=verbose)
          fd <- filterDetails(data, identifier(x))
          polygons <- fd$polygons
          lf <- length(polygons)
          ## we want to be able to use different colors for each population
          fill <- rep(gpar$fill, lf)
          col <- rep(gpar$col, lf)
          res <- vector(lf, mode="list")
          n <- names(data)[-1]
          for(i in 1:lf){
              tmp <- cbind(polygons[[i]]$x, polygons[[i]]$y)
              colnames(tmp) <- parameters(x)
              gpar$col <- col[i]
              gpar$fill <- fill[i]
              res[[i]] <- glpolygon(x=polygonGate(boundaries=tmp,
                                     filterId=n[i]),
                                     verbose=FALSE, gpar=gpar,
                                    names=ifelse(names, n[i], FALSE),
                                    ...)[[1]]
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



## ==========================================================================
## for complementFilter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We have to do this on the level of polygons or rectangles, so here we
## only set a flag and pass on the basic filter.
setMethod("glpolygon",
          signature(x="complementFilter", data="ANY"), 
          function(x, data, verbose=TRUE, inverse=FALSE, ...)
      {
          glpolygon(x@filters[[1]], data, verbose=verbose,
                    inverse=!inverse, ...)
      })
