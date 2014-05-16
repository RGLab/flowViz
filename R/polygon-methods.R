## ==========================================================================
## Draw gate regions as polygons on an existing plot. All graphical
## parameters will be passed on to graphics::polygon in the end, so
## these are basically extended polygon methods for gates and filters.
## Filters are only evaluated if that is needed for plotting of the
## regions.
## FIXME: Need to figure out how to make polygon a generic without masking
##        the default function



## ==========================================================================
## allow for dispatch to default polygon function in graphics
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#setMethod("polygon", signature(x="ANY"),
#          definition=function(x, ...) graphics::polygon(x, ...))




## ==========================================================================
## for all filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We only know how to add the gate boundaries if the definiton of that gate
## contains exactly two dimensions, and even then we need to guess that
## they match the plotted data, hence we warn. Supplying the names of the
## plotted channels via the "channels" argument always gets precendence.
## All downstream methods will eventually run across this method, so we only
## need to check once. 
setMethod("gpolygon",
          signature(x="filter", data="missing"), 
          function(x, data, verbose=TRUE, ...)
      {
          parms <- checkParameterMatch(parameters(x), verbose=verbose, ...)
          gpolygon(x, parms, ...)
      })


## Extract the filter definiton from a filterResult and pass that on
## along with it
setMethod("gpolygon",
          signature(x="filterResult", data="ANY"), 
          function(x, data, verbose=TRUE, ...)
      {
          fd <- filterDetails(x)
          filt <- fd[[length(fd)]]$filter
          if(!missing(data) && is.character(data) &&
             ! ("channels" %in% names(list(...))))
              gpolygon(filt, x, verbose=FALSE, channels=data, ...)
          else
             gpolygon(filt, x, verbose=FALSE, ...) 
      })

## We don't need the flowFrame if the filter is already evaluated, but we
## can check that the IDs are matching
setMethod("gpolygon",
          signature(x="filterResult", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...){
              checkIdMatch(x=x, f=data)
              gpolygon(x, verbose=verbose, ...)
              dropWarn("flowFrame", "filterResults", verbose=verbose)
          })




## ==========================================================================
## for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("gpolygon",
          signature(x="rectangleGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
          usr <- par("usr")
          if(length(parms)==1){## 1D rectangular gate (region gate)
              mt <- match(parms, data)
              if(mt==1)
                  rect(fixInf(x@min, usr[1]), usr[3], fixInf(x@max, usr[2]),
                       usr[4], ...)
              else if(mt==2)
                  rect(usr[1], fixInf(x@min, usr[3]), usr[2],
                       fixInf(x@max, usr[4]), ...)
              else
                  stop("How did you end up here????")
          }else{## 2D rectangular gate    
              bl <- x@min[data]
              tr <- x@max[data]
              rect(fixInf(bl[1], usr[1:2]), fixInf(bl[2], usr[3:4]),
                   fixInf(tr[1], usr[1:2]),  fixInf(tr[2], usr[3:4]), ...)
          }
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("gpolygon",
          signature(x="rectangleGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult", "rectangleGates", verbose=verbose)
          gpolygon(x, verbose=verbose, ...)
      })

## we can drop the dataFrame, don't need it for rectangleGates
setMethod("gpolygon",
          signature(x="rectangleGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("flowFrame", "rectangleGates", verbose=verbose)
          gpolygon(x, verbose=verbose, ...)
      })



## ==========================================================================
## for ellipsoidGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("gpolygon",
          signature(x="ellipsoidGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
	  ## We coerce to a polygon gate and plot that
          gpolygon(ell2Polygon(fd, parameters(x)), verbose=verbose, ...)      
      })


## We can ignore the filterResult, don't need it to plot the gate
setMethod("gpolygon",
          signature(x="ellipsoidGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult",
                   "ellipsoidGates", verbose=verbose)
          gpolygon(x, verbose=verbose, ...)
      })


## we can drop the dataFrame, don't need it for ellipsoidGates
setMethod("gpolygon",
          signature(x="ellipsoidGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("flowFrame", "ellipsoidGates", verbose=verbose)
          gpolygon(x, verbose=verbose, ...)
      })







## ==========================================================================
## for quadGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We plot this as four individual rectangle gates
setMethod("gpolygon",
          signature(x="quadGate", data="character"), 
          function(x, data, channels, verbose=TRUE, col, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
          if(missing(col)){
               col <-  col2rgb(colorRampPalette(brewer.pal(9, "Set1"))(4),
                              alpha=TRUE)
              col <- rgb(col[1,], col[2,], col[3,], 75, maxColorValue=255)
          }
          else
              col <- rep(col,4)
          v <- x@boundary[data[1]]
          h <- x@boundary[data[2]]
          mat <- matrix(c(-Inf, v, h, Inf, v, Inf, h, Inf, -Inf, v, -Inf,
                          h, v, Inf, -Inf, h), byrow=TRUE, ncol=4)              
          for(i in 1:4){
              rg <- rectangleGate(.gate=matrix(mat[i,], ncol=2,
                                  dimnames=list(c("min", "max"), data)))
              gpolygon(rg, verbose=FALSE, col=col[i], ...)
          }
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("gpolygon",
          signature(x="quadGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult", "quadGates", verbose=verbose)
          gpolygon(x, verbose=verbose, ...)
      })

## we can drop the dataFrame, don't need it for rectangleGates
setMethod("gpolygon",
          signature(x="quadGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("flowFrame", "quadGates", verbose=verbose)
          gpolygon(x, verbose=verbose, ...)
      })




## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("gpolygon",
          signature(x="polygonGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
          xp <- x@boundaries[,data[1]]
          yp <- x@boundaries[,data[2]]
          polygon(xp, yp, ...)
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("gpolygon",
          signature(x="polygonGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult", "polygonGates", verbose=verbose)
          gpolygon(x, verbose=verbose, ...)
      })

## we can drop the flowFrame, don't need it for polygonGates
setMethod("gpolygon",
          signature(x="polygonGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("flowFrame", "polygonGates", verbose=verbose)
          gpolygon(x, verbose=verbose, ...)
      })




# ==========================================================================
## for norm2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter
setMethod("gpolygon",
          signature(x="norm2Filter", data="ANY"), 
          function(x, data, verbose=TRUE, ...) evalError("norm2Filters"))

## Filter has been evaluated and the filterResult is provided
setMethod("gpolygon",
          signature(x="norm2Filter", data="logicalFilterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          checkFres(filter=x, fres=data, verbose=verbose)
          fd <- filterDetails(data, identifier(x))
          ## create a polygonGate and plot that
          gpolygon(norm2Polygon(fd, parameters(x)), verbose=verbose, ...)
          })

## Evaluate the filter and plot the filterResult
setMethod("gpolygon",
          signature(x="norm2Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          fres <- filter(data, x)
          gpolygon(x, fres, verbose=verbose, ...)
      })




## ==========================================================================
## for kmeansFilters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We don't know how to plot these, hence we warn
setMethod("gpolygon",
          signature(x="kmeansFilter", data="ANY"), 
          function(x, data, verbose=TRUE, ...){
              nnWarn("kmeansFilter", verbose=verbose)
          })
