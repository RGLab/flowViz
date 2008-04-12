## ==========================================================================
## Draw gate regions as polygons on an existing lattice plot. All graphical
## parameters will be passed on to lpolygon in the end, so
## these are basically extended lpolygon methods for gates and filters.
## Filters are only evaluated if that is needed for plotting of the
## regions.
## FIXME: Need to figure out how to make polygon a generic without masking
##        the default function



## ==========================================================================
## for all filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We only know how to add the gate boundaries if the definiton of that gate
## contains exactly two dimensions, and even then we need to guess that
## they match the plotted data, hence we warn
setMethod("glpolygon",
          signature(x="filter", data="missing"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(missing(channels))
              parms <- checkParameterMatch(parameters(x), verbose=verbose)
          else
              parms <- checkParameterMatch(channels, verbose=FALSE)
          glpolygon(x, parms, ...)
      })


## Extract the filter definiton from a filterResult and pass that on
## along with it
setMethod("glpolygon",
          signature(x="filterResult", data="missing"), 
          definition=function(x, data, verbose=TRUE, ...)
      {
          filt <- filterDetails(x)[[1]]$filter
          parms <- checkParameterMatch(parameters(filt))
          glpolygon(filt, x, verbose=FALSE, ...)
      })


## We don't need the flowFrame if the filter is already evaluated, but we
## can check that the IDs are matching
setMethod("glpolygon",
          signature(x="filterResult", data="flowFrame"), 
          definition=function(x, data, verbose=TRUE, channels, ...){
               checkIdMatch(x=x, f=data)
               glpolygon(x, verbose=verbose, channels=channels, ...)
               dropWarn("flowFrame", "filterResults", verbose=verbose)
          })




## ==========================================================================
## for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("glpolygon",
          signature(x="rectangleGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          data <- checkParameterMatch(data, verbose=FALSE)
          parms <- parameters(x)
          usr <- rep(c(-1e4, 1e4), 2)
          if(length(parms)==1){
              mt <- match(parms, data)## 1D rectangular gate (region gate)
              if(mt==1)
                  lrect(fixInf(x@min, usr[1:2]), usr[3],
                        fixInf(x@max, usr[1:2]), usr[4], ...)
              else if(mt==2)
                  lrect(usr[1], fixInf(x@min, usr[3:4]), usr[2],
                       fixInf(x@max, usr[3:4]), ...)
              else
                  stop("How did you end up here????")
          }else{## 2D rectangular gate   
              bl <- x@min[data]
              tr <- x@max[data]
              lrect(fixInf(bl[1], usr[1:2]), fixInf(bl[2], usr[3:4]),
                   fixInf(tr[1], usr[1:2]),  fixInf(tr[2], usr[3:4]), ...)
          }
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("glpolygon",
          signature(x="rectangleGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult", "rectangleGates", verbose=verbose)
          glpolygon(x, verbose=verbose, ...)
      })

## we can drop the flowFrame, don't need it for rectangleGates
setMethod("glpolygon",
          signature(x="rectangleGate", data="flowFrame"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          dropWarn("flowFrame", "rectangleGates", verbose=verbose)
          if(!missing(channels))
              glpolygon(x, channels, verbose=verbose, ...)
          else
              glpolygon(x, verbose=verbose, ...)
      })





## ==========================================================================
## for quadGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We could plot this as four individual rectangle gates, but just using
## abline makes more sense
setMethod("glpolygon",
          signature(x="quadGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          data <- checkParameterMatch(data, verbose=FALSE)
          parms <- parameters(x)
          v <- x@boundary[data[1]]
          h <- x@boundary[data[2]]
          panel.abline(v=v, h=h, ...)
          ## the previous doesn't allow for a fill color, but I am not sure
          ## that we actually need that. The following should work in case we do
          ##  mat <- matrix(c(-Inf, v, h, Inf, v, Inf, h, Inf, -Inf, v, -Inf,
          ##                           h, v, Inf, -Inf, h), byrow=TRUE, ncol=4)              
          ##           for(i in 1:4){
          ##               rg <- rectangleGate(.gate=matrix(mat[i,], ncol=2,
          ##                                   dimnames=list(c("min", "max"), parms)))
          ##               glpolygon(rg, verbose=FALSE, col=col[i], ...)
          ##           }
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("glpolygon",
          signature(x="quadGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult", "quadGates", verbose=verbose)
          glpolygon(x, verbose=verbose, ...)
      })

## we can drop the dataFrame, don't need it for rectangleGates
setMethod("glpolygon",
          signature(x="quadGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, channels, ...)
      {
          dropWarn("flowFrame", "quadGates", verbose=verbose)
          if(!missing(channels))
              glpolygon(x, channels, verbose=verbose, ...)
          else
              glpolygon(x, verbose=verbose, ...)
      })




## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("glpolygon",
          signature(x="polygonGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          data <- checkParameterMatch(data, verbose=FALSE)
          parms <- parameters(x)
          xp <- x@boundaries[,data[1]]
          yp <- x@boundaries[,data[2]]
          lpolygon(xp, yp, ...)
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("glpolygon",
          signature(x="polygonGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult", "polygonGates", verbose=verbose)
          glpolygon(x, verbose=verbose, ...)
      })

## we can drop the flowFrame, don't need it for polygonGates
setMethod("glpolygon",
          signature(x="polygonGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, channels, ...)
      {
          dropWarn("flowFrame", "polygonGates", verbose=verbose)
          if(!missing(channels))
              glpolygon(x, channels, verbose=verbose, ...)
          else
              glpolygon(x, verbose=verbose, ...)
      })




# ==========================================================================
## for norm2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter
setMethod("glpolygon",
          signature(x="norm2Filter", data="ANY"), 
          function(x, data, verbose=TRUE, ...) evalError("norm2Filters"))

## Filter has been evaluated and the filterResult is provided
setMethod("glpolygon",
          signature(x="norm2Filter", data="logicalFilterResult"), 
          definition=function(x, data, verbose=TRUE, ...){
              ## a lot of sanity checking up first
              fd <- filterDetails(data, identifier(x))
              if(!identical(identifier(x), identifier(data)) ||
                 class(x) != class(fd$filter))
                  stop("The 'filterResult' and the 'norm2Filter' ",
                       "don't match.", call.=FALSE)
              parms <- parameters(x)
              if(length(fd$filter@transformation) > 0 && verbose)
                  warning("'result' appears to have been applied on ",
                          "transformed data.\nThese are not supported yet.")
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
              ## create a polygonGate and plot that 
              glpolygon(polygonGate(boundaries=ans), verbose=verbose, ...)          
          })

## Evaluate the filter and plot the filterResult
setMethod("glpolygon",
          signature(x="norm2Filter", data="flowFrame"), 
          definition=function(x, data, verbose=TRUE, ...){
              fres <- filter(data, x)
               glpolygon(x, fres, verbose=verbose, ...)
          })



## ==========================================================================
## for curv2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter
setMethod("glpolygon",
          signature(x="curv2Filter", data="ANY"), 
          function(x, data, verbose=TRUE, ...)  evalError("curv2Filters")) 
         
## Filter has been evaluated and the filterResult is provided
setMethod("glpolygon",
          signature(x="curv2Filter", data="multipleFilterResult"), 
          function(x, data, verbose=TRUE, col, ...)
      {
          checkFres(filter=x, fres=data, verbose=verbose)
          fd <- filterDetails(data, identifier(x))
          parms <- parameters(x)
          polygons <- fd$polygons
          lf <- length(polygons)
          ## we want to use different colors for each population
          if(missing(col))
              col <-  colorRampPalette(brewer.pal(9, "Set1"))(lf)
          else
              col <- rep(col, lf)[1:lf]
          mapply(function(x, co, ...){
              tmp <- cbind(x$x, x$y)
              colnames(tmp) <- parms
              glpolygon(polygonGate(boundaries=tmp), col=co, ...)
          }, x=polygons, co=col, MoreArgs=list(verbose=FALSE, ...))
          return(invisible(NULL))
      })

## Evaluate the filter and plot the filterResult
setMethod("glpolygon",
          signature(x="curv2Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          fres <- filter(data, x)
          glpolygon(x, fres, verbose=verbose, ...)
      })



## ==========================================================================
## for curv1Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter
setMethod("glpolygon",
          signature(x="curv1Filter", data="ANY"), 
          function(x, data, verbose=TRUE, ...) evalError("curv1Filters"))
setMethod("glpolygon",
          signature(x="curv1Filter", data="missing"), 
          function(x, data, verbose=TRUE, ...) evalError("curv1Filters"))  

## Filter has been evaluated and the filterResult is provided
setMethod("glpolygon",
          signature(x="curv1Filter", data="multipleFilterResult"), 
          function(x, data, channels, verbose=TRUE, col, ...)
      {
          checkFres(filter=x, fres=data, verbose=verbose)
          fd <- filterDetails(data, identifier(x))
          parms <- parameters(x)
          if(missing(channels))
              channels <- parms
          channels <- checkParameterMatch(channels, verbose=FALSE)
          bounds <- fd$boundaries
          lb <- length(bounds)
          ## we want to use different colors for each population
          if(missing(col))
              col <-  colorRampPalette(brewer.pal(9, "Set1"))(lb)
          else
              col <- rep(col, lb)[1:lb]
          mapply(function(x, co, ...){
              tmp <- matrix(x, nrow=2)
              colnames(tmp) <- parms
              glpolygon(rectangleGate(.gate=tmp), channels, col=co, ...)
          }, x=bounds, co=col, MoreArgs=list(verbose=FALSE, ...))
          return(invisible(NULL))
      })
              
## Evaluate the filter and plot the filterResult
setMethod("glpolygon",
          signature(x="curv1Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          fres <- filter(data, x)
          glpolygon(x, fres, verbose=verbose, ...)
      })



## ==========================================================================
## for kmeansFilters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We don't know how to plot these, hence we warn
setMethod("glpolygon",
          signature(x="kmeansFilter", data="ANY"), 
          function(x, data, verbose=TRUE, ...) nnWarn("kmeansFilter", verbose))
