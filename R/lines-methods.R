## ==========================================================================
## Draw gate boundaries as lines on an existing plot. All graphical
## parameters will be passed on to graphics::lines in the end, so
## these are basically extended lines methods for gates and filters.
## Filters are only evaluated if that is needed for plotting of the
## boundaries.
## All auxiliary functions are defined in the file "gateplotting_utils.R"
##
## FIXME: Need to figure out how to make lines a generic without masking
##        the default function
## FIXME: How do we plot transform filters and combined filters?




## ==========================================================================
## for all filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We only know how to add the gate boundaries if the definiton of that gate
## contains exactly two dimensions, and even then we need to guess that
## they match the plotted data, hence we warn. Supplying the names of the
## plotted channels via the "channels" argument always gets precendence.
## All downstream methods will eventually run across this method, so we only
## need to check once.
setMethod("glines",
          signature(x="filter", data="missing"), 
          function(x, data, verbose=TRUE, ...)
      {
          parms <- checkParameterMatch(parameters(x), verbose=verbose, ...)
          glines(x, parms, ...)
      })

## Extract the filter definiton from a filterResult and pass that on
## along with it
setMethod("glines",
          signature(x="filterResult", data="ANY"), 
          function(x, data, verbose=TRUE, ...)
      {
          fd <- filterDetails(x)
          filt <- filterDetails(x)[[length(fd)]]$filter
          if(!missing(data) && is.character(data) &&
             ! ("channels" %in% names(list(...))))
              glines(filt, x, verbose=FALSE, channels=data, ...)
          else
              glines(filt, x, verbose=FALSE, ...)
      })

## We don't need the flowFrame if the filter is already evaluated, but we
## can check that the IDs are matching
setMethod("glines",
          signature(x="filterResult", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          checkIdMatch(x=x, f=data)
          glines(x, verbose=verbose, ...)
          dropWarn("flowFrame", "filterResults", verbose=verbose)
      })




## ==========================================================================
## for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("glines",
          signature(x="rectangleGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
          if(length(parms)==1){ ## ID rectangular gate (region gate)
              mt <- match(parms, data)
              if(mt == 1)
                  abline(v=c(x@min, x@max), ...)
              else if(mt == 2)
                  abline(h=c(x@min, x@max), ...)
              else
                  stop("How did you end up here????")
          }else{ ## 2D rectangular gate
              usr <- par("usr")
              bounds <- cbind(x@min, x@max)[data,]
              bounds[1,] <- fixInf(bounds[1,], usr[1:2])
              bounds[2,] <- fixInf(bounds[2,], usr[3:4])
              bl <- bounds[,1]
              tr <- bounds[,2]
              lines(c(bl[1], tr[1], tr[1], bl[1], bl[1]),
                    c(rep(c(bl[2], tr[2]), each=2), bl[2]), ...) 
          }
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("glines",
          signature(x="rectangleGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult", "rectangleGates", verbose=verbose)
          glines(x, verbose=verbose, ...)
      })

## we can drop the flowFrame, don't need it for rectangleGates
setMethod("glines",
          signature(x="rectangleGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("flowFrame", "rectangleGates", verbose=verbose)
          glines(x, verbose=verbose, ...)
      })




## ==========================================================================
## for quadGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We plot this as ablines
setMethod("glines",
          signature(x="quadGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
          v <- x@boundary[data[1]]
          h <- x@boundary[data[2]]
          abline(v=v,h=h, ...)
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("glines",
          signature(x="quadGate", data="filterResult"), 
          function(x, data, verbose=FALSE, ...)
      {
          dropWarn("filterResult", "quadGates", verbose=verbose)
          glines(x, verbose=verbose, ...)
      })

## we can drop the dataFrame, don't need it for quadGates
setMethod("glines",
          signature(x="quadGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("flowFrame", "quadGates", verbose=verbose)
          glines(x, verbose=verbose, ...)
      })




## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("glines",
          signature(x="polygonGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
          xp <- x@boundaries[,data[1]]
          yp <- x@boundaries[,data[2]]
          lines(c(xp, xp[1]), c(yp,yp[1]), ...)
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("glines",
          signature(x="polygonGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult", "polygonGates", verbose=verbose)
          glines(x, verbose=verbose, ...)
      })

## we can drop the flowFrame, don't need it for polygonGates
setMethod("glines",
          signature(x="polygonGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("flowFrame", "polygonGates", verbose=verbose)
          glines(x, verbose=verbose, ...)
      })




## ==========================================================================
## for norm2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter
setMethod("glines",
          signature(x="norm2Filter", data="ANY"), 
          function(x, data, verbose=TRUE, ...) evalError("norm2Filters"))      

## Filter has been evaluated and the filterResult is provided
setMethod("glines",
          signature(x="norm2Filter", data="logicalFilterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          checkFres(filter=x, fres=data, verbose=verbose)
          fd <- filterDetails(data, identifier(x))
          ## create a polygonGate and plot that
          glines(norm2Polygon(fd, parameters(x)), verbose=verbose, ...)      
      })

## Evaluate the filter and plot the filterResult
setMethod("glines",
          signature(x="norm2Filter", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          fres <- filter(data, x)
          glines(x, fres, verbose=verbose, ...)
      })



## ==========================================================================
## for ellipsoidGatess
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Plotting parameters are specified as a character vector
setMethod("glines",
          signature(x="ellipsoidGate", data="character"), 
          function(x, data, channels, verbose=TRUE, ...)
      {
          if(!missing(channels))
              data <- channels
          parms <- parameters(x)
	  ## We coerce to a polygon gate and plot that
          glines(ell2Polygon(x, parameters(x)), verbose=verbose, ...)      
      })


## We can ignore the filterResult, don't need it to plot the gate
setMethod("glines",
          signature(x="ellipsoidGate", data="filterResult"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("filterResult", "ellipsoidGates", verbose=verbose)
          glines(x, verbose=verbose, ...)
      })

## we can drop the flowFrame, don't need it for ellipsoidGates
setMethod("glines",
          signature(x="ellipsoidGate", data="flowFrame"), 
          function(x, data, verbose=TRUE, ...)
      {
          dropWarn("flowFrame", "ellipsoidGates", verbose=verbose)
          glines(x, verbose=verbose, ...)
      })






## ==========================================================================
## for kmeansFilters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We don't know how to plot these, hence we warn
setMethod("glines",
          signature(x="kmeansFilter", data="ANY"), 
          function(x, data, verbose=TRUE, ...) nnWarn("kmeansFilter", verbose))
