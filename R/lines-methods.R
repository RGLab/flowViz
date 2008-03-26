## ==========================================================================
## Draw gate boundaries as lines on an existing plot. All graphical
## parameters will be passed on to graphics::lines in the end, so
## these are basically extended lines methods for gates and filters.
## Filters are only evaluated if that is needed for plotting of the
## boundaries.
## FIXME: Need to figure out how to make lines a generic without masking
##        the default function



## ==========================================================================
## for all filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We only know how to add the gate boundaries if the definiton of that gate
## contains exactly two dimensions, and even then we need to guess that
## they match the plotted data, hence we warn
setMethod("glines",
          signature(x="filter", data="missing"), 
          definition=function(x, data, ...)
      {
          parms <- parameters(x)
          if(length(parms)!=2)
              stop("The filter definition contains the following parameters:\n",
                   paste(parms, collapse=", "), "\nDon't know how to match to",
                   " the plotted data.\nPlease specify plotting parameters as ",
                   "an additional argument.", call.=FALSE)
          warning("The filter is defined for parameters '",
                  paste(parms, collapse="' and '"), "'.\nPlease make sure ",
                  "that they match the plotting parameters.", call.=FALSE)
          glines(x, parms, ...)
      })




## ==========================================================================
## for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("glines",
          signature(x="rectangleGate", data="character"), 
          definition=function(x, data, ...)
      {
          parms <- parameters(x)
          if(length(data) != 2)
              stop("Plotting parameters need to be character vector of ",
                   "length 2", call.=FALSE)
          mt <- match(parms, data)
          if(any(is.na(mt)))
              stop("The filter is not defined for the following ",
                   "parameter(s):\n", paste(data[is.na(mt)], collapse=", "),
                   call.=FALSE)
          if(length(parms)==1){
              ## one-dimensional rectangular gate (region gate)
              if(mt == 1)
                  abline(v=c(x@min, x@max), ...)
              else if(mt == 2)
                  abline(h=c(x@min, x@max), ...)
              else
                  stop("How did you end up here????")
          }else{
              ## two-dimensional rectangular gate  
              bl <- x@min[parms]
              tr <- x@max[parms]
              lines(c(bl[1], tr[1], tr[1], bl[1], bl[1]),
                    c(rep(c(bl[2], tr[2]), each=2), bl[2]), ...) 
          }
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("glines",
          signature(x="rectangleGate", data="filterResult"), 
          definition=function(x, data, ...)
      {
          warning("No 'filterResult' needed to plot 'rectangleGates'.\n",
                  "Argument is ignored.", call.=FALSE)
          glines(x, ...)
      })



## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("glines",
          signature(x="polygonGate", data="character"), 
          definition=function(x, data, ...)
      {
          parms <- parameters(x)
          if(length(data) != 2)
              stop("Plotting parameters need to be character vector of ",
                   "length 2", call.=FALSE)
          mt <- match(parms, data)
          if(any(is.na(mt)))
              stop("The filter is not defined for the following ",
                   "parameter(s):\n", paste(data[is.na(mt)], collapse=", "),
                   call.=FALSE)
          xp <- x@boundaries[,mt[1]]
          yp <- x@boundaries[,mt[2]]
          lines(c(xp, xp[1]), c(yp,yp[1]), ...)
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("glines",
          signature(x="polygonGate", data="filterResult"), 
          definition=function(x, data, ...)
      {
          warning("No 'filterResult' needed to plot 'polygonGates'.\n",
                  "Argument is ignored.", call.=FALSE)
          glines(x, ...)
      })




## ==========================================================================
## for norm2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter
setMethod("glines",
          signature(x="norm2Filter", data="ANY"), 
          definition=function(x, data, ...)  
          stop("'norm2Filters' need to be evaluated for plotting.\n",
               "Either provide a 'flowFrame' or an appropriate ",
               "'filterResult' \nas second argument.", call.=FALSE) 
          )

## Filter has been evaluated and the filterResult is provided
setMethod("glines",
          signature(x="norm2Filter", data="logicalFilterResult"), 
          definition=function(x, data, ...){
              ## a lot of sanity checking up first
              fd <- filterDetails(data, identifier(x))
              if(!identical(identifier(x), identifier(data)) ||
                 class(x) != class(fd$filter))
                  stop("The 'filterResult' and the 'norm2Filter' ",
                       "don't match.", call.=FALSE)
              parms <- parameters(x)
              if(length(fd$filter@transformation) > 0)
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
              glines(polygonGate(boundaries=ans), ...)          
          })

## Evaluate the filter and plot the filterResult
setMethod("glines",
          signature(x="norm2Filter", data="flowFrame"), 
          definition=function(x, data, ...){
              fres <- filter(data, x)
               glines(x, fres, ...)
          })


## ==========================================================================
## for curv2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter
setMethod("glines",
          signature(x="curv2Filter", data="ANY"), 
          definition=function(x, data, ...)  
          stop("'curv2Filters' need to be evaluated for plotting.\n",
               "Either provide a 'flowFrame' or an appropriate ",
               "'filterResult' \nas second argument.", call.=FALSE) 
          )

## Filter has been evaluated and the filterResult is provided
setMethod("glines",
          signature(x="curv2Filter", data="multipleFilterResult"), 
          definition=function(x, data, ...){
              ## a lot of sanity checking up first
              fd <- filterDetails(data, identifier(x))
              if(!identical(identifier(x), identifier(data)) ||
                 class(x) != class(fd$filter))
                  stop("The 'filterResult' and the 'norm2Filter' ",
                       "don't match.", call.=FALSE)
              parms <- parameters(x)
              polygons <- fd$polygons
              oo <- options(warn=-1)
              sapply(polygons,
                     function(x, ...){
                         tmp <- cbind(x$x, x$y)
                         colnames(tmp) <- parms
                         glines(polygonGate(boundaries=tmp), ...)
                     }, ...)
              options(oo)
              warning("The filter is defined for parameters '",
                      paste(parms, collapse="' and '"), "'.\nPlease make sure ",
                      "that they match the plotting parameters.", call.=FALSE)
              return(NULL)
          })
              
## Evaluate the filter and plot the filterResult
setMethod("glines",
          signature(x="curv2Filter", data="flowFrame"), 
          definition=function(x, data, ...){
              fres <- filter(data, x)
               glines(x, fres, ...)
          })
