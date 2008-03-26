## ==========================================================================
## Draw gate regions as polygons on an existing plot. All graphical
## parameters will be passed on to graphics::polygon in the end, so
## these are basically extended polygon methods for gates and filters.
## Filters are only evaluated if that is needed for plotting of the
## regions.
## FIXME: Need to figure out how to make polygon a generic without masking
##        the default function


## replace infinite values by something useful
fixInf <- function(x, replacement){
    if(is.infinite(x) && x<0)
        replacement[1]
    else if(is.infinite(x) && x>0)
        replacement[2]
    else
        x
}
    

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
## they match the plotted data, hence we warn
setMethod("gpolygon",
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
          gpolygon(x, parms, ...)
      })



## ==========================================================================
## for rectangleGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("gpolygon",
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
          usr <- par("usr")
          if(length(parms)==1){
              ## one-dimensional rectangular gate
              if(mt==1)
                  rect(fixInf(x@min, usr[1]), usr[3], fixInf(x@max, usr[2]),
                       usr[4], ...)
              else if(mt==2)
                  rect(usr[1], fixInf(x@min, usr[3]), usr[2],
                       fixInf(x@max, usr[4]), ...)
              else
                  stop("How did you end up here????")
          }else{
              ## two-dimensional rectangular gate    
              bl <- x@min[parms]
              tr <- x@max[parms]
              rect(fixInf(bl[1], usr[1]), fixInf(bl[2], usr[3]),
                   fixInf(tr[1], usr[2]),  fixInf(tr[2], usr[4]), ...)
          }
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("gpolygon",
          signature(x="rectangleGate", data="filterResult"), 
          definition=function(x, data, ...)
      {
          warning("No 'filterResult' needed to plot 'rectangleGates'.\n",
                  "Argument is ignored.", call.=FALSE)
          gpolygon(x, ...)
      })



## ==========================================================================
## for polygonGates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plotting parameters are specified as a character vector
setMethod("gpolygon",
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
          polygon(xp, yp, ...)
      })

## We can ignore the filterResult, don't need it to plot the gate
setMethod("gpolygon",
          signature(x="polygonGate", data="filterResult"), 
          definition=function(x, data, ...)
      {
          warning("No 'filterResult' needed to plot 'polygonGates'.\n",
                  "Argument is ignored.", call.=FALSE)
          gpolygon(x, ...)
      })


# ==========================================================================
## for norm2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter
setMethod("gpolygon",
          signature(x="norm2Filter", data="ANY"), 
          definition=function(x, data, ...)  
          stop("'norm2Filters' need to be evaluated for plotting.\n",
               "Either provide a 'flowFrame' or an appropriate ",
               "'filterResult' \nas second argument.", call.=FALSE) 
          )

## Filter has been evaluated and the filterResult is provided
setMethod("gpolygon",
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
              gpolygon(polygonGate(boundaries=ans), ...)          
          })

## Evaluate the filter and plot the filterResult
setMethod("gpolygon",
          signature(x="norm2Filter", data="flowFrame"), 
          definition=function(x, data, ...){
              fres <- filter(data, x)
               gpolygon(x, fres, ...)
          })



## ==========================================================================
## for curv2Filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## An error if we can't evaluate the filter
setMethod("gpolygon",
          signature(x="curv2Filter", data="ANY"), 
          definition=function(x, data, ...)  
          stop("'curv2Filters' need to be evaluated for plotting.\n",
               "Either provide a 'flowFrame' or an appropriate ",
               "'filterResult' \nas second argument.", call.=FALSE) 
          )

## Filter has been evaluated and the filterResult is provided
setMethod("gpolygon",
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
                         gpolygon(polygonGate(boundaries=tmp), ...)
                     }, ...)
              options(oo)
              warning("The filter is defined for parameters '",
                      paste(parms, collapse="' and '"), "'.\nPlease make sure ",
                      "that they match the plotting parameters.", call.=FALSE)
              return(NULL)
          })
              
## Evaluate the filter and plot the filterResult
setMethod("gpolygon",
          signature(x="curv2Filter", data="flowFrame"), 
          definition=function(x, data, ...){
              fres <- filter(data, x)
               gpolygon(x, fres, ...)
          })
