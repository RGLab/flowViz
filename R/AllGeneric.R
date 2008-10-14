## =========================================================================##
## =========================================================================##
##                       Generic methods definition                         ##
## =========================================================================##
## =========================================================================##





## ===========================================================================
## Our own generic for polygons:
## ---------------------------------------------------------------------------
setGeneric("gpolygon", function(x, data, ...)
           standardGeneric("gpolygon"))
setGeneric("glpolygon", function(x, data, ...)
           standardGeneric("glpolygon"))
#setGeneric("polygon", function(x, ...) standardGeneric("polygon"), 
#           useAsDefault=FALSE)


## ===========================================================================
## Our own generic for lines:
## ---------------------------------------------------------------------------
setGeneric("glines", function(x, data, ...) standardGeneric("glines"))
setGeneric("gllines", function(x, data, ...) standardGeneric("gllines"))


## ===========================================================================
## Our own generic for points:
## ---------------------------------------------------------------------------
setGeneric("gpoints", function(x, data, channels, ...)
           standardGeneric("gpoints"))
setGeneric("glpoints", function(x, data, channels, ...)
           standardGeneric("glpoints"))


## ===========================================================================
## Our own generic to add gate names
## ---------------------------------------------------------------------------
setGeneric("addName", function(x, name, ...)
           standardGeneric("addName"))



## ===========================================================================
## timeLinePlot:
## ---------------------------------------------------------------------------
setGeneric("timeLinePlot",
           function(x, channel, ...) standardGeneric("timeLinePlot"))
          


## ===========================================================================
## contour plots:
## ---------------------------------------------------------------------------
setGeneric("contour", function(x,...) standardGeneric("contour"))



## ===========================================================================
## Generics for all plot types defined in lattice:
## ---------------------------------------------------------------------------
setGeneric("xyplot", function(x, data, ...)
           standardGeneric("xyplot"))
setGeneric("densityplot", function(x, data, ...)
           standardGeneric("densityplot"))
setGeneric("ecdfplot", function(x, data, ...)
           standardGeneric("ecdfplot"))
setGeneric("levelplot", function(x, data, ...)
           standardGeneric("levelplot"))
setGeneric("qqmath", function(x, data, ...)
           standardGeneric("qqmath"))
setGeneric("splom", function(x, data, ...)
           standardGeneric("splom"))
setGeneric("parallel", function(x, data, ...)
           standardGeneric("parallel"))
