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


## ===========================================================================
## Our own generic for lines:
## ---------------------------------------------------------------------------
setGeneric("glines", function(x, data, ...) standardGeneric("glines"))


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
## flowPlot:
## ---------------------------------------------------------------------------
setGeneric("flowPlot",function(x,...) standardGeneric("flowPlot"))


## ===========================================================================
## contour plots:
## ---------------------------------------------------------------------------
setGeneric("contour")


## ===========================================================================
## Generics for all plot types defined in lattice:
## ---------------------------------------------------------------------------
setGeneric("xyplot")
setGeneric("densityplot")
setGeneric("ecdfplot")
setGeneric("levelplot")
setGeneric("qqmath")
setGeneric("splom")
setGeneric("parallel")
