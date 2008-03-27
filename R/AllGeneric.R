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

