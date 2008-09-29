
hexbin <- function(x, y, bins = 50) {
    stopifnot(length(x) == length(y))
    complete <- !is.na(x) & !is.na(y)
    x <- x[complete]
    y <- y[complete]
    if( any(x > 1 ) | any(x < 0 ) ) stop("x values must be in 0,1")
    if( any(y > 1 ) | any(y < 0 ) ) stop("y values must be in 0,1")
    nB = (bins+20)^2
    res = .Call("binHex", as.integer(length(x)),  x = as.double(x), 
       y = as.double(y), as.integer(bins), as.integer(nB), 
       counts = vector("integer", length = nB), 
       xbin = vector("double", length=nB),
       ybin = vector("double", length=nB), PACKAGE="flowViz")

   return(list(counts=res$counts, xbin = res$xbin, ybin=res$ybin))
}
