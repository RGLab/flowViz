
hexbin <- function(x, y, bins = 50) {
    stopifnot(length(x) == length(y))
    complete <- !is.na(x) & !is.na(y)
    x <- x[complete]
    y <- y[complete]
    if( any(x > 1 ) | any(x < 0 ) ) stop("x values must be in 0,1")
    if( any(y > 1 ) | any(y < 0 ) ) stop("y values must be in 0,1")
    nB = (bins+20)^2
    xbin = rep(0, nB)
    ybin = rep(0, nB)
    counts = rep(0L, nB)
    res = .C("binHex", as.integer(length(x)),  x = as.double(x), 
       y = as.double(y), as.integer(bins), as.integer(nB), 
       counts = as.integer(counts),
       xbin = xbin, ybin = ybin, PACKAGE="flowViz")

   return(list(counts=res$counts, xbin = res$xbin, ybin=res$ybin))
}
