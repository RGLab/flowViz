## ==========================================================================
## a simple plot method without strange plot parameter and friends. It does
## the most intuitive thing: take a flowFrame and do a pairs plot if there
## are more than 2 parameters. Do a smoothScatter plot for exactly 2
## parameters and a histogram for exactly one.
## If you want to plot specific columns, subset the frame before
## plotting or specify them in the y argument.
## We record parameters of the plot in the internal flowViz.state environment
## which will later be useful to figure out the plotting parameters when
## adding gate outlines to a plot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## helper function to do the actual plotting
fplot <- function(x, smooth, pch, xlim, ylim = NULL, ...)
{
    values <- exprs(x)
    sel <- tolower(colnames(values)) != "time"
#    flowViz.state[["plotted"]] <- TRUE
    if(missing(pch))
        pch <- "."
    l <- ncol(values)
    if(l==1){
        if (missing(xlim))
            xlim <- unlist(range(x[,1]))
        hist(values, xlab=colnames(x), xlim = xlim, ylim = ylim, ...)
#        flowViz.state[["type"]] <- "hist"
    }
    else if (l==2) {
        if (missing(xlim))
            xlim <- unlist(range(x[,1]))
        if (is.null(ylim))
            ylim <- unlist(range(x[,2]))
        if (smooth) {
            smoothScatter(values, pch=pch, xlim=xlim, ylim=ylim, ...)
#            flowViz.state[["type"]] <- "smooth"
        }else{
            plot(values, pch=pch, xlim=xlim, ylim=ylim, ...)
#            flowViz.state[["type"]] <- "dot"
        }
    } else {
        if(smooth) {
            x <- x[,sel]
            print(splom(x, pch=pch, ...))
#            flowViz.state[["type"]] <- "splom"
        }
        else {
            pairs(values[,sel], pch=pch, ...)
#            flowViz.state[["type"]] <- "pairs"
        }
    }
#    flowViz.state[["parameters"]] <- colnames(values)
}


## only one argument: the flowFrame
setMethod("plot", signature(x="flowFrame", y="missing"),
          function(x, y, smooth=TRUE, ...)
      {
          fplot(x, smooth=smooth, ...)
          return(invisible(NULL))
      })

## second argument contains the parameters(s) to plot
setMethod("plot",signature(x="flowFrame", y="character"),
          function(x, y, smooth=TRUE, pch, ...)
      {
          if(!all(y %in% colnames(x)))
              stop("subset out of bounds", call.=FALSE)
          x <- x[,y]
          fplot(x, smooth=smooth, pch=pch, ...)
          return(invisible(NULL))
      })


