## Store state info in this internal environment
flowViz.state <- new.env(hash = FALSE)
flowViz.state[["plotted"]] <- FALSE

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
## only one argument: the flowFrame
setMethod("plot", signature(x="flowFrame", y="missing"),
          function(x, y, smooth=TRUE, pch, ...)
      {
          l <- ncol(x)
          values <- exprs(x)
          sel <- tolower(colnames(values)) != "time"
          flowViz.state[["plotted"]] <- TRUE
          if(missing(pch))
              pch <- "."
          if(l==1){
              hist(values, xlab=colnames(x), ...)
              flowViz.state[["type"]] <- "hist"
          }else if (l==2){
              if(smooth){
                  smoothScatter(values, pch=pch, ...)
                  flowViz.state[["type"]] <- "smooth"
              }else{
                  plot(values, pch=pch, ...)
                  flowViz.state[["type"]] <- "smooth"
              }
          }else{
              if(smooth){
                  x <- x[,sel]
                  print(splom(x, pch=pch, ...))
                  flowViz.state[["type"]] <- "splom"
              }
              else{
                  pairs(values[,sel], pch=pch, ...)
                  flowViz.state[["type"]] <- "pairs"
              }
          }
          flowViz.state[["parameters"]] <- colnames(x)
          return(invisible(NULL))
      })

## second argument contains the parameters(s) to plot
setMethod("plot",signature(x="flowFrame", y="character"),
          function(x, y, smooth=TRUE, ...){
              if(!all(y %in% colnames(x)))
                  stop("subset out of bounds", call.=FALSE)
              plot(x[,y], smooth=smooth, ...)
              #callGeneric(x[,y], smooth=smooth, ...)
      })
         

