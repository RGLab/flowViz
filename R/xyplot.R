



setMethod("xyplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab, ylab,
                   as.table = TRUE,
                   pch = ".", smooth = TRUE,
                   filter = NULL,
                   filterResults = NULL,
                   ...)
      {

          
          if (!is.null(filter) && is.null(filterResults))
              filterResults <- filter(data, filter)

          

          pd <- pData(phenoData(data))
          uniq.name <- createUniqueColumnName(pd)
          ## ugly hack to suppress warnings about coercion introducing
          ## NAs (needs to be `undone' inside prepanel and panel
          ## functions):
          pd[[uniq.name]] <- factor(sampleNames(data)) 
          channel.y <- x[[2]]
          channel.x <- x[[3]]
          if (length(channel.x) == 3)
          {
              channel.x <- channel.x[[2]]
              x[[3]][[2]] <- as.name(uniq.name)
              x[[2]] <- NULL
          }
          else
          {
              x[[3]] <- as.name(uniq.name)
              x[[2]] <- NULL
          }
          ## channel.x <- as.character(channel.x)
          ## channel.y <- as.character(channel.y)
          channel.x <- as.expression(channel.x)
          channel.y <- as.expression(channel.y)

          prepanel.xyplot.flowset <- 
              function(x, 
                       frames, channel.x, channel.y,
                       ...)
              {
                  if (length(nm <- as.character(x)) > 1)
                      stop("must have only one flow frame per panel")
                  if (length(nm) == 1)
                  {
                      xx <- evalInFlowFrame(channel.x, frames[[nm]])
                      yy <- evalInFlowFrame(channel.y, frames[[nm]])
                      list(xlim = range(xx, finite = TRUE),
                           ylim = range(yy, finite = TRUE),
                           dx = diff(xx), dy = diff(yy))
                  }
                  else list()
              }

          panel.xyplot.flowset <-
              function(x, 
                       frames, channel.x, channel.y,
                       filterResults = NULL,

                       pch, smooth,
                       ...)
              {
                  x <- as.character(x)
                  if (length(x) > 1) stop("must have only one flow frame per panel")
                  if (length(x) < 1) return()

                  nm <- x
                  xx <- evalInFlowFrame(channel.x, frames[[nm]])
                  yy <- evalInFlowFrame(channel.y, frames[[nm]])
##                   xx <- exprs(frames[[nm]])[, channel.x]
##                   yy <- exprs(frames[[nm]])[, channel.y]
                  

                  if (!is.null(filterResults))
                  {
                      this.filter.result <- filterResults[[nm]]
                      groups <- this.filter.result@subSet

                      ## do something special if there's a natural
                      ## representation of the filter/gate in this
                      ## panel.  When is that true?  First of all, the
                      ## filter has to be a "2-D" filter with no holes
                      ## (topologically speaking), and it's "x" and
                      ## "y" variables must match the current
                      ## channel.x and channel.y

                      ## to complicate things, both the filter and the
                      ## data being plotted may be transformed, not
                      ## necessarily the same way.  One option, that
                      ## should work most of the time, is to simply
                      ## draw a convex hull around 'groups=TRUE' when
                      ## the first set of requirements are met.
                      
                      if (TRUE)
                      {
                          hull <- chull(xx[groups], yy[groups])
                      }

                      
                  }
                  else groups <- NULL

                  if (smooth) {
                      panel.smoothScatter(xx, yy, ...)
                      if (!is.null(groups))
                      {
                          lpolygon(xx[groups][hull],
                                   yy[groups][hull],
                                   border = "yellow")
                      }
                  }
                  else panel.xyplot(xx, yy, pch = pch,
                                    groups = groups,
                                    subscripts = seq_along(groups),
                                    ...)
              }

          if (missing(xlab)) xlab <- as.character(channel.x)
          if (missing(ylab)) ylab <- as.character(channel.y)

          densityplot(x, data = pd, 

                      prepanel = prepanel.xyplot.flowset,
                      panel = panel.xyplot.flowset,

                      frames = data@frames,
                      channel.x = channel.x,
                      channel.y = channel.y,
                      filterResults = filterResults,
                      as.table = as.table,

                      xlab = xlab,
                      ylab = ylab,
                      pch = pch, smooth = smooth,

                      ...)
      })






## xyplot for flowFrames, plotting columns against time

setMethod("xyplot",
          signature(x = "flowFrame", data = "missing"),
          function(x, data, xlab = time, ylab = "", time = "Time",
                   layout,
                   type = "l", smooth = FALSE,
                   ...)
      {
          expr <- exprs(x)
          if (!(time %in% colnames(expr)))
              stop("Name of the Time (X) variable must be specified as the 'time' argument")
          time.x <- expr[, time]
          fakedf <- data.frame(channel = setdiff(colnames(expr), time), time = 1)
          prepanel.xyplot.flowframe <- 
              function(x, y, time.x, expr, ...)
              {
                  xx <- time.x
                  yy <- expr[, as.character(y)]
                  list(xlim = range(xx, finite = TRUE),
                       ylim = range(yy, finite = TRUE),
                       dx = diff(xx), dy = diff(yy))
              }
          panel.xyplot.flowframe <- 
              function(x, y, time.x, expr, smooth = FALSE, ...)
              {
                  xx <- time.x
                  yy <- expr[, as.character(y)]
                  if (smooth) panel.smoothScatter(xx, yy, ...)
                  else panel.xyplot(xx, yy, ...)
              }
          if (missing(layout)) layout <- c(1, ncol(expr) - 1)
          xyplot(channel ~ time | channel, data = fakedf,
                 prepanel = prepanel.xyplot.flowframe,
                 panel = panel.xyplot.flowframe,
                 type = type, smooth = smooth,
                 time.x = time.x, expr = expr,
                 xlab = xlab, ylab = ylab,
                 layout = layout,
                 default.scales = list(y = list(relation = "free", rot = 0)),
                 ...)
      })


## xyplot with a formula.  We'll make this very simple; the drawback
## being that the expression matrix will be copied, the upshot being
## that all the fancy xyplot formula stuff will be valid.

setMethod("xyplot",
          signature(x = "formula", data = "flowFrame"),
          function(x, data,
                   smooth = TRUE, 
                   panel = if (smooth) panel.smoothScatter else panel.xyplot,
                   ...)
      {
          xyplot(x, data = as.data.frame(exprs(data)),
                 smooth = smooth, 
                 panel = panel,
                 ...)
      })


