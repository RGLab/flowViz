#' @importFrom MASS kde2d
#' @export 
#' @rdname lattice-methods
setMethod("levelplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab, ylab,
                   as.table = TRUE, contour = TRUE, labels = FALSE,
                   n = 50, 
                   ...)
      {
          samples <- sampleNames(data)
          channel.y <- x[[2]]
          channel.x <- x[[3]]
          if (length(channel.x) == 3) ## cy ~ cx | foo+bar
          {
              channel.x <- channel.x[[2]]
              my.formula <- x
              my.formula[[2]] <- as.name("z")
              my.formula[[3]][[2]] <- (z ~ x * y)[[3]]
          }
          else ## cy ~ cx
          {
              my.formula <- z ~ x * y
          }
          channel.x.name <- expr2char(channel.x)
          channel.y.name <- expr2char(channel.y)
          channel.x <- as.expression(channel.x)
          channel.y <- as.expression(channel.y)

          range.x <- range(unlist(lapply(samples, function(nm) range(evalInFlowFrame(channel.x, data@frames[[nm]]), finite = TRUE))))
          range.y <- range(unlist(lapply(samples, function(nm) range(evalInFlowFrame(channel.y, data@frames[[nm]]), finite = TRUE))))
          ## range.y <- range(unlist(eapply(data@frames, function(x) range(x[, channel.y], finite = TRUE))))

          range.x <- lattice:::extend.limits(range.x, prop = 0.05)
          range.y <- lattice:::extend.limits(range.y, prop = 0.05)
          computeKde <- function(nm)
          {
              xx <- evalInFlowFrame(channel.x, data@frames[[nm]])
              yy <- evalInFlowFrame(channel.y, data@frames[[nm]])
              temp <- kde2d(xx, yy, n = n, lims = c(range.x, range.y))
              data.frame(x = rep(temp$x, n),
                         y = rep(temp$y, each = n),
                         z = as.vector(temp$z))
          }
          ft <- do.call(make.groups, sapply(samples, computeKde, simplify = FALSE))
          ## The next step is not efficient, as it replicates the
          ## phenodata many times over.  However, bypassing this is a
          ## task for another day
          ft <- cbind(ft, pData(phenoData(data))[as.character(ft$which), , drop = FALSE])
          if (missing(xlab)) xlab <- channel.x.name
          if (missing(ylab)) ylab <- channel.y.name
          levelplot(my.formula, data = ft,
                    contour = contour, labels = labels,
                    xlab = xlab,
                    ylab = ylab,
                    as.table = as.table, 
                    ...)
      })


##Example:
##
## require(colorspaces)
## YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
## colori <- colorRampPalette(YlOrBr, space = "Lab")
## levelplot(`SSC-H` ~ `FSC-H` | name, samp, col.regions=colori(50), main="Contour Plot")
##



## QQ plot

## setGeneric("qqmath")


#' Methods implementing Lattice displays for flow data
#' 
#' 
#' Various methods implementing multipanel visualizations for flow data using
#' infrastructure provided in the lattice package.  The original generics for
#' these methods are defined in lattice, and these S4 methods (mostly) dispatch
#' on a formula and the \code{data} argument which must be of class
#' \code{flowSet} or \code{flowFrame}.  The formula has to be fairly basic:
#' conditioning can be done using phenodata variables and channel names (the
#' \code{colnames} slot) can be used as panel variables. See examples below for
#' sample usage.
#' 
#' Not all standard lattice arguments will have the intended effect, but many
#' should.  For a fuller description of possible arguments and their effects,
#' consult documentation on lattice (Trellis docs would also work for the
#' fundamentals).
#' 
#' @name lattice-methods
#' @param x a formula describing the structure of the plot and the variables to
#' be used in the display.
#' @param data a \code{flowSet} object that serves as a source of data.
#' @param xlab,ylab Labels for data axes, with suitable defaults taken from the
#' formula
#' @param f.value,distribution number of points used in Q-Q plot, and the
#' reference distribution used.  See \code{\link[lattice:qqmath]{qqmath}} for
#' details.
#' @param n the number of bins on each axis to be used when evaluating the
#' density
#' @param as.table,contour,labels These arguments are passed unchanged to the
#' corresponding methods in lattice, and are listed here only because they
#' provide different defaults.  See documentation for the original methods for
#' details.
#' @param time A character string giving the name of the column recording time.
#' @param exclude.time logical, specifying whether to exclude the time variable
#' from a scatter plot matrix or parallel coordinates plot.  It is rarely
#' meaningful not to do so.
#' @param reorder.by a function, which is applied to each column.  The columns
#' are ordered by the results.  Reordering can be suppressed by setting this to
#' \code{NULL}.
#' @param \dots more arguments, usually passed on to the underlying lattice
#' methods.
#' @section Methods: \describe{
#' 
#' \item{qqmath}{\code{signature(x = "formula", data = "flowSet")}: creates
#' theoretical quantile plots of a given channel, with one or more samples per
#' panel }
#' 
#' \item{levelplot}{\code{signature(x = "formula", data = "flowSet")}: similar
#' to the \code{xyplot} method, but plots estimated density (using
#' \code{\link[MASS:kde2d]{kde2d}}) with a common z-scale and an optional color
#' key.  }
#' 
#' \item{parallel}{\code{signature(x = "flowFrame", data = "missing")}: draws a
#' parallel coordinates plot of all channels (excluding time, by default) of a
#' \code{flowFrame} object.  This is rarely useful without transparency, but
#' that is currently only possible with the \code{\link{pdf}} device (and
#' perhaps the aqua device as well).  }
#' 
#' }
#' @keywords methods dplot
#' @examples
#' 
#' 
#' data(GvHD)
#' 
#' qqmath( ~ `FSC-H` | factor(Patient), GvHD,
#'        grid = TRUE, type = "l",
#'        f.value = ppoints(100))
#' 
#' 
#' ## contourplot of bivariate density:
#' 
#' require(colorspace)
#' YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
#' colori <- colorRampPalette(YlOrBr)
#' levelplot(asinh(`SSC-H`) ~ asinh(`FSC-H`) | Visit + Patient, GvHD, n = 20,
#'           col.regions = colori(50), main = "Contour Plot")
#' 
#' 
#' 
#' ## parallel coordinate plots
#' 
#' parallel(GvHD[["s6a01"]])
#' 
#' \dontrun{
#' 
#' ## try with PDF device
#' parallel(GvHD[["s7a01"]], alpha = 0.01)
#' 
#' }
#' 
#' @importFrom lattice densityplot histogram levelplot make.groups panel.abline panel.grid panel.lines panel.parallel panel.points panel.smoothScatter panel.text panel.xyplot parallel prepanel.default.parallel qqmath splom trellis.par.get trellis.par.set which.packet xyplot standard.theme
#' @export 
#' @rdname lattice-methods
#' @aliases qqmath,formula,flowSet-method
setMethod("qqmath",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab, ylab,
                   f.value = function(n) ppoints(ceiling(sqrt(n))),
                   distribution = qnorm,
                   ...)
      {
          pd <- pData(phenoData(data))
          uniq.name <- createUniqueColumnName(pd)
          ## ugly hack to suppress warnings about coercion introducing NAs
          pd[[uniq.name]] <- factor(sampleNames(data))
          channel <- x[[2]]
          if (length(channel) == 3)
          {
              channel <- channel[[2]]
              x[[2]][[2]] <- as.name(uniq.name)
          }
          else x[[2]] <- as.name(uniq.name)
          channel.name <- expr2char(channel)
          channel <- as.expression(channel)

          prepanel.qqmath.flowset <- 
              function(x, frames, channel,
                       f.value, distribution,
                       ...)
              {
                  distribution <- if (is.function(distribution)) 
                      distribution
                  else if (is.character(distribution)) 
                      get(distribution)
                  else eval(distribution)

                  nobs <- max(unlist(eapply(frames, function(x) nrow(exprs(x)))))

                  qrange <-
                      if (is.null(f.value)) 
                          ppoints(sqrt(nobs))
                      else if (is.numeric(f.value))
                            f.value
                      else f.value(nobs)
                  qrange <- range(qrange)

                  xx <- distribution(qrange)
                  yy <-
                      unlist(eapply(frames,
                                    function(x) {
                                        quantile(evalInFlowFrame(channel, x),
                                                 qrange,
                                                 na.rm = TRUE)
                                    }))

                  dx <- rep(diff(distribution(c(0.25, 0.75))), length(ls(frames)))
                  dy <-
                      unlist(eapply(frames,
                                    function(x) {
                                        diff(quantile(evalInFlowFrame(channel, x),
                                                      c(0.25, 0.75),
                                                      na.rm = TRUE))
                                    }))

                  list(xlim = range(xx, finite = TRUE),
                       ylim = range(yy, finite = TRUE),
                       dx = dx, dy = dy)
              }

          panel.qqmath.flowset <-
              function(x, 
                       frames, channel,
                       f.value, distribution,
                       grid = FALSE,

                       col = superpose.symbol$col,
                       col.points = col,
                       pch = superpose.symbol$pch,
                       cex = superpose.symbol$cex,
                       alpha = superpose.symbol$alpha,
                       col.line = col,
                       lty = superpose.line$lty,
                       lwd = superpose.line$lwd,
                       ...)
              {
                  if (grid) panel.grid(h = -1, v = -1)
                  superpose.symbol <- trellis.par.get("superpose.symbol")
                  superpose.line <- trellis.par.get("superpose.line")

                  nx <- length(x)
                  col.points <- rep(col.points, length = nx)
                  col.line <- rep(col.line, length = nx)
                  pch <- rep(pch, length = nx)
                  cex <- rep(cex, length = nx)
                  lty <- rep(lty, length = nx)
                  lwd <- rep(lwd, length = nx)
                  alpha <- rep(alpha, length = nx)

                  distribution <- if (is.function(distribution)) 
                      distribution
                  else if (is.character(distribution)) 
                      get(distribution)
                  else eval(distribution)

                  nobs <- max(unlist(eapply(frames, function(x) nrow(exprs(x)))))
                  qq <-
                      if (is.null(f.value)) 
                          ppoints(sqrt(nobs))
                      else if (is.numeric(f.value))
                          f.value
                      else f.value(nobs)
                  xx <- distribution(qq)
                  x <- as.character(x)
                  for (i in seq_along(x))
                  {
                      nm <- x[i]
                      yy <-
                          quantile(evalInFlowFrame(channel, frames[[nm]]),
                                   qq, na.rm = TRUE)
                      panel.xyplot(xx, yy,
                                   col.line = col.line[i],
                                   col = col[i],
                                   cex = cex[i],
                                   pch = pch[i],
                                   lty = lty[i],
                                   lwd = lwd[i],
                                   alpha = alpha[i],
                                   ...)
                  }
              }

          if (missing(xlab)) xlab <- deparse(substitute(distribution))
          if (missing(ylab)) ylab <- channel.name
          qqmath(x, data = pd,
                 f.value = f.value, distribution = distribution,
                 prepanel = prepanel.qqmath.flowset,
                 panel = panel.qqmath.flowset,
                 frames = data@frames,
                 channel = channel,
                 xlab = xlab,
                 ylab = ylab,
                      
                 ...)
      })




