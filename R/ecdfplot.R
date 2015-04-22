#' @export 
#' @rdname ecdfplot
prepanel.ecdfplot.flowset <- 
    function(x, frames, channel,
             f.value,
             ...)
{
    distribution <- qunif
    nobs <- max(unlist(eapply(frames, function(x) nrow(exprs(x)))))

    qrange <-
        if (is.null(f.value)) 
            ppoints(sqrt(nobs))
        else if (is.numeric(f.value))
            f.value
        else f.value(nobs)
    qrange <- range(qrange)

    yy <- distribution(qrange)
    xx <-
        unlist(eapply(frames,
                      function(x) {
                          quantile(evalInFlowFrame(channel, x),
                                   qrange,
                                   na.rm = TRUE)
                      }))

    dy <- rep(diff(distribution(c(0.25, 0.75))), length(ls(frames)))
    dx <-
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



#' @export 
#' @rdname ecdfplot
panel.ecdfplot.flowset <-
    function(x, 
             frames, channel,
             f.value, 
             ref = TRUE,

             groups = NULL,
             subscripts,

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
    if (ref)
    {
        reference.line <- trellis.par.get("reference.line")
        do.call(panel.abline, c(list(h = c(0, 1)), reference.line))
    }
    superpose.symbol <- trellis.par.get("superpose.symbol")
    superpose.line <- trellis.par.get("superpose.line")

    if (is.null(groups))
    {
        nx <- length(x)
        col.points <- rep(col.points, length = nx)
        col.line <- rep(col.line, length = nx)
        pch <- rep(pch, length = nx)
        cex <- rep(cex, length = nx)
        lty <- rep(lty, length = nx)
        lwd <- rep(lwd, length = nx)
        alpha <- rep(alpha, length = nx)
    }
    else
    {
        groups <- as.factor(groups)[subscripts]
        stopifnot(length(groups) == length(x))
        ## goal: make colors etc vectors as before, but
        ## associate by group

        ng <- nlevels(groups)
        gcode <- as.numeric(groups)
        col.points <- rep(col.points, length = ng)[gcode]
        col.line <- rep(col.line, length = ng)[gcode]
        pch <- rep(pch, length = ng)[gcode]
        cex <- rep(cex, length = ng)[gcode]
        lty <- rep(lty, length = ng)[gcode]
        lwd <- rep(lwd, length = ng)[gcode]
        alpha <- rep(alpha, length = ng)[gcode]
    }

    distribution <- qunif
    nobs <- max(unlist(eapply(frames, function(x) nrow(exprs(x)))))
    qq <-
        if (is.null(f.value)) 
            ppoints(sqrt(nobs))
        else if (is.numeric(f.value))
            f.value
        else f.value(nobs)
    yy <- distribution(qq)
    x <- as.character(x)
    for (i in seq_along(x))
    {
        nm <- x[i]
        xx <-
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




#' Method implementing Lattice ECDF plots for flow data
#' 
#' 
#' This function creates Trellis displays of Empirical Cumulative Distribution
#' Functions from flow cytometry data using a formula interface.
#' 
#' 
#' @name ecdfplot
#' @aliases ecdfplot ecdfplot,formula,flowSet-method panel.ecdfplot.flowset
#' prepanel.ecdfplot.flowset
#' @docType methods
#' @param x a formula describing the structure of the plot and the variables to
#' be used in the display.  For the prepanel and panel functions, a vector of
#' names for the flow frames to be used in the panel.
#' @param data a \code{flowSet} object that serves as a source of data
#' @param xlab Labels for data axes, with suitable defaults taken from the
#' formula
#' @param f.value determines the number of points used in the plot
#' \code{\link[latticeExtra:ecdfplot]{ecdfplot}} for details.
#' @param panel,prepanel the panel and prepanel functions.
#' @param type type of rendering; by default lines are drawn
#' @param as.table logical; whether to draw panels from top left
#' @param ref logical; whether to add reference lines at 0 and 1
#' @param frames environment containing frame-specific data
#' @param channel expression involving names of columns in the data
#' @param groups,subscripts grouping variable, if specified, and subscripts
#' indexing which frames are being used in the panel.  See
#' \code{\link[lattice:xyplot]{xyplot}} for details.
#' @param col,col.points,pch,cex,alpha,col.line,lty,lwd vector of graphical
#' parameters that are replicated for each group
#' @param \dots more arguments, usually passed on to the underlying lattice
#' methods and the panel function.
#' @section Methods: \describe{
#' 
#' \item{ecdfplot}{\code{signature(x = "formula", data = "flowSet")}: plote
#' empirical CDF for a given channel, with one or more samples per panel } }
#' @seealso Not all standard lattice arguments will have the intended effect,
#' but many should.  For a fuller description of possible arguments and their
#' effects, consult documentation on lattice.
#' @keywords methods dplot
#' @examples
#' 
#' 
#' data(GvHD)
#' 
#' ecdfplot(~ `FSC-H` | Patient, GvHD, f.value = ppoints(100))
#' 
#' ecdfplot(~ asinh(`FSC-H`) | Patient, GvHD,
#'          strip = strip.custom(strip.names = TRUE),
#'          ref = FALSE)
#' 
#' ecdfplot(~ asinh(`FSC-H`) | Patient, GvHD, groups = Visit,
#'          strip = strip.custom(strip.names = TRUE),
#'          ref = FALSE, auto.key = list(columns = 4))
#' 
#' 
#' @export 
setMethod("ecdfplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab,
                   f.value = function(n) ppoints(ceiling(sqrt(n))),
                   prepanel = prepanel.ecdfplot.flowset,
                   panel = panel.ecdfplot.flowset,
                   type = "l", as.table = TRUE,
                   ...)
      {
          ocall <- sys.call(sys.parent())
          ccall <- match.call(expand.dots = FALSE)
          ccall <- manipulate.call(ocall, ccall)
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
          if (missing(xlab)) xlab <- channel.name
          ccall$x <- x
          ccall$data <- pd
          ccall$f.value <- f.value
          ccall$prepanel <- prepanel
          ccall$panel <- panel
          ccall$type <- type
          ccall$as.table <- as.table
          ccall$xlab <- xlab
          ccall$frames <- data@frames
          ccall$channel <- channel
          ccall[[1]] <- quote(latticeExtra::ecdfplot)
          ans <- eval.parent(ccall)
          ans$call <- ocall
          ans
      })



