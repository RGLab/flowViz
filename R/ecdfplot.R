


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



