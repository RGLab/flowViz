

setMethod("densityplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab, ...)
      {
          pd <- phenoData(data)@data
          ## ugly hack to suppress warnings about coercion introducing NAs
          pd$name <- factor(pd$name) 
          channel <- x[[2]]
          if (length(channel) == 3)
          {
              channel <- channel[[2]]
              x[[2]][[2]] <- as.name("name")
          }
          else x[[2]] <- as.name("name")
          channel <- as.character(channel)

          prepanel.densityplot.flowset <- 
              function(x, darg = list(n = 30), frames, channel, ...)
              {
                  xl <- numeric(0)
                  yl <- 0
                  dxl <- numeric(0)
                  dyl <- numeric(0)
                  for (nm in as.character(x))
                  {
                      xx <- exprs(frames[[nm]])[, channel]
                      h <- do.call(density, c(list(x = xx), darg))
                      xl <- c(xl, range(h$x))
                      yl <- c(yl, range(h$y))
                      quants <- quantile(xx, prob = c(0.15, 0.85), 
                                         names = FALSE, na.rm = TRUE)
                      ok <- h$x > quants[1] & h$x < quants[2]
                      dxl <- c(dxl, diff(h$x[ok]))
                      dyl <- c(dyl, diff(h$y[ok]))
                  }
                  list(xlim = range(xl, finite = TRUE),
                       ylim = range(yl, finite = TRUE),
                       dx = dxl, dy = dyl)
              }

          panel.densityplot.flowset <-
              function(x, darg = list(n = 30), ref = FALSE,
                       frames, channel,
                       col = superpose.line$col,
                       lty = superpose.line$lty,
                       lwd = superpose.line$lwd,
                       alpha = superpose.line$alpha,
                       ...)
              {
                  superpose.line <- trellis.par.get("superpose.line")
                  nx <- length(x)
                  col <- rep(col, length = nx)
                  lty <- rep(lty, length = nx)
                  lwd <- rep(lwd, length = nx)
                  alpha <- rep(alpha, length = nx)

                  if (ref)
                  {
                      reference.line <- trellis.par.get("reference.line")
                      panel.abline(h = 0,
                                   col = reference.line$col,
                                   lty = reference.line$lty,
                                   lwd = reference.line$lwd,
                                   alpha = reference.line$alpha)
                  }

                  x <- as.character(x)
                  for (i in seq_along(x))
                  {
                      nm <- x[i]
                      xx <- exprs(frames[[nm]])[, channel]
                      h <- do.call(density, c(list(x = xx), darg))
                      llines(h,
                             col = col[i], lty = lty[i],
                             lwd = lwd[i], alpha = alpha[i],
                             ...)
                  }
              }

          if (missing(xlab)) xlab <- channel
          densityplot(x, data = pd, 

                      prepanel = prepanel.densityplot.flowset,
                      panel = panel.densityplot.flowset,

                      frames = data@frames,
                      channel = channel,

                      xlab = xlab,
                      
                      ...)
      })

## QQ plot

## setGeneric("qqmath")



setMethod("qqmath",
          signature(x = "formula", data = "flowSet"),
          function(x, data, ylab,
                   f.value = function(n) ppoints(ceiling(sqrt(n))),
                   distribution = qnorm,
                   ...)
      {
          pd <- phenoData(data)@data
          ## ugly hack to suppress warnings about coercion introducing NAs
          pd$name <- factor(pd$name) 
          channel <- x[[2]]
          if (length(channel) == 3)
          {
              channel <- channel[[2]]
              x[[2]][[2]] <- as.name("name")
          }
          else x[[2]] <- as.name("name")
          channel <- as.character(channel)

          prepanel.densityplot.flowset <- 
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
                                        quantile(exprs(x)[, channel],
                                                 qrange,
                                                 na.rm = TRUE)
                                    }))

                  dx <- rep(diff(distribution(c(0.25, 0.75))), length(ls(frames)))
                  dy <-
                      unlist(eapply(frames,
                                    function(x) {
                                        diff(quantile(exprs(x)[, channel],
                                                      c(0.25, 0.75),
                                                      na.rm = TRUE))
                                    }))

                  list(xlim = range(xx, finite = TRUE),
                       ylim = range(yy, finite = TRUE),
                       dx = dx, dy = dy)
              }

          panel.densityplot.flowset <-
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
                          quantile(exprs(frames[[nm]])[, channel],
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

          if (missing(ylab)) ylab <- channel
          qqmath(x, data = pd,

                 f.value = f.value, distribution = distribution,

                 prepanel = prepanel.densityplot.flowset,
                 panel = panel.densityplot.flowset,

                 frames = data@frames,
                 channel = channel,

                 ylab = ylab,
                      
                 ...)
      })



