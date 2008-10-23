






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




