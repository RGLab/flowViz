
setMethod("densityplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab,
                   as.table = TRUE, overlap = 0.3, 
                   ...)
      {
          pd <- phenoData(data)@data
          ## ugly hack to suppress warnings about coercion introducing
          ## NAs (needs to be `undone' inside prepanel and panel
          ## functions):
          pd$name <- factor(pd$name) 
          channel <- x[[3]]  
          if (length(channel) == 3)
          {
              channel <- channel[[2]]
              x[[3]][[2]] <- as.name("name")
          }
          else x[[3]] <- as.name("name")
          channel <- as.character(channel)

          prepanel.densityplot.flowset <- 
              function(x, y, darg = list(n = 30),
                       frames, channel,
                       overlap = 0.3,
                       ...)
              {
                  xl <- numeric(0)
                  for (nm in as.character(x))
                  {
                      xx <- exprs(frames[[nm]])[, channel]
                      xl <- c(xl, range(xx))
                  }
                  list(xlim = range(xl, finite = TRUE))
              }

          panel.densityplot.flowset <-
              function(x, y, darg = list(n = 30), ref = FALSE,
                       frames, channel,
                       overlap = 0.3,

                       col = superpose.polygon$col,
                       lty = superpose.polygon$lty,
                       lwd = superpose.polygon$lwd,
                       alpha = superpose.polygon$alpha,
                       border = superpose.polygon$border,
                       ...)
              {
                  superpose.line <- trellis.par.get("superpose.line")
                  superpose.polygon <- trellis.par.get("superpose.polygon")
                  reference.line <- trellis.par.get("reference.line")
                  ycode <- as.numeric(y)
                  if (any(duplicated(ycode)))
                  {
                      warning("Some combinations seem to have multiple samples.  \n  Only one will be used.")
                  }
                  ny <- nlevels(y)
                  col <- rep(col, length = ny)
                  lty <- rep(lty, length = ny)
                  lwd <- rep(lwd, length = ny)
                  alpha <- rep(alpha, length = ny)
                  border <- rep(border, length = ny)
                  x <- as.character(x)
                  height <- (1 + overlap)
                  for (i in rev(seq_len(ny)))
                      if (i %in% ycode)
                      {
                          nm <- x[match(i, ycode)]
                          xx <- exprs(frames[[nm]])[, channel]
                          h <- do.call(density, c(list(x = xx), darg))
                          n <- length(h$x)
                          max.d <- max(h$y)
                          panel.polygon(x = h$x[c(1, 1:n, n)],
                                        y = i + height * c(0, h$y, 0) / max.d,
                                        col = col[i], border = border[i],
                                        lty = lty[i], lwd = lwd[i], alpha = alpha[i])
                          if (ref)
                          {
                              panel.abline(h = i,
                                           col = reference.line$col,
                                           lty = reference.line$lty,
                                           lwd = reference.line$lwd,
                                           alpha = reference.line$alpha)
                          }
                      }
              }

          if (missing(xlab)) xlab <- channel
          xyplot(x, data = pd, 

                 prepanel = prepanel.densityplot.flowset,
                 panel = panel.densityplot.flowset,

                 frames = data@frames,
                 channel = channel,
                 as.table = as.table,
                 overlap = overlap,

                 xlab = xlab,

                 lattice.options = list(axis.padding = list(factor = c(0.6, 1 + 2 * overlap))),
                 horizontal = TRUE, ...)
      })








setMethod("xyplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab, ylab,
                   as.table = TRUE,
                   pch = ".", smooth = TRUE,
                   ...)
      {
          pd <- phenoData(data)@data
          ## ugly hack to suppress warnings about coercion introducing
          ## NAs (needs to be `undone' inside prepanel and panel
          ## functions):
          pd$name <- factor(pd$name) 
          channel.y <- x[[2]]
          channel.x <- x[[3]]
          if (length(channel.x) == 3)
          {
              channel.x <- channel.x[[2]]
              x[[3]][[2]] <- as.name("name")
              x[[2]] <- NULL
          }
          else
          {
              x[[3]] <- as.name("name")
              x[[2]] <- NULL
          }
          channel.x <- as.character(channel.x)
          channel.y <- as.character(channel.y)

          prepanel.xyplot.flowset <- 
              function(x, 
                       frames, channel.x, channel.y,
                       ...)
              {
                  if (length(x) > 1) stop("must have only one flow frame per panel")
                  xl <- numeric(0)
                  yl <- numeric(0)
                  for (nm in as.character(x))
                  {
                      xx <- exprs(frames[[nm]])[, channel.x]
                      yy <- exprs(frames[[nm]])[, channel.y]
                      xl <- c(xl, range(xx))
                      yl <- c(yl, range(yy))
                  }
                  list(xlim = range(xl, finite = TRUE),
                       ylim = range(yl, finite = TRUE))
              }

          panel.xyplot.flowset <-
              function(x, 
                       frames, channel.x, channel.y,

                       pch, smooth,
                       ...)
              {
                  x <- as.character(x)
                  if (length(x) > 1) stop("must have only one flow frame per panel")
                  
                  for (nm in x)
                  {
                      xx <- exprs(frames[[nm]])[, channel.x]
                      yy <- exprs(frames[[nm]])[, channel.y]
                      
                      if (smooth) panel.smoothScatter(xx, yy, ...)
                      else panel.xyplot(xx, yy, pch = pch, ...)
                  }
              }

          if (missing(xlab)) xlab <- channel.x
          if (missing(ylab)) ylab <- channel.y
          densityplot(x, data = pd, 

                      prepanel = prepanel.xyplot.flowset,
                      panel = panel.xyplot.flowset,

                      frames = data@frames,
                      channel.x = channel.x,
                      channel.y = channel.y,
                      as.table = as.table,

                      xlab = xlab,
                      ylab = ylab,
                      pch = pch, smooth = smooth,

                      ...)
      })




setMethod("levelplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab, ylab,
                   as.table = TRUE, contour = TRUE, labels = FALSE,
                   n = 50, 
                   ...)
      {
          pd <- phenoData(data)@data
          ## ugly hack to suppress warnings about coercion introducing
          ## NAs (needs to be `undone' inside prepanel and panel
          ## functions):
          pd$name <- factor(pd$name) 
          channel.y <- x[[2]]
          channel.x <- x[[3]]
          if (length(channel.x) == 3)
          {
              channel.x <- channel.x[[2]]
              cond <- paste("|", paste(as.character(x[[3]][[3]]), collapse = " "))
          }
          else cond <- ""
          channel.x <- as.character(channel.x)
          channel.y <- as.character(channel.y)

          computeKde <- function(nm)
          {
              xx <- exprs(data@frames[[nm]])[, channel.x]
              yy <- exprs(data@frames[[nm]])[, channel.y]
              temp <- kde2d(xx, yy, n = n)
              d <- data.frame(x = rep(temp$x, n),
                              y = rep(temp$y, each = n),
                              z = as.vector(temp$z))
          }
          ft <- do.call(make.groups, sapply(as.character(pd$name), computeKde, simplify = FALSE))

          rownames(pd) <- pd$name 
          ft <- cbind(ft, pd[as.character(ft$which), , drop = FALSE])

          if (missing(xlab)) xlab <- channel.x
          if (missing(ylab)) ylab <- channel.y

          my.formula <- as.formula(paste("z ~ x * y", cond))

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





### old version with formulas like ~x and ~x | a

## setMethod("densityplot",
##           signature(x = "formula", data = "flowSet"),
##           function(x, data, xlab,
##                    stack = FALSE, overlap = 0.3,
##                    ...)
##       {
##           pd <- phenoData(data)@data
##           ## ugly hack to suppress warnings about coercion introducing
##           ## NAs (needs to be `undone' inside prepanel and panel
##           ## functions):
##           pd$name <- factor(pd$name) 
##           channel <- x[[2]]
##           if (length(channel) == 3)
##           {
##               channel <- channel[[2]]
##               x[[2]][[2]] <- as.name("name")
##           }
##           else x[[2]] <- as.name("name")
##           channel <- as.character(channel)

##           prepanel.densityplot.flowset <- 
##               function(x, darg = list(n = 30),
##                        frames, channel,
##                        stack = FALSE, overlap = 0.3,
##                        ...)
##               {
##                   xl <- numeric(0)
##                   yl <- if (stack) c(1, length(x) + 1 + overlap) else 0
##                   dxl <- if (stack) 1 else numeric(0)
##                   dyl <- if (stack) 1 else numeric(0)
                      
##                   for (nm in as.character(x))
##                   {
##                       xx <- exprs(frames[[nm]])[, channel]
##                       if (!stack)
##                       {
##                           h <- do.call(density, c(list(x = xx), darg))
##                           xl <- c(xl, range(h$x))
##                           yl <- c(yl, max(h$y))
##                           quants <- quantile(xx, prob = c(0.15, 0.85), 
##                                              names = FALSE, na.rm = TRUE)
##                           ok <- h$x > quants[1] & h$x < quants[2]
##                           dxl <- c(dxl, diff(h$x[ok]))
##                           dyl <- c(dyl, diff(h$y[ok]))
##                       }
##                       else
##                       {
##                           xl <- c(xl, range(xx))
##                       }
##                   }
##                   list(xlim = range(xl, finite = TRUE),
##                        ylim = range(yl, finite = TRUE),
##                        dx = dxl, dy = dyl)
##               }

##           panel.densityplot.flowset <-
##               function(x, darg = list(n = 30), ref = FALSE,
##                        frames, channel,
##                        stack = FALSE, overlap = 0.3,

##                        col = if (stack) plot.polygon$col else superpose.line$col,
##                        lty = if (stack) plot.polygon$lty else superpose.line$lty,
##                        lwd = if (stack) plot.polygon$lwd else superpose.line$lwd,
##                        alpha = if (stack) plot.polygon$alpha else superpose.line$alpha,
##                        border = plot.polygon$border,
##                        ...)
##               {
##                   superpose.line <- trellis.par.get("superpose.line")
##                   plot.polygon <- trellis.par.get("plot.polygon")
##                   reference.line <- trellis.par.get("reference.line")
##                   nx <- length(x)

##                   if (stack)
##                   {
##                       height <- (1 + overlap)
##                       infolist <-
##                           lapply(as.character(x),
##                                  function(nm) {
##                                      xx <- exprs(frames[[nm]])[, channel]
##                                      h <- do.call(density, c(list(x = xx), darg))
##                                      list(med = median(xx, na.rm = TRUE),
##                                           dens = h)
##                                  })
##                       ord.med <- order(sapply(infolist, "[[", "med"))
##                       for (i in rev(seq_along(ord.med)))
##                       {
##                           h <- infolist[[ord.med[i]]][["dens"]]
##                           n <- length(h$x)
##                           max.d <- max(h$y)
##                           panel.polygon(x = h$x[c(1, 1:n, n)],
##                                         y = i + height * c(0, h$y, 0) / max.d,
##                                         col = col, border = border,
##                                         lty = lty, lwd = lwd, alpha = alpha)
##                           if (ref)
##                           {
##                               panel.abline(h = i,
##                                            col = reference.line$col,
##                                            lty = reference.line$lty,
##                                            lwd = reference.line$lwd,
##                                            alpha = reference.line$alpha)
##                           }
##                       }
##                   }
##                   else
##                   {
##                       col <- rep(col, length = nx)
##                       lty <- rep(lty, length = nx)
##                       lwd <- rep(lwd, length = nx)
##                       alpha <- rep(alpha, length = nx)

##                       if (ref)
##                       {
##                           panel.abline(h = 0,
##                                        col = reference.line$col,
##                                        lty = reference.line$lty,
##                                        lwd = reference.line$lwd,
##                                        alpha = reference.line$alpha)
##                       }

##                       x <- as.character(x)
##                       for (i in seq_along(x))
##                       {
##                           nm <- x[i]
##                           xx <- exprs(frames[[nm]])[, channel]
##                           h <- do.call(density, c(list(x = xx), darg))
##                           llines(h,
##                                  col = col[i], lty = lty[i],
##                                  lwd = lwd[i], alpha = alpha[i],
##                                  ...)
##                       }
##                   }
##               }

##           if (missing(xlab)) xlab <- channel
##           default.scales <-
##               if (stack) list(y = list(draw = FALSE))
##               else list()
##           densityplot(x, data = pd, 

##                       prepanel = prepanel.densityplot.flowset,
##                       panel = panel.densityplot.flowset,

##                       frames = data@frames,
##                       channel = channel,
##                       stack = stack, overlap = overlap,

##                       xlab = xlab,
##                       default.scales = default.scales,
                      
##                       ...)
##       })

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

                 prepanel = prepanel.qqmath.flowset,
                 panel = panel.qqmath.flowset,

                 frames = data@frames,
                 channel = channel,

                 ylab = ylab,
                      
                 ...)
      })

