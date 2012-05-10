setMethod("contour",
          signature("flowFrame"),
          function(x, y=1:2, nlevels=10, bw, grid.size=c(65,65), add=FALSE,
                   xlab, ylab, xlim, ylim, lwd=1, lty=1, col=par("fg"),
                   fill="transparent", ...)
      {
          if(length(y) != 2)
              stop("You must specify two dimensions!")
          if(is.character(y))
              y <- match(y[1:2], colnames(x))
          if (missing(bw))
              bw <- diff(apply(exprs(x[,y]), 2, quantile, probs=c(0.05, 0.95),
                               na.rm=TRUE)) / 25
          cnx <- colnames(x)[y[1]]
          if(missing(xlab))
              xlab <- cnx
          cny <- colnames(x)[y[2]]
          if(missing(ylab))
              ylab <- cny
          plotType("contour", c(cnx, cny))
          exp <- exprs(x[,y])
          ## compute the bivariate density and the contour lines
          xr <- range(exp[,1], na.rm=TRUE)
          yr <- range(exp[,2], na.rm=TRUE)
          range <- list(xr+c(-1,1)*bw[1]*2.5, yr+c(-1,1)*bw[2]*2.5)
          z <- bkde2D(exp, bw, grid.size, range.x=range)
          ll <- contourLines(z$x1, z$x2, z$fhat, nlevels=nlevels)
          ## plot everything as polygons
          if(missing(xlim))
              xlim <- unlist(range(x[,y[1]]))
           if(missing(ylim))
              ylim <- unlist(range(x[,y[2]]))
          if(!add) 
              plot(z$x1, z$x2, type="n", xlab=xlab, ylab=ylab,
                   xlim=xlim, ylim=ylim, ...)
          for(ct in ll) {
              polygon(ct$x, ct$y, border=col, col=fill, lwd=lwd, lty=lty)
          }
      })
	
	
setMethod("contour",
          signature("flowSet"),
          function(x, y=1:2, add=FALSE, xlab, ylab, lwd=1, lty=1,
                   col=par("fg"), fill="transparent", ...)
      {
          localPlot <- function(..., bw, nlevels, grid.size)
              plot(...)
          if(length(y) != 2)
              stop("You must specify two dimensions!")
          if(is.character(y))
              y <- match(y[1:2], colnames(x))
           if(missing(xlab))
              xlab <- colnames(x)[y[1]]
          if(missing(ylab))
              ylab <- colnames(x)[y[2]]
          if(!add) {
              mins <- apply(fsApply(x,each_col,min),2,min)
              maxs <- apply(fsApply(x,each_col,max),2,max)
              localPlot(c(mins[y[1]], maxs[y[1]]),
                        c(mins[y[2]], maxs[y[2]]), type="n",
                        xlab=xlab, ylab=ylab, ...)		
          }
          fill <- rep(fill, length.out=length(x))
          col <- rep(col, length.out=length(x))
          lwd <- rep(lwd, length.out=length(x))
          lty <- rep(lty, length.out=length(x))
          for(i in 1:length(x))
              contour(x[[i]], y, col=col[i], fill=fill[i], lwd=lwd[i],
                      lty=lty[i], add=TRUE, ...)
      })
