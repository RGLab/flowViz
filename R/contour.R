#' Contour plots for flow data
#' 
#' Basic contour plots for both
#' \code{\link[flowCore:flowFrame-class]{flowFrame}}s and
#' \code{\link[flowCore:flowFrame-class]{flowSet}}s. The densities for the
#' contours are estimated using the fast kernel density estimation algorithm
#' \code{\link[KernSmooth]{bkde2D}}.
#' 
#' 
#' @name contour-methods
#' @aliases contour contour-methods contour,ANY-method contour,flowFrame-method
#' contour,flowSet-method
#' @docType methods
#' @param x An object of class
#' \code{\link[flowCore:flowFrame-class]{flowFrame}} or
#' \code{\link[flowCore:flowFrame-class]{flowSet}}.
#' @param y Numeric or character vector of length 2 indicating the channels to
#' plot.
#' @param nlevels The approximate number of contour line levels, see
#' \code{\link[graphics]{contour}} for details.
#' @param bw The bandwidth factor used for the kernel density estimation, see
#' \code{\link[KernSmooth]{bkde2D}} for details.
#' @param grid.size The grid size used for the kernel density estimation, see
#' \code{\link[KernSmooth]{bkde2D}} for details.
#' @param add Logical, indicating whether contour lines should be superimposed
#' on an existing plot.
#' @param xlab,ylab The axis annotation.
#' @param xlim,ylim The plotting ranges.
#' @param lwd,lty,col,fill The usual plotting parameters, i.e. the line width,
#' line type, line color and fill color. When using a fill color you should
#' consider alpha blending to improve the results.
#' @param \dots Parameters that are passed on to the plotting functions.
#' @section Methods: \describe{
#' 
#' \item{x = "flowFrame"}{ A regular contour plot of the flow data in the
#' frame. It can be added on top of an existing plot using the \code{add}
#' argument.}
#' 
#' \item{x = "flowSet"}{ Overlay of contours of densities for each individual
#' frame in the set. You should consider using differnt colors and alpha
#' blending to improve the result. This is only useful for a very limited
#' number of frames in a set (~5), for larger sets you should consider a
#' panelled lattice-type plot. Not that \code{bw}, \code{gridSize} and
#' \code{nlevels} are passed on via the \dots{} argument.}
#' 
#' }
#' @author F. Hahne
#' @seealso
#' 
#' \code{\link[KernSmooth]{bkde2D}}, \code{\link[graphics]{contour}},
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowFrame-class]{flowSet}}
#' @keywords methods
#' @examples
#' 
#' library(flowCore)
#' data(GvHD)
#' 
#' ## simple contour plot
#' contour(GvHD[[1]])
#' 
#' ## overlay with existing plot
#' plot(GvHD[[1]], c("FSC-H", "SSC-H"))
#' contour(GvHD[[1]], add=TRUE, col="lightgray", lty=3)
#' 
#' ## colored contours
#' contour(GvHD[[1]], fill="red")
#' cols <- rainbow(3, alpha=0.1)
#' contour(GvHD[[1]], fill=cols, col=cols)
#' 
#' ## overlay of multiple flowFrames in a flowSet
#' contour(GvHD[1:3], col=cols, fill=cols)
#' 
#' @importFrom KernSmooth bkde2D dpik
#' @export 
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
