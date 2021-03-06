% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contour.R
\docType{methods}
\name{contour-methods}
\alias{contour-methods}
\alias{contour}
\alias{contour,ANY-method}
\alias{contour,flowFrame-method}
\alias{contour,flowSet-method}
\title{Contour plots for flow data}
\usage{
\S4method{contour}{flowFrame}(x, y = 1:2, nlevels = 10, bw,
  grid.size = c(65, 65), add = FALSE, xlab, ylab, xlim, ylim,
  lwd = 1, lty = 1, col = par("fg"), fill = "transparent", ...)
}
\arguments{
\item{x}{An object of class
\code{\link[flowCore:flowFrame-class]{flowFrame}} or
\code{\link[flowCore:flowFrame-class]{flowSet}}.}

\item{y}{Numeric or character vector of length 2 indicating the channels to
plot.}

\item{nlevels}{The approximate number of contour line levels, see
\code{\link[graphics]{contour}} for details.}

\item{bw}{The bandwidth factor used for the kernel density estimation, see
\code{\link[KernSmooth]{bkde2D}} for details.}

\item{grid.size}{The grid size used for the kernel density estimation, see
\code{\link[KernSmooth]{bkde2D}} for details.}

\item{add}{Logical, indicating whether contour lines should be superimposed
on an existing plot.}

\item{xlab, ylab}{The axis annotation.}

\item{xlim, ylim}{The plotting ranges.}

\item{lwd, lty, col, fill}{The usual plotting parameters, i.e. the line width,
line type, line color and fill color. When using a fill color you should
consider alpha blending to improve the results.}

\item{\dots}{Parameters that are passed on to the plotting functions.}
}
\description{
Basic contour plots for both
\code{\link[flowCore:flowFrame-class]{flowFrame}}s and
\code{\link[flowCore:flowFrame-class]{flowSet}}s. The densities for the
contours are estimated using the fast kernel density estimation algorithm
\code{\link[KernSmooth]{bkde2D}}.
}
\section{Methods}{
 \describe{

\item{x = "flowFrame"}{ A regular contour plot of the flow data in the
frame. It can be added on top of an existing plot using the \code{add}
argument.}

\item{x = "flowSet"}{ Overlay of contours of densities for each individual
frame in the set. You should consider using differnt colors and alpha
blending to improve the result. This is only useful for a very limited
number of frames in a set (~5), for larger sets you should consider a
panelled lattice-type plot. Not that \code{bw}, \code{gridSize} and
\code{nlevels} are passed on via the \dots{} argument.}

}
}

\examples{

library(flowCore)
data(GvHD)

## simple contour plot
contour(GvHD[[1]])

## overlay with existing plot
plot(GvHD[[1]], c("FSC-H", "SSC-H"))
contour(GvHD[[1]], add=TRUE, col="lightgray", lty=3)

## colored contours
contour(GvHD[[1]], fill="red")
cols <- rainbow(3, alpha=0.1)
contour(GvHD[[1]], fill=cols, col=cols)

## overlay of multiple flowFrames in a flowSet
contour(GvHD[1:3], col=cols, fill=cols)

}
\seealso{
\code{\link[KernSmooth]{bkde2D}}, \code{\link[graphics]{contour}},
\code{\link[flowCore:flowFrame-class]{flowFrame}},
\code{\link[flowCore:flowFrame-class]{flowSet}}
}
\author{
F. Hahne
}
\keyword{methods}
