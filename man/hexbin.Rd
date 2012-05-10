\name{hexbin}
\alias{hexbin}
\title{ Fast hexagon binning. }
\description{
 Yet another hexagon binning routine - this one adapted from the
scagnostics package.  This is very preliminary and may not be retained.
}
\usage{
hexbin(x, y, bins = 50)
}
\arguments{
  \item{x}{ x values }
  \item{y}{ y values }
  \item{bins}{ The approximate number of bins, in the x direction. }
}
\details{
 This implementation maps x and y to the unit square and does the binning
 in that region.
}
\value{
 A list, with componentes \code{counts}, \code{xbin} and \code{ybin}.
}
\references{ The scagnostics package.}
\examples{
 x = runif(100)
 y = runif(100)
 hb = hexbin(x,y)
}
