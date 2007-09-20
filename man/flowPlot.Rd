\name{flowPlot}
\alias{flowPlot}
\alias{flowPlot,flowFrame-method}
\docType{methods}
\title{Standard Plots for Flow Cytometry Data
}
\description{
  A method that makes standard plots from a \code{flowFrame}. The user may also
  provide various \code{filter} or \code{filterResult} arguments to customize
  the plot.
}
\usage{

flowPlot(x, ...)

\S4method{flowPlot}{flowFrame}(x, child, filter = NULL,
         plotParameters = c("FSC-H","SSC-H"),
         logx = FALSE, logy = FALSE,
         parent, colParent="Grey", colChild="Blue",
         showFilter = TRUE, gate.fill = "transparent",
         gate.border = "black", xlab, ylab, xlim, ylim,
         \dots)
}

\arguments{
	\item{x}{An object of class \code{flowFrame} that contains the data to be plotted.}
	\item{child}{An optional argument of class \code{filterResult} that specifies
	a subset of the data that are included in the \code{filterResult}}
	\item{filter}{}
	\item{plotParameters}{}
	\item{logx}{}
	\item{logy}{}
	\item{parent}{}
	\item{colParent}{}
	\item{colChild}{}
	\item{showFilter}{}
	\item{gate.fill}{}
	\item{gate.border}{}
	\item{xlab}{}
	\item{ylab}{}
	\item{xlim}{}
	\item{ylim}{}
	\item{\dots}{}
}
\details{
  The plot that is most commonly used in flow cytometry data analysis is usuall called
  a "dot plot". In common statistical language, we would call this a scatter plot. The
  basic idea is a 2-dimensional plot that shows the location of every cell in regard
  to the measurements made on it, for example, forward scatter vs side scatter. Most
  applications will, in addition to the data, want to show information about one or
  more filters (gates). Since there can be a very large number of cells in a sample,
  it is common to show a smoothed version of the data that doesn't involve registering
  every point on the graph.

}
\author{P. Haaland}

\keyword{methods}
\seealso{
\code{\link[flowCore:flowCore-package]{flowCore}}
}
\examples{

data(GvHD)
flowPlot(GvHD[["s5a01"]])
flowPlot(transform("SSC-H"=asinh,"FSC-H"=asinh) \%on\% GvHD[["s5a01"]])

}
