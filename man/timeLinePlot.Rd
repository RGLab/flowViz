\name{timeLinePlot}
\alias{timeLinePlot}
\title{flowSet time domain}
\description{
  Plot values of one parameter for each flowFrame in a flowSet against time.
}
\usage{

timeLinePlot(x, par, type=c("stacked", "scaled", "native"), col, ylab=names(x), ...)

}

\arguments{
	\item{x}{An object of class \code{flowSet} containing the data to be plotted.}
	\item{par}{The parameter for which the data is to be plotted}
	\item{type}{One in 'stacked', 'scaled' or 'native'. 'stacked'
	will plot the measurements for the frames on top of each
	other. 'scaled' will align the median values around zero and
	'native' will plot the values in the original dimensions of the
	measurement range.}
	\item{col}{Optional color parameter}
	\item{ylab}{The axis annotation to add on hte y-axis for stacked
	plots}
	\item{\dots}{Further arguments that are passed on to the base
	plotting functions}
}
\details{
Plotting flow cytometry data against the time domain can help to
identify problems with the fuidics or drifts in the instrument setting
during measurement runs.

This function creates plots for all flowFrames in a flowSet for a given
parameter against time. A barplot legend indicates the deviation from
the median for each sample.
}

\author{F. Hahne}

\keyword{dplot}
\seealso{
\code{\link[flowCore:flowCore-package]{flowCore}}
}
\examples{

data(GvHD)
par(ask=TRUE)
timeLinePlot(GvHD, "FL2-A")

timeLinePlot(GvHD, "FL2-A", type="scaled")

timeLinePlot(GvHD, "FL2-A", type="native")
}
