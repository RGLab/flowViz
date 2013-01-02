\name{timeLinePlot-methods}
\docType{methods}
\alias{timeLinePlot}
\alias{timeLinePlot-methods}
\alias{timeLinePlot,ANY,missing-method}
\alias{timeLinePlot,flowFrame,character-method}
\alias{timeLinePlot,flowSet,character-method}

\title{ Plot channel values against time }

\description{ Plots values of one parameter for each flowFrame in a
  flowSet against time.  }

\section{Methods}{

  General usage:
     
  \code{timeLinePlot(x, channel, type=c("stacked", "scaled",
    "native"), col, ylab=names(x), binSize, varCut=1, ...)}
  
  \describe{
    
    \item{x}{An object of class \code{\link[flowCore:flowFrame-class]{flowFrame}} or
      \code{\link[flowCore:flowSet-class]{flowSet}} containing the data to be plotted.}
	
    \item{channel}{The parameter for which the data is to be plotted}
	
    \item{type}{One in 'stacked', 'scaled' or 'native'. 'stacked' will
      plot the measurements for the frames on top of each
      other. 'scaled' will align the median values around zero and
      'native' will plot the values in the original dimensions of the
      measurement range.}
      
    \item{col}{Optional color parameter.}
  
    \item{ylab}{The axis annotation to add on the y-axis for stacked
      plots.}
	
    \item{binSize}{The number of events per bin. If not set, a reasonable
      default is computed.}
    
    \item{varCut}{The cutoff in the adjusted variance to which the
      quality score is computed. Basically, all values that are outside
      of the confidence interval defined by \code{\[my - signma *
      varCut, my + sigma * varCut\]} will contribute to a positive
      quality score value.}
 
    \item{\dots}{Further arguments that are passed on to the base plotting
      functions.}
  
  }
  
  \describe{

    \item{x = "ANY", channel = "missing"}{ This casts a useful error
      message about missing inputs. }

    \item{x = "flowFrame", channel = "character"}{ For a
      \code{\link[flowCore:flowFrame-class]{flowFrame}} we only plot a single parameter
      against time. }
    
    \item{x = "flowSet", channel = "character"}{ For a
      \code{\link[flowCore:flowSet-class]{flowSet}}, we plot the lines for a  single parameter
      against time for each \code{\link[flowCore:flowFrame-class]{flowFrame}} in the set. }
  
  }
}

\value{A numeric vector of quality scores.}
      
\details{

  Plotting flow cytometry data against the time domain can help to
  identify problems with the fluidics or drifts in the instrument
  setting during measurement runs.

  This function creates plots for all flowFrames in a flowSet for a
  given parameter against time. A barplot legend indicates the deviation
  from the median for each sample. There is also a flowFrame method,
  which will create a plot for a single flowFrame only.

  In addition, the function computes a quality score for each frame,
  which essentially is the sum of the positive distances of each bin
  mean from a frame-specific confidence interval, divided by the number
  of bins. Values larger than zero indicate a problem.

}

\author{F. Hahne}

\keyword{dplot}
\keyword{methods}

\seealso{
  
  \code{\link[flowCore:flowFrame-class]{flowFrame}}, \code{\link[flowCore:flowSet-class]{flowSet}}
  
}

\examples{

data(GvHD)
opar <- par(ask=TRUE)

res <- timeLinePlot(GvHD[[1]], "SSC-H")
res

res <- timeLinePlot(GvHD, "SSC-H")

res <- timeLinePlot(GvHD, "SSC-H", type="scaled", varCut=4)

res <- timeLinePlot(GvHD[1:4], "SSC-H", type="native", binSize=50)

par(opar)

}


