\name{splom}

\docType{methods}

\alias{splom}
\alias{splom,flowFrame,missing-method}
\alias{panel.splom.flowframe}


\title{ Method implementing Lattice scatter plot matrices for flow data. }


\description{
  
  This function create Trellis scatter plots matrices (splom) from flow
  cytometry data.

}

\usage{

\S4method{splom}{flowFrame,missing}(
    x,
    data, 
    pscales,
    time,
    exclude.time=TRUE,
    names=FALSE,
    \dots)

panel.splom.flowframe(
    x,
    frame,
    \dots)

}




\arguments{

  \item{x}{ A formula describing the structure of the plot and the
    variables to be used in the display. }

  
  \item{data, frame}{ A \code{\link[flowCore:flowFrame-class]{flowFrame}}
    object that serves as the source of data. }

    
  \item{pscales}{ This arguments is passed unchanged to the
    corresponding methods in lattice, and is listed here only because it
    provides a different default.  See documentation for the original
    methods for details. }

  
  \item{time}{ A character string giving the name of the data column
    recording time. If not provided, we try to guess from the available
    parameters. }

  
  \item{exclude.time}{ Logical, specifying whether to exclude the time
    variable from a scatter plot matrix. Defaults to \code{TRUE}. }

  
  \item{names}{ Logical specifying wether gate names should be added to
    the plot. Currently, this feature is not supported for splom plots. }

  
  \item{\dots}{ More arguments, usually passed on to the underlying
    lattice methods. }

}


\section{Methods}{
  \describe{

    \item{splom}{\code{signature(x = "flowFrame", data = "missing")}:
      Creates a scatter plot matrix for a whole \code{flowFrame}. The
      actual plotting is done by \code{\link{panel.xyplot.flowframe}},
      so most additional arguments from the regular \code{xyplot} are
      also valid here. }

  }
}


\details{
  
  The function draws a scatter plot matrix of the data for each flow
  parameter in a \code{flowFrame}. For the most, one can think about
  this as a rectangular arrangement of separate
  \code{\link[xyplot]{xyplots}}, and most of that functionality is also
  available here. To be more precise, the function repeatedly calls
  \code{panel.xyplot.flowframe} to do the actual plotting. Please see
  its documentation for details.
  
}

\seealso{
  
  Not all standard lattice arguments will have the intended
  effect, but many should.  For a fuller description of possible
  arguments and their effects, consult documentation on lattice.

}


\author{ F. Hahne, D. Sarkar }

\examples{

data(GvHD)

tf <- transformList(colnames(GvHD)[3:7], asinh)
dat <- tf \%on\% GvHD[[3]]


## scatter plot matrix of individual flowFrames
lattice.options(panel.error=NULL)
splom(dat)

splom(dat[,1:3], smooth = FALSE)


## displaying filters
rg <- rectangleGate("FSC-H"=c(200,400), "SSC-H"=c(300,700),
"FL1-H"=c(2,4), "FL2-A"=c(4,7))
splom(dat, filter=rg)

splom(dat, filter=rectangleGate("FSC-H"=c(400,800)))

splom(dat[,1:4], smooth = FALSE, filter=norm2Filter("FSC-H", "SSC-H", scale=1.5))


}

\keyword{methods}
\keyword{dplot}
