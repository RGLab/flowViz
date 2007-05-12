\name{ecdfplot}
\docType{methods}
\alias{ecdfplot}
\alias{ecdfplot,formula,flowSet-method}

\title{ Method implementing Lattice ECDF plots for flow data }

\description{

  This function creates Trellis displays of Empirical Cumulative
  Distribution Functions from flow cytometry data using a formula
  interface.

}

\usage{

## methods for 'flowSet' objects

\S4method{ecdfplot}{formula,flowSet}(x, data, xlab, ylab,
         f.value = function(n) ppoints(ceiling(sqrt(n))),
         type = "l", 
         \dots)

}

\arguments{
  \item{x}{ a formula describing the structure of the plot and the
    variables to be used in the display.   }
  \item{data}{ a \code{flowSet} object that serves as a source of data }
  \item{xlab, ylab}{ Labels for data axes, with suitable defaults taken
    from the formula }
  \item{f.value}{ determines the number of points used in the plot
    \code{\link[latticeExtra:ecdfplot]{ecdfplot}} for details.  }
  \item{type}{ type of rendering; by default lines are drawn }
%   \item{ref}{ logical; whether to add reference lines at 0 and 1 }
  \item{\dots}{ more arguments, usually passed on to the underlying
    lattice methods and the panel function. }
}


\section{Methods}{
  \describe{

    \item{ecdfplot}{\code{signature(x = "formula", data = "flowSet")}:
      plote empirical CDF for a given channel, with one or more samples
      per panel
    }
  }
}

\seealso{
  Not all standard lattice arguments will have the intended
  effect, but many should.  For a fuller description of possible
  arguments and their effects, consult documentation on lattice.
}


\examples{

data(GvHD)

ecdfplot(~ `FSC-H` | Patient, GvHD, f.value = ppoints(100))

ecdfplot(~ asinh(`FSC-H`) | Patient, GvHD,
         strip = strip.custom(strip.names = TRUE),
         ref = FALSE)

ecdfplot(~ asinh(`FSC-H`) | Patient, GvHD, groups = Visit,
         strip = strip.custom(strip.names = TRUE),
         ref = FALSE, auto.key = list(columns = 4))

}

\keyword{methods}
\keyword{dplot}
