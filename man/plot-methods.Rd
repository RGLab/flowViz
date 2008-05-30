\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,flowFrame,ANY-method}
\alias{plot,flowFrame,character-method}
\alias{plot,flowFrame,missing-method}

\title{ Very basic plotting of flowFrames}

\description{ A basic method to plot \code{\link[flowCore:flowFrame-class]{flowFrame}}
  objects. Depending on the number of dimensions, different types of
  plots are generates. See below for details.  }

\section{Methods}{
  
  \describe{

    \item{x = "flowFrame", y = "ANY"}{ We decide on the number of parameters
      in the  \code{\link[flowCore:flowFrame-class]{flowFrame}}which plot type to use.}
      
    \item{x = "flowFrame", y = "missing"}{ see above }
    
    \item{x = "flowFrame", y = "character"}{ The parameters to plot are
      given as a second argument in the form of a character vector. }
    
  }
}

\details{

  Basic plots for \code{flowFrame} objects. If the object has only a
  single parameter this produces a
  \code{\link[graphics:hist]{histogram}}. For exactly two parameters we
  plot a bivariate density map (see
  \code{\link[geneplotter:smoothscatter]{smoothScatter}}) and for more
  than two parameters we produce a simple \code{\link[lattice]{splom}}
  plot. To select specific parameters from a \code{flowFrame} for
  plotting, either subset the object or specify the parameters as a
  character vector in the second argument to \code{plot}. The
  \code{smooth} parameters lets you toggle between density-type
  \code{\link[geneplotter:smoothscatter]{smoothScatter}} plots and
  regular scatter- or pairs plots.  For far more sophisticated plotting
  of flow cytometry data, see the lattice-style plot methods provided by
  this package.  }

\author{ F. Hahne }

\seealso{
  
  \code{\link{xyplot}},
  \code{\link[flowCore:flowFrame-class]{flowFrame}}, \code{\link{densityplot}}
  
}

\keyword{methods}

