\name{lattice-methods}
\docType{methods}
\alias{qqmath}
\alias{levelplot}
\alias{parallel}
\alias{qqmath,formula,flowSet-method}
\alias{levelplot,formula,flowSet-method}
\alias{parallel,flowFrame,missing-method}
\alias{parallel,formula,flowSet-method}

\title{ Methods implementing Lattice displays for flow data }

%% things documented here:
% setMethod("qqmath",
%           signature(x = "formula", data = "flowSet"),
% setMethod("levelplot",
%           signature(x = "formula", data = "flowSet"),

% setMethod("parallel",
%           signature(x = "flowFrame", data = "missing"),


\description{

  Various methods implementing multipanel visualizations for flow data
  using infrastructure provided in the lattice package.  The original
  generics for these methods are defined in lattice, and these S4
  methods (mostly) dispatch on a formula and the \code{data} argument
  which must be of class \code{flowSet} or \code{flowFrame}.  The
  formula has to be fairly basic: conditioning can be done using
  phenodata variables and channel names (the \code{colnames} slot) can
  be used as panel variables. See examples below for sample usage.

}

\usage{

## methods for 'flowSet' objects
\S4method{qqmath}{formula,flowSet}(
    x,
    data,
    xlab,
    ylab,
    f.value = function(n) ppoints(ceiling(sqrt(n))),
    distribution = qnorm,
    \dots)


\S4method{levelplot}{formula,flowSet}(
    x,
    data,
    xlab,
    ylab,
    as.table = TRUE,
    contour = TRUE,
    labels = FALSE,
    n = 50,
    \dots)


## methods for 'flowFrame' objects
\S4method{parallel}{flowFrame,missing}(
    x,
    data, 
    reorder.by = function(x) var(x, na.rm = TRUE),
    time = "Time",
    exclude.time = TRUE,
    \dots)

}

\arguments{
  
  \item{x}{ a formula describing the structure of the plot and the
    variables to be used in the display. }
  
  \item{data}{ a \code{flowSet} object that serves as a source of data. }
  
  \item{xlab, ylab}{ Labels for data axes, with suitable defaults taken
    from the formula }
  
  \item{f.value, distribution}{ number of points used in Q-Q plot, and
    the reference distribution used.  See
    \code{\link[lattice:qqmath]{qqmath}} for details.  }
  
  \item{n}{ the number of bins on each axis to be used when evaluating
    the density }

  \item{as.table, contour, labels}{ These arguments are passed unchanged
    to the corresponding methods in lattice, and are listed here only
    because they provide different defaults.  See documentation for the
    original methods for details.}
  
  \item{time}{ A character string giving the name of the column
    recording time. }
  
  \item{exclude.time}{ logical, specifying whether to exclude the time
    variable from a scatter plot matrix or parallel coordinates plot.
    It is rarely meaningful not to do so. }
  
  \item{reorder.by}{ a function, which is applied to each column.  The
    columns are ordered by the results.  Reordering can be suppressed by
    setting this to \code{NULL}.  }
  
  \item{\dots}{ more arguments, usually passed on to the underlying
    lattice methods. }
}


\section{Methods}{
  \describe{

    \item{qqmath}{\code{signature(x = "formula", data = "flowSet")}:
      creates theoretical quantile plots of a given channel, with one or
      more samples per panel
    }

    \item{levelplot}{\code{signature(x = "formula", data = "flowSet")}:
      similar to the \code{xyplot} method, but plots estimated density
      (using \code{\link[MASS:kde2d]{kde2d}}) with a common z-scale and
      an optional color key.
    }

    \item{parallel}{\code{signature(x = "flowFrame", data = "missing")}:
      draws a parallel coordinates plot of all channels (excluding time,
      by default) of a \code{flowFrame} object.  This is rarely useful
      without transparency, but that is currently only possible with the
      \code{\link{pdf}} device (and perhaps the aqua device as well).
    }

  }
}

\details{
  Not all standard lattice arguments will have the intended
  effect, but many should.  For a fuller description of possible
  arguments and their effects, consult documentation on lattice (Trellis
  docs would also work for the fundamentals).
}


\examples{

data(GvHD)

qqmath( ~ `FSC-H` | factor(Patient), GvHD,
       grid = TRUE, type = "l",
       f.value = ppoints(100))


## contourplot of bivariate density:

require(colorspace)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
colori <- colorRampPalette(YlOrBr)
levelplot(asinh(`SSC-H`) ~ asinh(`FSC-H`) | Visit + Patient, GvHD, n = 20,
          col.regions = colori(50), main = "Contour Plot")



## parallel coordinate plots

parallel(GvHD[["s6a01"]])

\dontrun{

## try with PDF device
parallel(GvHD[["s7a01"]], alpha = 0.01)

}

}

\keyword{methods}
\keyword{dplot}
