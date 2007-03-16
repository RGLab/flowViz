\name{lattice-methods}
\docType{methods}
\alias{densityplot,formula,flowSet-method}
\alias{xyplot,formula,flowSet-method}
\alias{qqmath,formula,flowSet-method}
\alias{levelplot,formula,flowSet-method}
\title{ Methods implementing Lattice displays for flow data }

\description{
  Various methods implementing multipanel visualizations for flow data
  using infrastructure provided in the lattice package.  The original
  generics for these methods are defined in lattice, and these S4
  methods dispatch on a formula and the \code{data} argument which must
  be of class \code{flowSet}.  The formula has to be fairly basic:
  conditioning can be done using phenodata variables and channel names
  (the \code{colnames} slot) can be used as panel variables.  In the
  case of the \code{densityplot} method, a phenodata variable must be
  used on the left hand side of the formula as a panel variable to stack
  the densities.  See examples below for sample usage.
}

\usage{

\S4method{densityplot}{formula,flowSet}(x, data, 
            xlab, as.table = TRUE, overlap = 0.3, 
            \dots)

\S4method{xyplot}{formula,flowSet}(x, data,
       xlab, ylab, as.table = TRUE,
       pch = ".", smooth = TRUE, \dots)

\S4method{qqmath}{formula,flowSet}(x, data, xlab, ylab,
       f.value = function(n) ppoints(ceiling(sqrt(n))),
       distribution = qnorm,
       \dots)

\S4method{levelplot}{formula,flowSet}(x, data, xlab, ylab,
          as.table = TRUE, contour = TRUE, labels = FALSE,
          n = 50, 
          \dots)

}

\arguments{
  \item{x}{ a formula describing the structure of the plot and the
    variables used to create it }
  \item{data}{ a \code{flowSet} object that serves as a source of data }
  \item{xlab, ylab}{ Labels for data axes, with suitable defaults taken
    from the formula }
  \item{overlap}{ the amount of overlap between stacked density plots }
  \item{pch}{ the plotting character used when \code{smooth=FALSE}}
  \item{smooth}{ logical.  If \code{TRUE},
    \code{\link[geneplotter:panel.smoothScatter]{panel.smoothScatter}}
    is used to display a partially smoothed version of the data.
    Otherwise, points are plotted as in a standard scatter plot.
  }
  \item{f.value, distribution}{ number of points used in Q-Q plot, and
    the reference distribution used.  See
    \code{\link[lattice:qqmath]{qqmath}} for details.  }
  \item{n}{ the number of bins on each axis to be used when evaluating
    the density }
  \item{as.table, contour, labels}{ These arguments are passed unchanged
    to the corresponding methods in lattice, and are listed here only
    because they provide different defaults.  See documentation for the
    original methods for details.
  }


}


\section{Methods}{
  \describe{

    \item{densityplot}{\code{signature(x = "formula", data =
	"flowSet")}: creates density plots for a given channel, with
      samples stacked according to a phenodata variable.  Colors are
      used to indicate common values of this covariate across panels.
    }

    \item{xyplot}{\code{signature(x = "formula", data = "flowSet")}:
      creates scatter plots (a.k.a. dot plots) of a pair of channels
    }

    \item{qqmath}{\code{signature(x = "formula", data = "flowSet")}:
      creates theoretical quantile plots of a given channel, with one or
      more samples per panel
    }

    \item{levelplot}{\code{signature(x = "formula", data = "flowSet")}:
      similar to the \code{xyplot} method, but plots estimated density
      (using \code{\link[MASS:kde2d]{kde2d}}) with a common z-scale and
      an optional color key.
    }

  }
}

\details{
  Transformations of raw data (such as taking logarithms) are not
  supported.  Not all standard lattice arguments will have the intended
  effect, but many should.
}


\examples{

data(GvHD)
densityplot(factor(Visit) ~ `FSC-H` | factor(Patient), GvHD)

qqmath( ~ `FSC-H` | factor(Patient), GvHD,
       grid = TRUE, type = "l",
       f.value = ppoints(100))


phenoData(GvHD)$Patient <- factor(phenoData(GvHD)$Patient)
phenoData(GvHD)$Visit <- factor(phenoData(GvHD)$Visit)

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD, 
       layout = c(7, 5))

xyplot(`SSC-H` ~ `FSC-H` | Visit + Patient, data = GvHD, 
       smooth = FALSE, strip = strip.custom(strip.names = TRUE))


require(colorspace)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
colori <- colorRampPalette(YlOrBr, space = "Lab")
levelplot(`SSC-H` ~ `FSC-H` | Visit + Patient, GvHD, n = 20,
          col.regions = colori(50), main = "Contour Plot")

}

\keyword{methods}
\keyword{dplot}
