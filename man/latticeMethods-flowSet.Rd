\name{lattice-methods}
\docType{methods}
\alias{densityplot,formula,flowSet-method}
\alias{xyplot,formula,flowSet-method}
\alias{qqmath,formula,flowSet-method}
\alias{levelplot,formula,flowSet-method}
\alias{splom,flowFrame,missing-method}
\alias{parallel,flowFrame,missing-method}
\alias{xyplot,flowFrame,missing-method}
\alias{xyplot,formula,flowFrame-method}

\title{ Methods implementing Lattice displays for flow data }

%% things documented here:

% setMethod("densityplot",
%           signature(x = "formula", data = "flowSet"),
% setMethod("xyplot",
%           signature(x = "formula", data = "flowSet"),
% setMethod("qqmath",
%           signature(x = "formula", data = "flowSet"),
% setMethod("levelplot",
%           signature(x = "formula", data = "flowSet"),

% setMethod("splom",
%           signature(x = "flowFrame", data = "missing"),
% setMethod("parallel",
%           signature(x = "flowFrame", data = "missing"),
% setMethod("xyplot",
%           signature(x = "flowFrame", data = "missing"),
% setMethod("xyplot",
%           signature(x = "formula", data = "flowFrame"),



\description{
  Various methods implementing multipanel visualizations for flow data
  using infrastructure provided in the lattice package.  The original
  generics for these methods are defined in lattice, and these S4
  methods (mostly) dispatch on a formula and the \code{data} argument
  which must be of class \code{flowSet} or \code{flowFrame}.  The formula
  has to be fairly basic: conditioning can be done using phenodata
  variables and channel names (the \code{colnames} slot) can be used as
  panel variables.  In the case of the \code{densityplot} method, a
  phenodata variable must be used on the left hand side of the formula
  as a panel variable to stack the densities.  See examples below for
  sample usage.
}

\usage{

## methods for 'flowSet' objects

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


## methods for 'flowFrame' objects

\S4method{splom}{flowFrame,missing}(x, data, 
      pscales = 0,
      smooth = TRUE, pch = ".",
      time = "Time", exclude.time = TRUE,
      panel = function(x, y, smooth, ...) {
          if (smooth) panel.smoothScatter(x, y, ...)
          else panel.xyplot(x, y, ...)
      },
      \dots)


\S4method{parallel}{flowFrame,missing}(x, data, 
         reorder.by = function(x) var(x, na.rm = TRUE),
         time = "Time", exclude.time = TRUE,
         \dots)

\S4method{xyplot}{flowFrame,missing}(x, data,
       xlab = time, ylab = "", time = "Time",
       layout,
       type = "l", smooth = FALSE,
       \dots)

\S4method{xyplot}{formula,flowFrame}(x, data,
       smooth = TRUE, 
       panel = if (smooth) panel.smoothScatter else panel.xyplot,
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

  \item{as.table, contour, labels, pscales, layout, type}{ These
    arguments are passed unchanged to the corresponding methods in
    lattice, and are listed here only because they provide different
    defaults.  See documentation for the original methods for details.
  }
  \item{time}{ A character string giving the name of the column
    recording time. }
  \item{exclude.time}{ logical, specifying whether to exclude the time
    variable from a scatter plot matrix or parallel coordinates plot.
    It is rarely meaningful not to do so. }
  \item{reorder.by}{ a function, which is applied to each column.  The
    columns are ordered by the results.  Reordering can be suppressed by
    setting this to \code{NULL}.  }
  \item{panel}{ panel function used. }
  \item{\dots}{ more arguments, usually passed on to the underlying
    lattice methods. }
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

    \item{splom}{\code{signature(x = "flowFrame", data = "missing")}:
      draws a scatter plot matrix of all channels (excluding time, by
      default) of a \code{flowFrame} object.
    }

    \item{parallel}{\code{signature(x = "flowFrame", data = "missing")}:
      draws a parallel coordinates plot of all channels (excluding time,
      by default) of a \code{flowFrame} object.  This is rarely useful
      without transparency, but that is currently only possible with the
      \code{\link{pdf}} device (and perhaps the aqua device as well).
    }

    \item{xyplot}{\code{signature(x = "flowFrame", data = "missing")}:
      produces diagnostic time series plots of all channels against
      time.
    }

    \item{xyplot}{\code{signature(x = "formula", data = "flowFrame")}:
      more general scatter plots from a \code{flowFrame} object.  The
      formula can be fairly general.
    }
  }
}

\details{
  Transformations of raw data (such as taking logarithms) are not
  supported, except for the \code{xyplot("formula", "flowFrame")}
  method.  Not all standard lattice arguments will have the intended
  effect, but many should.  For a fuller description of possible
  arguments and their effects, consult documentation on lattice (Trellis
  docs would also work for the fundamentals). }


\examples{

data(GvHD)
densityplot(factor(Visit) ~ `FSC-H` | factor(Patient), GvHD)

qqmath( ~ `FSC-H` | factor(Patient), GvHD,
       grid = TRUE, type = "l",
       f.value = ppoints(100))


phenoData(GvHD)$Patient <- factor(phenoData(GvHD)$Patient)
phenoData(GvHD)$Visit <- factor(phenoData(GvHD)$Visit)


## contourplot of bivariate density:

require(colorspace)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
colori <- colorRampPalette(YlOrBr)
levelplot(`SSC-H` ~ `FSC-H` | Visit + Patient, GvHD, n = 20,
          col.regions = colori(50), main = "Contour Plot")



## simple bivariate scatter plot (a.k.a. dot plot)
## by default ('smooth=TRUE') panel.smoothScatter is used

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD, 
       layout = c(7, 5))

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit,
       data = transform("SSC-H"=asinh,"FSC-H"=asinh) \%on\% GvHD,
       layout = c(7, 5))

## version without smooth

xyplot(`SSC-H` ~ `FSC-H` | Visit + Patient, data = GvHD, 
       smooth = FALSE, strip = strip.custom(strip.names = TRUE))

## several examples with time on the X axis

xyplot(`FSC-H` ~ Time | Visit, GvHD, 
       smooth = FALSE, type = "l", 
       subset = (Patient == 5))

xyplot(`FSC-H` ~ Time | Patient:Visit, GvHD, 
       smooth = FALSE, type = "l",
       strip = FALSE, strip.left = TRUE,
       aspect = "xy")


## combine plots for two channels

ssc.time <- 

    xyplot(`SSC-H` ~ Time | factor(Patient):factor(Visit), GvHD, 
           smooth = FALSE, type = "l",
           strip = FALSE,
           strip.left = strip.custom(horizontal = TRUE),
           par.strip.text = list(lines = 3),
           between = list(y = rep(c(0, 0.5), c(6, 1))),
           scales = list(x = list(axs = "i"), y = list(draw = FALSE)),
           layout = c(1, 35))


fsc.time <- 

    xyplot(`FSC-H` ~ Time | factor(Patient):factor(Visit), GvHD, 
           smooth = FALSE, type = "l",
           strip = FALSE,
           strip.left = strip.custom(horizontal = TRUE),
           par.strip.text = list(lines = 3),
           between = list(y = rep(c(0, 0.5), c(6, 1))),
           scales = list(x = list(axs = "i"), y = list(draw = FALSE)),
           layout = c(1, 35))


plot(fsc.time, split = c(1, 1, 2, 1))
plot(ssc.time, split = c(2, 1, 2, 1), newpage = FALSE)


## saving plots as variables allows more manipulation

plot(update(fsc.time[29:35], layout = c(1, 7)),
     split = c(1, 1, 1, 2))

plot(update(ssc.time[29:35], layout = c(1, 7)),
     split = c(1, 2, 1, 2), newpage = FALSE)


## scatter plot matrix of individual flowFrames

splom(GvHD[["s10a07"]], smooth = FALSE)

splom(GvHD[["s10a07"]], smooth = FALSE,
      prepanel.limits = function(x) quantile(as.numeric(x), c(0.0, 0.99)))


## parallel coordinate plots

parallel(GvHD[["s6a01"]])

\dontrun{

## try with PDF device
parallel(GvHD[["s7a01"]], alpha = 0.01)

}

## time plot of individual flowFrames

xyplot(GvHD[["s5a07"]], time = "Time")

xyplot(`FSC-H` + `SSC-H` + `FL1-H` ~ Time, GvHD[["s5a07"]],
       aspect = "xy", outer = TRUE,
       scales = list(y = "free"),
       type = "l", smooth = FALSE)

xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a01"]], nbin = 100)

}

\keyword{methods}
\keyword{dplot}
