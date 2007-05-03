\name{xyplot}
\docType{methods}
\alias{xyplot}
\alias{splom}
\alias{xyplot,formula,flowSet-method}
\alias{xyplot,flowFrame,missing-method}
\alias{xyplot,formula,flowFrame-method}
\alias{splom,flowFrame,missing-method}

\title{ Method implementing Lattice ECDF plots for flow data }

%% things documented here:

% setMethod("xyplot",
%           signature(x = "formula", data = "flowSet"),
% setMethod("xyplot",
%           signature(x = "flowFrame", data = "missing"),
% setMethod("xyplot",
%           signature(x = "formula", data = "flowFrame"),
% setMethod("splom",
%           signature(x = "flowFrame", data = "missing"),


\description{

  These functions create Trellis scatter plots (a.k.a. dot plots in the
  Flow Cytometry community) from flow cytometry data.

}

\usage{

## methods for 'flowSet' objects

\S4method{xyplot}{formula,flowSet}(x, data,
       xlab, ylab, as.table = TRUE,
       pch = ".", smooth = TRUE,
       filter = NULL,
       filterResults = NULL,
       displayFilter = TRUE,
       \dots)

## methods for 'flowFrame' objects

\S4method{xyplot}{flowFrame,missing}(x, data,
       xlab = time, ylab = "", time = "Time",
       layout,
       type = "l", smooth = FALSE,
       \dots)

\S4method{xyplot}{formula,flowFrame}(x, data,
       smooth = TRUE, 
       panel = if (smooth) panel.smoothScatter else panel.xyplot,
       \dots)

\S4method{splom}{flowFrame,missing}(x, data, 
      pscales = 0,
      smooth = TRUE, pch = ".",
      time = "Time", exclude.time = TRUE,
      panel = function(x, y, smooth, ...) {
          if (smooth) panel.smoothScatter(x, y, ...)
          else panel.xyplot(x, y, ...)
      },
      \dots)


}

\arguments{
  \item{x}{ a formula describing the structure of the plot and the
    variables to be used in the display.   }
  \item{data}{ a \code{flowSet} or \code{flowFrame} object that serves
    as the source of data }
  \item{smooth}{ logical.  If \code{TRUE},
    \code{\link[geneplotter:panel.smoothScatter]{panel.smoothScatter}}
    is used to display a partially smoothed version of the data.
    Otherwise, points are plotted as in a standard scatter plot.
  }
  \item{xlab, ylab}{ Labels for data axes, with suitable defaults taken
    from the formula }
  \item{pch}{ the plotting character used when \code{smooth=FALSE}}
  \item{type}{ type of rendering; see
    \code{\link[lattice:panel.xyplot]{panel.xyplot}} for details. }

  \item{panel}{ panel function used. }
  \item{filter}{ A \link[flowCore:filter-class]{filter} object.  Ignored
    when \code{filterResults} is not \code{NULL}. }
  \item{filterResults}{ A
    \link[flowCore:filterResult-class]{filterResult} object, that
    is the result of applying a filter on each frame in \code{data}.  If
    this is \code{NULL} but \code{filter} is not, it will be created on
    the fly by applying \code{filter} to each frame.  }
  \item{displayFilter}{ logical or list.  If \code{TRUE}, an attempt is
    made to display the geometric boundaries of \code{filter} (if
    specified) in the plot.  The same holds if \code{displayFilter} is a
    list, but in that case the components are interpreted as graphical
    parameters to be passed on to
    \code{\link[lattice:llines]{panel.polygon}}
  }
  \item{time}{ A character string giving the name of the column
    recording time. }
  \item{exclude.time}{ logical, specifying whether to exclude the time
    variable from a scatter plot matrix. }
  \item{as.table, pscales, layout}{ These arguments are passed unchanged
    to the corresponding methods in lattice, and are listed here only
    because they provide different defaults.  See documentation for the
    original methods for details.
  }


  \item{\dots}{ more arguments, usually passed on to the underlying
    lattice methods. }
}


\section{Methods}{
  \describe{

    \item{xyplot}{\code{signature(x = "formula", data = "flowSet")}:
      creates scatter plots (a.k.a. dot plots) of a pair of channels
    }

    \item{xyplot}{\code{signature(x = "flowFrame", data = "missing")}:
      produces diagnostic time series plots of all channels against
      time.
    }

    \item{xyplot}{\code{signature(x = "formula", data = "flowFrame")}:
      more general scatter plots from a \code{flowFrame} object.  The
      formula can be fairly general.
    }

    \item{splom}{\code{signature(x = "flowFrame", data = "missing")}:
      draws a scatter plot matrix of all channels (excluding time, by
      default) of a \code{flowFrame} object.
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

## time plot of individual flowFrames

xyplot(GvHD[["s5a07"]], time = "Time")

xyplot(`FSC-H` + `SSC-H` + `FL1-H` ~ Time, GvHD[["s5a07"]],
       aspect = "xy", outer = TRUE,
       scales = list(y = "free"),
       type = "l", smooth = FALSE)

xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a01"]], nbin = 100)

## displaying filters

n2gate <- norm2Filter("SSC-H", "FSC-H")

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD,
       filter = n2gate, ## displayFilter = TRUE,
       layout = c(7, 5))

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit,
       data = transform("SSC-H"=asinh,"FSC-H"=asinh) \%on\% GvHD,
       smooth = FALSE,
       filter = n2gate, ## displayFilter = list(border = "black"),
       layout = c(7, 5))

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit,
       data = transform("SSC-H"=asinh,"FSC-H"=asinh) \%on\% GvHD,
       subset = Patient == "6",
       filter = rectangleGate("FSC-H"=c(5.5, 7),"SSC-H"=c(4, 6.5)),
       displayFilter = list(lwd = 3))

n2gate.results <- filter(GvHD, n2gate \%on\% transform("SSC-H"=asinh,"FSC-H"=asinh))

xyplot(`SSC-H` ~ `FSC-H` | Visit, data = GvHD,
       subset = Patient == "6",
       filterResults = n2gate.results, smooth = FALSE)

xyplot(`SSC-H` ~ `FSC-H` | Visit,
       data = GvHD,
       subset = Patient == "6",
       smooth = FALSE,
       filter = n2gate,
       filterResults = n2gate.results,
       displayFilter = list(lwd = 0.5))


xyplot(`SSC-H` ~ `FSC-H` | Visit,
       data = transform("SSC-H"=asinh,"FSC-H"=asinh) \%on\% GvHD,
       subset = Patient == "6",
       smooth = FALSE,
       filter = rectangleGate("FSC-H"=c(5.5, 7),"SSC-H"=c(4, 6.5)),
       filterResults = n2gate.results,
       displayFilter = list(border = 'red'))


}

\keyword{methods}
\keyword{dplot}
