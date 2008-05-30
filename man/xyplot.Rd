\name{xyplot}
\docType{methods}
\alias{xyplot}
\alias{splom}
\alias{panel.xyplot.flowset}
\alias{prepanel.xyplot.flowset}
\alias{panel.xyplot.flowframe}
\alias{panel.xyplot.flowframe.time}
\alias{prepanel.xyplot.flowframe.time}
\alias{xyplot,formula,flowSet-method}
\alias{xyplot,flowFrame,missing-method}
\alias{xyplot,formula,flowFrame-method}
\alias{splom,flowFrame,missing-method}

\title{ Method implementing Lattice xyplots for flow data }

\description{
  
  These functions create Trellis scatter plots (a.k.a. dot
  plots in the Flow Cytometry community) from flow cytometry data.

}

\usage{

## methods for 'flowSet' objects
\S4method{xyplot}{formula,flowSet}(
    x,
    data,
    xlab,
    ylab,
    as.table=TRUE,
    prepanel=prepanel.xyplot.flowset,
    panel=panel.xyplot.flowset,
    pch = ".",
    smooth = TRUE,
    filter = NULL,
    \dots)


## methods for 'flowFrame' objects
\S4method{xyplot}{flowFrame,missing}(
    x,
    data,
    xlab=time,
    ylab="",
    time="Time",
    layout,
    type="l",
    smooth=FALSE,
    \dots)

\S4method{xyplot}{formula,flowFrame}(
    x,
    data,
    smooth=TRUE, 
    panel=if(smooth) panel.smoothScatter else panel.xyplot,
    \dots)

\S4method{splom}{flowFrame,missing}(
    x,
    data, 
    pscales=0,
    smooth=TRUE,
    pch=".",
    time="Time",
    exclude.time=TRUE,
    panel=function(x, y, smooth, ...) {
         if (smooth) panel.smoothScatter(x, y, ...)
         else panel.xyplot(x, y, ...)
      },
    \dots)


## panel and prepanel functions for flowSets
prepanel.xyplot.flowframe.time(
    x,
    y,
    time.x,
    expr,
    \dots)

panel.xyplot.flowframe.time(
    x, 
    y,
    time.x,
    expr,
    type="l",
    nrpoints=0,
    binSize=100,
    ...)

panel.xyplot.flowframe(
    x,
    y,
    frame,
    filter=NULL,
    smooth=TRUE,
    channel.x.name,
    channel.y.name,
    pch=".",
    gpar,
    ...)

## panel and prepanel functions for flowSets
prepanel.xyplot.flowset(
    x,
    frames,
    channel.x,
    channel.y,
    \dots)

panel.xyplot.flowset(
    x,
    frames,
    channel.x,
    channel.y,
    channel.x.name,
    channel.y.name, 
    filter=NULL,
    pch=".",
    smooth=TRUE,
    gpar,
    \dots)

}

\arguments{

  \item{x}{ a formula describing the structure of the plot and the
    variables to be used in the display.  In
    \code{prepanel.xyplot.flowset} and \code{panel.xyplot.flowset}, also
    the names of \code{\link[flowCore:flowFrame-class]{flowFrame}}s. }
  
  \item{data, y}{ a \code{\link[flowCore:flowSet-class]{flowSet}} or
    \code{\link[flowCore:flowFrame-class]{flowFrame}} object that serves as the source of
    data. }
  
  \item{smooth}{ logical. If \code{TRUE},
    \code{\link[geneplotter:panel.smoothScatter]{panel.smoothScatter}}
    is used to display a partially smoothed version of the data.
    Otherwise, points are plotted as in a standard scatter plot. }

  \item{nrpoints} {the number of points plotted on the smoothed plot in
    sparse regions. This is only listed here because we use a different
    default. See
    \code{\link[geneplotter:panel.smoothScatter]{panel.smoothScatter}}
    for details. }
  
  \item{xlab, ylab}{ Labels for data axes, with suitable defaults taken
    from the formula. }
  
  \item{pch}{ the plotting character used when \code{smooth=FALSE}. }
  
  \item{type}{ type of rendering; see
    \code{\link[lattice:panel.xyplot]{panel.xyplot}} for details. For
    the basic \code{\link[flowCore:flowFrame-class]{flowFrame}} method without a detailed
    formula, the addtional type \code{discrete} is available, which
    plots a smoothed average of the flow cytometry values against
    time. }
  
  \item{binSize} {the size of a bin (i.e., the number of events within a
    bin) used for the smoothed average timeline plots. }

  \item{prepanel}{ the prepanel function.  See
    \code{\link[lattice:xyplot]{xyplot}}. }
  
  \item{panel}{ the panel function.  See
    \code{\link[lattice:xyplot]{xyplot}}. }
  
  \item{filter}{ a \code{\link[flowCore:filter-class]{filter}} or
    \code{\link[flowCore:filterResult-class]{filterResult}} object or a list of such objects
    of the same length as the \code{\link[flowCore:flowSet-class]{flowSet}} to be
    superimposed on the plot. The software will figure out whether the
    \code{\link[flowCore:filter-class]{filter}} needs to be evaluated (in which case
    providing a \code{\link[flowCore:filterResult-class]{filterResult}} can speed up things
    considerably). }
  
  \item{time}{ a character string giving the name of the column
    recording time. If not provided, we try to guess from the parameter
    names. }

  \item{time.x}{ a numeric vector of time values. }

  \item{expr}{ a numeric matrix of flow cytometry values, basically a
    copy of the \code{exprs} slot of a \code{\link[flowCore:flowFrame-class]{flowFrame}}. }
  
  \item{exclude.time}{ logical, specifying whether to exclude the time
    variable from a scatter plot matrix. Defaults to \code{TRUE}. }
  
  \item{as.table, pscales, layout}{ These arguments are passed unchanged
    to the corresponding methods in lattice, and are listed here only
    because they provide different defaults.  See documentation for the
    original methods for details. }

  \item{channel.x, channel.y}{ expressions defining the x and y
    variables in terms of columns in the data.  Can involve functions or
    multiple columns from the data. }
  
  \item{channel.x.name, channel.y.name}{ character strings giving
    corresponding names, used to match filter parameters if
    applicable. }
  
  \item{frame, frames}{an environment containing frame-specific data. }

  \item{gpar}{an list of graphical parameters that are passed down to
    the low level grid functions. See \code{\link[grid:gpar]{gpar}} for details. }

  \item{\dots}{ more arguments, usually passed on to the underlying
    lattice methods. }

}


\section{Methods}{
  \describe{

    \item{xyplot}{\code{signature(x = "formula", data = "flowSet")}:
      creates scatter plots (a.k.a. dot plots) of a pair of channels }

    \item{xyplot}{\code{signature(x = "flowFrame", data = "missing")}:
      produces diagnostic time series plots of all channels against
      time.  }

    \item{xyplot}{\code{signature(x = "formula", data = "flowFrame")}:
      more general scatter plots from a \code{flowFrame} object.  The
      formula can be fairly general.  }

    \item{splom}{\code{signature(x = "flowFrame", data = "missing")}:
      draws a scatter plot matrix of all channels (excluding time, by
      default) of a \code{flowFrame} object.  }

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
       smooth = FALSE, type = "a",
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


## a bivariate scatterplot
xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a01"]], nbin = 100)


## displaying filters
n2gate <- norm2Filter("SSC-H", "FSC-H")

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD,
       filter=n2gate, layout=c(7, 5))

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit,
       data=transform("SSC-H"=asinh,"FSC-H"=asinh) \%on\% GvHD,
       smooth=FALSE, filter=n2gate, layout=c(7, 5))

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit,
       data = transform("SSC-H"=asinh,"FSC-H"=asinh) \%on\% GvHD,
       subset = Patient == "6",
       filter = rectangleGate("FSC-H"=c(5.5, 7),"SSC-H"=c(4, 6.5)))

n2gate.results <- filter(GvHD, n2gate)

xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD,
       subset=Patient == "6",
       filter=n2gate.results, smooth=FALSE)

xyplot(`SSC-H` ~ `FSC-H` | Visit,
       data = GvHD,
       subset = Patient == "6",
       smooth = FALSE,
       filter = n2gate,
       filterResults = n2gate.results,
       displayFilter = list(lwd = 0.5))

}

\keyword{methods}
\keyword{dplot}
