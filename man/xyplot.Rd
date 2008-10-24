\name{xyplot}

\docType{methods}

\alias{panel.xyplot.flowframe.time}
\alias{panel.xyplot.flowframe}
\alias{panel.xyplot.flowset}
\alias{prepanel.xyplot.flowframe.time}
\alias{prepanel.xyplot.flowframe}
\alias{prepanel.xyplot.flowset}
\alias{xyplot,flowFrame,missing-method}
\alias{xyplot,formula,flowFrame-method}
\alias{xyplot,formula,flowSet-method}
\alias{xyplot,formula,gateView-method}
\alias{xyplot,formula,view-method}
\alias{xyplot,view,missing-method}
\alias{xyplot}


\title{ Methods implementing Lattice xyplots for flow data. }


\description{
  
  These functions create Trellis scatter plots (a.k.a. dot
  plots in the Flow Cytometry community) from flow cytometry data.

}

\usage{
## Method for 'flowFrame' objects without a formula.
## This creates plots of all flow parameters agains
## time. 
\S4method{xyplot}{flowFrame,missing}(
    x,
    data,
    time,
    xlab,
    ylab="",
    layout,
    prepanel=prepanel.xyplot.flowframe.time,
    panel=panel.xyplot.flowframe.time,
    type="discrete",
    \dots)


## prepanel function for time line plots of flowFrames
prepanel.xyplot.flowframe.time(
    x,
    y,
    frame,
    time,
    \dots)


## panel function for time line plots of flowFrames
panel.xyplot.flowframe.time(
    x, 
    y,
    frame,
    time,
    type="discrete",
    nrpoints=0,
    binSize=100,
    \dots)



## method for formulae with 'flowFrame' objects
\S4method{xyplot}{formula,flowFrame}(
    x,
    data,
    smooth=TRUE,
    prepanel=prepanel.xyplot.flowframe,
    panel=panel.xyplot.flowframe,
    \dots)


## prepanel function for generic xyplots of flowFrames
prepanel.xyplot.flowframe(
    frame,
    channel.x.name,
    channel.y.name,
    \dots)
   

## panel function for generic xyplots of flowFrames
panel.xyplot.flowframe(
    x,
    y,
    frame,
    filter=NULL,
    smooth=TRUE,
    margin=TRUE,
    outline=FALSE,
    channel.x.name,
    channel.y.name,
    pch=gpar$flow.symbol$pch,
    alpha=gpar$flow.symbol$alpha,
    cex=gpar$flow.symbol$cex,
    col=gpar$flow.symbol$col,
    gp,
    \dots)



## method for 'flowSet' objects
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
    par.settings=NULL,
    \dots)


## prepanel function for generic xyplots of flowSets
prepanel.xyplot.flowset(
    x,
    frames,
    channel.x.name,
    channel.y.name,
    \dots)


## panel function for generic xyplots of flowSets
panel.xyplot.flowset(
    x,
    frames,
    filter=NULL,
    channel.x,
    channel.y,
    \dots)


## method for various workflow objects
\S4method{xyplot}{formula,view}(
    x,
    data,
    \dots)


\S4method{xyplot}{view,missing}(
    x,
    data,
    \dots) 


\S4method{xyplot}{formula,gateView}(
    x,
    data,
    filter=NULL,
    par.settings,
    \dots)
      
}

\arguments{

  \item{x}{ A formula describing the structure of the plot and the
    variables to be used in the display.  In the \code{prepanel} and
    \code{panel} functions, also the names of
    \code{\link[flowCore:flowFrame-class]{flowFrames}} or any of the
    annotation data columns in the \code{phenoData} slot. }

  
  \item{data, y, frame}{ a \code{flowSet} or
    \code{\link[flowCore:flowFrame-class]{flowFrame}} object that serves
    as the source of data. For the workflow methods, this can also be
    various \code{\link[flowCore]{view}} or
    \code{\link[flowCore]{actionItem}} objects.}

  
  \item{time}{ A character string giving the name of the data column
    recording time. If not provided, we try to guess from the available
    parameters. }

  
  \item{xlab, ylab}{ Labels for data axes, with suitable defaults taken
    from the formula. }

  
  \item{as.table, layout}{ These arguments are passed unchanged
    to the corresponding methods in lattice, and are listed here only
    because they provide different defaults.  See documentation for the
    original methods for details. }

  
  \item{type}{ type of rendering; see
    \code{\link[lattice:panel.xyplot]{panel.xyplot}} for details. For
    the basic \code{flowFrame} method without a detailed formula, the
    addtional type \code{discrete} is available, which plots a smoothed
    average of the flow cytometry values against time. }

  
  \item{nrpoints} {The number of points plotted on the smoothed plot in
    sparse regions. This is only listed here because we use a different
    default. See
    \code{\link[geneplotter:panel.smoothScatter]{panel.smoothScatter}}
    for details. }

  
  \item{binSize} {The size of a bin (i.e., the number of events within a
    bin) used for the smoothed average timeline plots. }

  
  \item{channel.x.name, channel.y.name}{ Character strings giving
    corresponding names used to match filter parameters if applicable. }

  
  \item{smooth}{ Logical. If \code{TRUE}, \code{panel.smoothScatter} is
    used to display a partially smoothed version of the data.
    Otherwise, events are plotted individually, as in a standard scatter
    plot. }

  
  \item{filter}{ A \code{\link[flowCore:filter-class]{filter}},
    \code{\link[flowCore:filterResult-class]{filterResult}} or
    \code{\link[flowCore:filterResult-class]{filterResultList}} object
    or a list of such objects of the same length as the
    \code{flowSet}. The appropriate spherical 2D representation of this
    filter will be superimposed on the plot if \code{smooth=TRUE}, or
    the resukt of the filtering operation will be indicated by grouping
    if \code{smooth=FALSE} . The software will figure out whether the
    \code{filter} needs to be evaluated in order to be plotted (in which
    case providing a \code{filterResult} can speed things up
    considerably). }

  
  \item{margin}{Logical indicating whether to truncate the density
    estimation on the margins of the measurement range and plot margin
    events as lines if \code{smooth=TRUE}. To avoid visual artifacts it
    is highly recommended to set this option to \code{TRUE}. }

  
  \item{outline}{ Logical, specifying whether to add the boundaries of a
    gate to the plot when \code{smooth=FALSE} in addition to the
    grouping. Defaults to \code{FALSE}. }

  
  \item{pch, cex, col, alpha}{ Graphical parameters used when
    \code{smooth=FALSE}. These mostly exist for conveniance and much
    more control is available throught the \code{lattice}-like
    \code{par.setting} and \code{flowViz.par.set} customization. See
    \code{\link{flowViz.par.set}} for details. }

  
  \item{par.settings}{ A list of lists of graphical parameters.  See
    \code{\link{flowViz.par.set}} for details.}

  
  \item{gp}{A list of graphical parameters that are passed down to the
    low level panel functions. This is for internal use only. The public
    user interface to set graphical parameters is either
    \code{par.settings} for customization of a single call or
    \code{flowViz.par.set} for customization of session-wide defaults. }

  
  \item{prepanel}{ The prepanel function. See
    \code{\link[lattice:xyplot]{xyplot}}. }

  
  \item{panel}{ The panel function. See
    \code{\link[lattice:xyplot]{xyplot}}. }

  
  \item{channel.x, channel.y}{ Expressions defining the x and y
    variables in terms of columns in the data.  Can involve functions or
    multiple columns from the data, however this usage is discouraged. }


  \item{frames}{An environment containing frame-specific data. }

  
  \item{\dots}{ More arguments, usually passed on to the underlying
    lattice methods. }
  
}


\section{Methods}{

  \describe{

    \item{xyplot}{\code{signature(x = "flowFrame", data = "missing")}:
      Creates diagnostic time series plots of flow parameter values
      against time. These plots are useful to detect quality issues in
      the raw data. If not provided explicitely via the \code{tine}
      argument, the time parameter will be automatically detected. The
      additional arguments \code{xlab}, \code{ylab}, \code{nrpoints},
      and \code{layout} are only listed because \code{flowViz} provides
      different defaults. Internally, they are directly passed on to the
      underlying \code{lattice} functions}. Argument \code{type} can be
      a combination of any of the types allowed in \code{lattice
      xyplots}, or \code{discrete}, in which case a smoothed average of
      the parameter against time is plotted. \code{binSize} controls the
      binning that is used for the smoothing procedure.

    
    \item{xyplot}{\code{signature(x = "formula", data = "flowFrame")}:
      Creates scatter plots (a.k.a. dot plots) of a pair of FCM
      channels. Depending on the setting of the \code{smooth} argument,
      the data will be rendered as a partially smoothed density estimate
      (\code{smooth=TRUE}, the default) or as a regular scatter plot
      with separate points for individual events. The formula interface
      allows for fairly general plotting, however there are certain
      limitations on the use of expressions as part of the
      formulae. Unless you are sure about what you are doing, you should
      transform the raw data in a separate step using one of the tools
      in the \code{\link[flowCore-package]{flowCore}} package rather
      than inline using the formula interface. The method allows to
      superimpose gating results though the \code{filter} argument. If
      \code{smooth=TRUE}, we try to add spherical 2D representations of
      the gates if applicable. For \code{smooth=FALSE}, gates are
      indicated by a grouping mechanism using different point shapes or
      colors (unless \code{outline} is also \code{TRUE}, in which case
      the gate outlines are superimposed in addition to the
      grouping). Argument \code{margins} controls how events on the
      margins of the measurement range are treated. The default
      (\code{TRUE}) is to discard them from any density estimation and
      later add them as separate glyphs. See
      \code{\link{flowViz.par.set}} for details on controlling graphical
      parameters in these plots. }


    \item{xyplot}{\code{signature(x = "formula", data = "flowSet")}:
      Scatter plots from a \code{flowSet} object. We allow for
      conditioning on variables in the \code{phenoData} slot of the
      \code{flowSet}. All additional arguments that apply to the
      \code{flowFrame} method are also valid for \code{flowSets}. }


    \item{xyplot}{\code{signature(x = "formula", data = "view")}:
      Scatter plots from a \code{view} object. Depending on the
      particulars of the view, the method tries to come up with
      reasonable defaults. Full customization is also available though
      all the before mentioned arguments. }


    \item{xyplot}{\code{signature(x = "view", data = "missing")}:
      Scatter plots from a \code{view} object. This is the fallback
      method in the case of absolutely no customization. It serves as
      the default plot method for all \code{views}.  }


    \item{xyplot}{\code{signature(x = "formula", data = "gateView")}:
      Scatter plots from a \code{gateView} object. This allows for
      customization of the \code{gateView} plots but still provides
      reasonable defaults, e.g., for added gates. }

} 
}



\details{
  
  The implementation of \code{xyplot} in \code{flowViz} is very close to
  the original \code{lattice} version. Concepts like conditioning and
  the use of panels apply directly to the flow cytometry data. The
  single fundamental difference is that conditioning variables are not
  evaluated in the context of the raw data, but rather in the
  \code{phenoData} slot environment (only for the \code{flowSet}
  methods. Thus, we can directly condition on pheotypic variables like
  sample groups, patients or treatments.

  In the formula interface, the primary and secondary variables
  (separated by the tilde) have to be valid parameter names. Please note
  that frequently used variants like \code{FSC-H} and \code{SSC-H} are
  not syntatically correct R symbols, and need to be wrapped in
  \code{`}. E.g., \code{`FSC-H`}. For \code{flowSets}, the use of a
  conditioning variable is optional. We implicitely condition on
  \code{flowFrames} and the default is to arrange panels by sample
  names.
  
  
}

\seealso{
  
  Not all standard lattice arguments will have the intended
  effect, but many should.  For a fuller description of possible
  arguments and their effects, consult documentation on lattice.

}


\author{ F. Hahne, D. Sarkar }

\examples{

data(GvHD)
GvHD <- GvHD[pData(GvHD)$Patient \%in\% 5:6]

## a bivariate scatterplot
## by default ('smooth=TRUE') panel.smoothScatter is used
xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], nbin = 100,
main="A single flowFrame")

## A non-smooth version of the same data
xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], nbin = 100,
main="A single flowFrame", smooth=FALSE)

## Visual artifacts created by the pileup of margin events
xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], nbin = 100,
main="A single flowFrame", margin=FALSE)


## simple bivariate scatter plot (a.k.a. dot plot)
## for the whole flowSet, conditioning on Patient and
## Visit
xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD)


## several examples with time on the X axis
## first for a flowFrame
xyplot(GvHD[[1]])

## and for flowSets
xyplot(`FSC-H` ~ Time | Visit, GvHD, 
       smooth = FALSE, type = "l", 
       subset = (Patient == 5))

xyplot(`FSC-H` ~ Time | Patient+Visit, GvHD, 
       smooth = FALSE, type = "a",
       strip = FALSE, strip.left = TRUE,
       aspect = "xy")


## combine plots for two channels
ssc.time <- 

    xyplot(`SSC-H` ~ Time | factor(Patient):factor(Visit), GvHD, 
           smooth = FALSE, type = "a",
           strip = FALSE,
           strip.left = strip.custom(horizontal = TRUE),
           par.strip.text = list(lines = 3),
           between = list(y = rep(c(0, 0.5), c(6, 1))),
           scales = list(x = list(axs = "i"), y = list(draw = FALSE)),
           layout = c(1, 14))

fsc.time <- 
    xyplot(`FSC-H` ~ Time | factor(Patient):factor(Visit), GvHD, 
           smooth = FALSE, type = "a",
           strip = FALSE,
           strip.left = strip.custom(horizontal = TRUE),
           par.strip.text = list(lines = 3),
           between = list(y = rep(c(0, 0.5), c(6, 1))),
           scales = list(x = list(axs = "i"), y = list(draw = FALSE)),
           layout = c(1, 14))

plot(fsc.time, split = c(1, 1, 2, 1))
plot(ssc.time, split = c(2, 1, 2, 1), newpage = FALSE)


## saving plots as variables allows more manipulation
plot(update(fsc.time[8:14], layout = c(1, 7)),
     split = c(1, 1, 1, 2))

plot(update(ssc.time[8:14], layout = c(1, 7)),
     split = c(1, 2, 1, 2), newpage = FALSE)


## displaying filters
n2gate <- norm2Filter("SSC-H", "FSC-H")

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD,
       filter=n2gate, subset=Patient==5)

xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit,
       data=transform("SSC-H"=asinh,"FSC-H"=asinh) \%on\% GvHD,
       smooth=FALSE, filter=n2gate, subset=Patient==5)


n2gate.results <- filter(GvHD, n2gate)

xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD,
       subset=Patient == "6",
       filter=n2gate.results, smooth=FALSE)

}

\keyword{methods}
\keyword{dplot}
