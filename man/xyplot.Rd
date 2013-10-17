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
    time, xlim, ylim,
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
	filter = NULL,    
    overlay=NULL,
    stats = FALSE,
    strip.text = NULL,
    \dots)


## prepanel function for generic xyplots of flowFrames
prepanel.xyplot.flowframe(
    frame,
    channel.x.name,
    channel.y.name,
    x, y, xlim, ylim,
    \dots)
   

## panel function for generic xyplots of flowFrames
panel.xyplot.flowframe(
    frame,
    filter=NULL,
    smooth=TRUE,
    margin=TRUE,
    outline=FALSE,
    channel.x.name,
    channel.y.name,
    pch=gp$flow.symbol$pch,
    alpha=gp$flow.symbol$alpha,
    cex=gp$flow.symbol$cex,
    col=gp$flow.symbol$col,
    gp
    ,xbins=0
    ,binTrans=sqrt
    ,stats=FALSE
    ,pos=0.5
    ,digits=2
    ,abs=FALSE
    ,overlay.x=NULL
	,overlay.y=NULL
	,checkName = TRUE
    ,\dots)



## method for 'flowSet' objects
\S4method{xyplot}{formula,flowSet}(
    x,
    data,
    \dots)


## prepanel function for generic xyplots of flowSets
prepanel.xyplot.flowset(
    x,
    frames,
    channel.x.name,
    channel.y.name,
    xlim, ylim,
    \dots)


## panel function for generic xyplots of flowSets
panel.xyplot.flowset(
    x,
    frames,
    filter=NULL,
    channel.x,
    channel.y,
    overlay=NULL,
    stats = FALSE,
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
    various \code{\link[flowCore:view-class]{view}} or
    \code{\link[flowCore:actionItem-class]{actionItem}} objects.}

  
  \item{time}{ A character string giving the name of the data column
    recording time. If not provided, we try to guess from the available
    parameters. }

  
  \item{xlab, ylab}{ Labels for data axes, with suitable defaults taken
    from the formula. }
	
   \item{xlim, ylim}{ limits for data axes. If not given, they are taken from the ranges stored in flowFrame }
  
  \item{layout}{ These arguments are passed unchanged
    to the corresponding methods in lattice, and are listed here only
    because they provide different defaults.  See documentation for the
    original methods for details. }

  
  \item{type}{ type of rendering; see
    \code{\link[lattice:panel.xyplot]{panel.xyplot}} for details. For
    the basic \code{flowFrame} method without a detailed formula, the
    addtional type \code{discrete} is available, which plots a smoothed
    average of the flow cytometry values against time. }

  
  \item{nrpoints}{The number of points plotted on the smoothed plot in
    sparse regions. This is only listed here because we use a different
    default. See
    \code{\link[lattice:panel.smoothScatter]{panel.smoothScatter}}
    for details. }

  
  \item{binSize}{The size of a bin (i.e., the number of events within a
    bin) used for the smoothed average timeline plots. }
	
   \item{overlay}{The extra cell events plotted on top of the current cell population.
   				It is a \code{flowSet} for \code{panel.xyplot.flowset} function and a \code{flowFrame} for \code{xyplot(c("formula","flowFrame"))} method.
   	}

  \item{checkName}{ A \code{logical} scalar indicating whether to validity check the channel name. Default is TRUE, which consider '(' as invalid character in channel names}


  \item{channel.x.name, channel.y.name}{ Character strings giving
    corresponding names used to match filter parameters if applicable. }

  
  \item{smooth}{ Logical. If \code{TRUE}, \code{panel.smoothScatter} is
    used to display a partially smoothed version of the data.
    Otherwise, events are plotted individually, as in a standard scatter
    plot. If \code{FALSE}, a graphical parameter \code{colramp} can be
    used to obtain a coloring of points that is indicative of their
    local density. } 
  
  \item{filter}{ A \code{\link[flowCore:filter-class]{filter}},
    \code{\link[flowCore:filterResult-class]{filterResult}} or
    \code{\link[flowCore:filterResult-class]{filterResultList}} object
    or a list of such objects of the same length as the
    \code{flowSet}.
    Also a \code{\link[flowCore:filters-class]{filters}}
    or A \code{\link[flowCore:filtersList-class]{filtersList}}
    can be passed to xyplot in order to plot multiple filters/gates(with the same x,y parameters)
    on one panel to represent multiple sub-populations. 
    The appropriate spherical 2D representation of this
    filter will be superimposed on the plot if \code{smooth=TRUE}, or
    the result of the filtering operation will be indicated by grouping
    if \code{smooth=FALSE}. The software will figure out whether the
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

    \item{xbins}{The argument passed to \code{\link[hexbin:hexbin]{hexbin}} 
		    	 ,which is the number of bins partitioning the range of xbnds.
		    	 It is set as 0 by default,which plots all the events without binning.
		    	 When it is larger than 0,hexbin plot engine is used for the faster plotting.
		    	 Note that it is only valid when \code{smooth} is set as FALSE .
		    	 }
	\item{binTrans}{The argument passed to \code{\link[hexbin:grid.hexagons]{grid.hexagons}}
					,which is a transformation function (or NULL) for the count. 
					It is \code{\link[base:sqrt]{sqrt}} by default.
		    		}    			 		
	
 	\item{stats,pos,digits,abs}{Arguments to control statistics that is associated with \code{\link[flowCore:filter-class]{filter}} to be plotted
 	 					Currently only population proportion/percentage is supported.
 	 					\code{stats} is a \code{logical} scalar indicating whether to display statistics. Default is FALSE.
						\code{pos} is the \code{numeric} scalar (range within c(0,1)) or vector(length of 2,first is for x-axis,second for y-axis)
								  to control the position of the statistics label. It is set as 0.5,which is the center.
						\code{digits} is an \code{integer} indicating the number of significant digits to be used when displaying the percentage of population statistics,Default is 2. see more details from \code{\link[base]{format}}
						\code{abs} is a \code{logical} scalar indicating whether the \code{pos} is relative to the gate boundary or the entire xy-axis(absolute position).
									By default it is set as FALSE,which indicates the position is relative to gate. 
  						}
  	\item{overlay.x,overlay.y}{ The x,y coordinates of overlay cell points.}  	
  	
  \item{strip.text}{A \code{character} that customizes the text in strip. Default is NULL, which does not display the strip box at all. It is only valid when plotting a \code{flowFrame}}
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
      in the \code{\link[flowCore:flowCore-package]{flowCore}} package rather
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
  not syntactically correct R symbols, and need to be wrapped in
  \code{` `}. E.g., \code{`FSC-H`}. For \code{flowSets}, the use of a
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

## A non-smooth version of the same data with customerized color scheme
require(IDPmisc)
colramp <- colorRampPalette(IDPcolorRamp(21))
xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], nbin = 100,
       main="A single flowFrame", smooth=FALSE,
       colramp=colramp, pch=20, cex=0.1)

## A hexbin version of non-smooth scatter plot  
xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], xbin = 128
       ,main="A single flowFrame", smooth=FALSE)


## Visual artifacts created by the pileup of margin events
xyplot(`FSC-H` ~ `SSC-H`, GvHD[["s5a05"]], nbin = 100,
       main="A single flowFrame", margin=FALSE)


## simple bivariate scatter plot (a.k.a. dot plot)
## for the whole flowSet, conditioning on Patient and
## Visit
xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD)

## Same bivariate scatter plot with replacing default color
require(IDPmisc)
cols <- colorRampPalette(IDPcolorRamp(21))
xyplot(`SSC-H` ~ `FSC-H` | Patient:Visit, data = GvHD, colramp=cols)

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
       smooth=FALSE, filter=n2gate, subset=Patient == 5)


## displaying filters with stats
n2gate.results <- filter(GvHD, n2gate)

xyplot(`SSC-H` ~ `FSC-H` | Visit, data=GvHD,
       subset=Patient == "6",
       filter=n2gate.results, smooth=FALSE
       ,stats=TRUE
       ,abs=TRUE
       ,digits=3 
       )

       
## displaying multiple filters in one panel with stats
recGate1<-rectangleGate("FL3-H"=c(2.3,4.1),"FL2-H"=c(6.8,9))
recGate2<-rectangleGate("FL3-H"=c(1,3),"FL2-H"=c(4,6))
filters1<-filters(list(recGate1,recGate2))
trans<-transform("FL2-H"=asinh,"FL3-H"=asinh)
trans_data<-transform(GvHD[[1]],trans)
xyplot(`FL2-H` ~ `FL3-H`
	   ,data=trans_data 
       ,filter=filters1
       ,stats=TRUE
       ,margin=FALSE
       )

#display recGate2 as a overlay 
overlay<-Subset(trans_data,recGate2)
xyplot(`FL2-H` ~ `FL3-H`
	   ,data=trans_data 
       ,filter=recGate1
       ,stats=TRUE
       ,margin=FALSE
       ,overlay=overlay
       )



}

\keyword{methods}
\keyword{dplot}
