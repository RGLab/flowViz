\name{densityplot}

\docType{methods}
\alias{densityplot}
\alias{densityplot,formula,flowSet-method}
\alias{densityplot,formula,flowFrame-method}
\alias{densityplot,formula,view-method}
\alias{densityplot,view,missing-method}
\alias{panel.densityplot.flowset}
\alias{prepanel.densityplot.flowset}

\title{ One-dimensional density plots for flow data }



\description{
  
  For \code{\link[flowCore:flowSet-class]{flowSets}} the idea is to
  horizontally stack plots of density estimates for all frames in the
  \code{flowSet} for one or several flow parameters. In the latter case,
  each parameter will be plotted in a separate panel, i.e., we
  implicitely condition on parameters.
  
}

\usage{

## method for 'flowSet' objects
\S4method{densityplot}{formula,flowSet}(
    x,
    data,
    \dots)

## method for 'flowFrame' objects
\S4method{densityplot}{formula,flowFrame}(
    x,
    data, 
    \dots)

prepanel.densityplot.flowset(
    x,
    y,
    frames, 
    overlap=0.3,
    subscripts,
    \dots,
    which.channel)

panel.densityplot.flowset(
    x,
    y,
    darg=list(n=50, na.rm=TRUE),
    frames,
    channel,
    overlap = 0.3,
    channel.name,
    filter=NULL,
    fill=superpose.polygon$col,
    lty=superpose.polygon$lty,
    lwd=superpose.polygon$lwd,
    alpha=superpose.polygon$alpha,
    col=superpose.polygon$border,
    groups=NULL,
    refline=NULL,
    margin=0.005
   ,stats=FALSE
	 ,pos=0.5
	 ,digits=2
	 ,abs=FALSE
	 ,fitGate=TRUE
	 ,checkName = TRUE
    ,gp,
    \dots)



## methods for various workflow objects
\S4method{densityplot}{formula,view}(
    x,
    data,
    \dots)


\S4method{densityplot}{view,missing}(
    x,
    data,
    channels,
    \dots)

}

\arguments{
  
  \item{x}{ A formula describing the structure of the plot and the
    variables to be used in the display. The structure of the formula is
    \code{factor ~ parameter}, where \code{factor} can be any of the
    phenotypic factors in the \code{phenoData} slot or an appropriate
    factor object and \code{parameter} is a flow parameter. Panels for
    multiple parameters are drawn if the formula structure is similar to
    \code{factor ~ parameter1 + parameter2}, and \code{factor} can be
    missing, in which case the sample names are used as y-variable. To
    facilitate programatic access, the formula can be of special
    structure \code{factor ~ .}, in which case the optional
    \code{channel} argument is considered for parameter selection. For
    the workflow methods, \code{x} can also be one of the several
    workflow objects. }
  
  \item{data}{ A flow data object that serves as a source of data,
    either a \code{\link[flowCore:flowFrame-class]{flowFrame}} or
    \code{\link[flowCore:flowSet-class]{flowSet}}, or one of the several
    workflow objects. }
  
  	\item{overlap}{The amount of overlap between stacked density plots. This argument is ignored for the \code{flowFrame} method.}    
	    
	\item{filter}{A \code{\link[flowCore:filter-class]{filter}},
	\code{\link[flowCore:filterResult-class]{filterResult}} or
	\code{\link[flowCore:filterResult-class]{filterResultList}} object
	or a list of such objects of the same length as the
	\code{flowSet}. If applicable, the gate region will be superiposed
	on the density curves using color shading. The software will figure
	out whether the \code{filter} needs to be evaluated in order to be
	plotted (in which case providing a \code{filterResult} can speed
	things up considerably). }
	
	\item{groups}{Use identical colors for grouping. The value of the
			    	argument is expected to be a phenotypic variable in the
			    	\code{flowSet}, or a factor.}
	\item{channels}{A character vector of parameters that are supposed to be plotted when the formula in \code{x} is of structure \code{factor ~ .}.}

  \item{subscripts, which.channel, channel.name, y}{ Internal indices
    necessary to map panels to parameters. }

  \item{frames}{ An environment containing frame-specific data. }

  \item{channel}{ The name of the currently plotted flow parameter. }
	\item{darg}{These arguments are passed unchanged to
									     the corresponding methods in lattice, and are listed here only
									    because they provide different defaults.  See documentation for the
									    original methods for details. \code{darg} gets passed on to
									    \code{\link[stats]{density}}.} 

  \item{col, fill, lty, lwd, alpha}{ Graphical parameters. These
    mostly exist for conveniance and much more control is available
    throught the \code{lattice}-like \code{par.setting} and
    \code{flowViz.par.set} customization. The relevant parameter
    category for density plots is \code{gate.density} with available
    parameters \code{col}, \code{fill}, \code{lwd}, \code{alpha} and
    \code{lty}. See
    \code{\link[flowViz:flowViz.par.get]{flowViz.par.set}} for details.} 
  
  
  \item{refline}{ Logical. Add one ore more vertical reference lines to
    the plot.  This argument is directly passed to
    \code{\link[lattice:panel.functions]{panel.abline}}. }

  \item{margin}{ Numeric in \code{[0,1]}. Indicate margin events by
    horizontal bars. The value of \code{margin} is interpreted as the
    proportion of events on the margin over which the bars are
    added. E.g., a value of \code{0,5} means to indicate margin events
    if there are more than \code{0.5} times the total number of
    events. \code{1} means to ignore margin events completetly. For
    \code{0} bars are added even if there is only a single margin event.}

   \item{stats,pos,digits,abs}{Arguments to control statistics that is associated with \code{\link[flowCore:filter-class]{filter}} to be plotted.
						see \link[flowViz]{xyplot} for details.}
	\item{fitGate}{A \code{logical} scalar indicating whether to display the gate as fitted 1d density gate region 
					or simply display the gate boundaries using vertical lines. The latter would be helpful to display the gate when the gated density region is too small to see. }
	
	\item{checkName}{ A \code{logical} scalar indicating whether to validity check the channel name. Default is TRUE, which consider '(' as invalid character in channel names}
	
  \item{gp}{ A list of graphical parameters that are passed down to
    the low level panel functions. This is for internal use only. The
    public user interface to set graphical parameters is either
    \code{par.settings} for customization of a single call or
    \code{flowViz.par.set} for customization of session-wide defaults. }
  
  \item{\dots}{ More arguments, usually passed on to the underlying
    lattice methods. 
    \itemize{
    	  
       
		 \item xlab: Label for data x axis, with suitable defaults taken from the formula 

 		 
  
  		\item prepanel: The prepanel function.  See \code{\link[lattice]{xyplot}} 
  
  		\item panel: the panel function.  See \code{\link[lattice]{xyplot}} 
		
    	
    	\item axis: axis function passed to lattice, default is \code{axis.grid}
		
		\item ... : other arguments passed to panel.densityplot.flowset
		
 
    	}
    }
}


\section{Methods}{
  \describe{

    \item{densityplot}{\code{signature(x = "formula", data =
	"flowSet")}: Creates density plots for one or several channels,
      with samples stacked according to a \code{phenoData} variable.
      Colors are used to indicate common values of this covariate
      across panels. Filters can be added as the optional
      \code{filter} arguments. See \code{\link{xyplot}} for details. }
    

    \item{densityplot}{\code{signature(x = "formula", data = "view")}: A
      method to create density plots for workspace
      \code{\link[flowCore:view-class]{view}} objects. This still allows
      for some level of customization, but most defaults will be set depending
      on the input object. }

    \item{densityplot}{\code{signature(x = "view", data = "missing")}:
      The default method for view objects. All defaults will be set. }
    
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
GvHD <- GvHD[pData(GvHD)$Patient \%in\% 6:7]

densityplot(~ `FSC-H`, GvHD)

densityplot(~ `FSC-H` + `SSC-H`, GvHD)

densityplot(~ ., GvHD[1:3])

## include a filter
densityplot(~ `FSC-H`, GvHD, filter=curv1Filter("FSC-H"))

#display the gate by its boundaries with statistics 
densityplot(~ `FSC-H`, GvHD[1:2], filter=curv1Filter("FSC-H"),fitGate=FALSE,stats=TRUE)

## plot a single flowFrame
densityplot(~ `SSC-H`, GvHD[[1]], margin=FALSE)

}

\keyword{methods}
\keyword{dplot}
