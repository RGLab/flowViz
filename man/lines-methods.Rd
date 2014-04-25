\name{glines-methods}
\docType{methods}
\alias{glines}
\alias{glines-methods}
\alias{glines,curv1Filter,ANY-method}
\alias{glines,curv1Filter,flowFrame-method}
\alias{glines,curv1Filter,missing-method}
\alias{glines,curv1Filter,multipleFilterResult-method}
\alias{glines,curv2Filter,ANY-method}
\alias{glines,curv2Filter,flowFrame-method}
\alias{glines,curv2Filter,multipleFilterResult-method}
\alias{glines,filter,missing-method}
\alias{glines,filterResult,flowFrame-method}
\alias{glines,filterResult,ANY-method}
\alias{glines,kmeansFilter,ANY-method}
\alias{glines,norm2Filter,ANY-method}
\alias{glines,norm2Filter,flowFrame-method}
\alias{glines,norm2Filter,logicalFilterResult-method}
\alias{glines,polygonGate,character-method}
\alias{glines,polygonGate,filterResult-method}
\alias{glines,polygonGate,flowFrame-method}
\alias{glines,quadGate,character-method}
\alias{glines,quadGate,filterResult-method}
\alias{glines,quadGate,flowFrame-method}
\alias{glines,rectangleGate,character-method}
\alias{glines,rectangleGate,filterResult-method}
\alias{glines,rectangleGate,flowFrame-method}
\alias{glines,ellipsoidGate,character-method}
\alias{glines,ellipsoidGate,filterResult-method}
\alias{glines,ellipsoidGate,flowFrame-method}

\title{ Drawing filter boundaries }

\description{ These methods extend the basic graphics
  \code{\link{lines}} methods for drawing of
  \code{\link[flowCore:filter-class]{filter}} boundaries. They allow for
  multiple dispatch, since not all \code{\link[flowCore:filter-class]{filter}}
  types need to be evaluated for plotting, but this decision should be
  made internally.  }

\section{Methods}{
    
  \describe{
     
    \item{x = "filter", data = "missing"}{ General method for all
      objects inheriting from
      \code{\link[flowCore:filter-class]{filter}}. This is used as the
      default when no more explicit method is found. It tries to find
      the plotted parameters from the internal \code{flowViz.state}
      environment. This only works if the flow data has been plotted
      using the \code{plot} or \code{xyplot} methods provided by this
      \code{flowViz} package. }

    \item{x = "filterResult", data = "ANY"}{ General method for all
      \code{\link[flowCore:filterResult-class]{filterResult}}
      object. This basically extracts the
      \code{\link[flowCore:filter-class]{filter}} from the
      \code{\link[flowCore:filterResult-class]{filterResult}} and
      dispatches on that. }

    \item{x = "filterResult", data = "flowFrame"}{ For some
      \code{\link[flowCore:filter-class]{filter}} types we need the raw
      data to re-evaluate the filter. }

    \item{x = "curv1Filter", data = "ANY"}{ We either need a
      \code{\link[flowCore:filterResult-class]{filterResult}} or the raw
      data as a \code{\link[flowCore:flowFrame-class]{flowFrame}} for
      \code{\link[flowStats:curv1Filter-class]{curv1Filter}}. }

    \item{x = "curv1Filter", data = "flowFrame"}{ see above }
    
    \item{x = "curv1Filter", data = "missing"}{ see above }
    
    \item{x = "curv1Filter", data = "multipleFilterResult"}{ see above }
    
    \item{x = "curv2Filter", data = "ANY"}{ We either need a
      \code{\link[flowCore:filterResult-class]{filterResult}} or the raw
      data as a \code{\link[flowCore:flowFrame-class]{flowFrame}} for
      \code{\link[flowStats:curv2Filter-class]{curv2Filter}}.}

    \item{x = "curv2Filter", data = "flowFrame"}{ see above }

    \item{x = "curv2Filter", data = "multipleFilterResult"}{ see above }

    \item{x = "kmeansFilter", data = "ANY"}{ We don't know how to plot
      outlines of a
      \code{\link[flowCore]{kmeansFilter}}, hence we
      warn. }

    \item{x = "norm2Filter", data = "ANY"}{ We either need a
      \code{\link[flowCore:filterResult-class]{filterResult}} or the raw
      data as a \code{\link[flowCore:flowFrame-class]{flowFrame}} for
      \code{\link[flowCore:norm2Filter-class]{norm2Filter}}.}

    \item{x = "norm2Filter", data = "flowFrame"}{ see above }

    \item{x = "norm2Filter", data = "logicalFilterResult"}{ see above }

    \item{x = "polygonGate", data = "character"}{ We can plot a
      \code{\link[flowCore]{polygonGate}} directly from the gate
      definition. }

    \item{x = "polygonGate", data = "filterResult"}{ see above }

    \item{x = "polygonGate", data = "flowFrame"}{ see above }

    \item{x = "quadGate", data = "character"}{ We can plot a
      \code{\link[flowCore]{quadGate}} directly from the gate
      definition. }

    \item{x = "quadGate", data = "filterResult"}{ see above }

    \item{x = "quadGate", data = "flowFrame"}{ see above }

    \item{x = "rectangleGate", data = "character"}{ We can plot a
      \code{\link[flowCore]{rectangleGate}} directly from the gate
      definition. }

    \item{x = "rectangleGate", data = "filterResult"}{ see above }

    \item{x = "rectangleGate", data = "flowFrame"}{ see above }

    \item{x = "ellipsoidGate", data = "character"}{ We can plot a
      \code{\link[flowCore]{rectangleGate}} directly from the gate
      definition. }

    \item{x = "ellipsoidGate", data = "filterResult"}{ see above }

    \item{x = "ellipsoidGate", data = "flowFrame"}{ see above }
    
    
  }
}

\details{ When plotting
  \code{\link[flowCore:flowFrame-class]{flowFrames}} using the
  \code{plot} or \code{xyplot} methods provided by \code{flowViz}, the
  plotted parameters are recorded, which makes it possible to correctly
  overlay the outlines of \code{\link[flowCore:filter-class]{filter}}
  assuming that they are defined for the respective parameters. Warnings
  and error will be cast for the cases where the parameters are
  non-distinct or ambigious.

  The flow parameters plotted can be passed on to any of the methods
  through the optional \code{channels} argument, which always gets
  precedence over automatically detected parameters.

  The methods support all plotting parameters that are available for the
  \code{base} \code{lines} functions.
  
}

\author{ F. Hahne }

\seealso{
  
  \code{\link[flowCore:filter-class]{filter}},
  \code{\link[flowCore:flowFrame-class]{flowFrame}},
  \code{\link[flowViz:points-methods]{gpoints}}
  
}

\keyword{methods}
