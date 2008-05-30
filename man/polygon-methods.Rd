\name{gpolygon-methods}
\docType{methods}
\alias{gpolygon}
\alias{gpolygon-methods}
\alias{gpolygon,curv1Filter,ANY-method}
\alias{gpolygon,curv1Filter,flowFrame-method}
\alias{gpolygon,curv1Filter,missing-method}
\alias{gpolygon,curv1Filter,multipleFilterResult-method}
\alias{gpolygon,curv2Filter,ANY-method}
\alias{gpolygon,curv2Filter,flowFrame-method}
\alias{gpolygon,curv2Filter,multipleFilterResult-method}
\alias{gpolygon,filter,missing-method}
\alias{gpolygon,filterResult,flowFrame-method}
\alias{gpolygon,filterResult,missing-method}
\alias{gpolygon,kmeansFilter,ANY-method}
\alias{gpolygon,norm2Filter,ANY-method}
\alias{gpolygon,norm2Filter,flowFrame-method}
\alias{gpolygon,norm2Filter,logicalFilterResult-method}
\alias{gpolygon,polygonGate,character-method}
\alias{gpolygon,polygonGate,filterResult-method}
\alias{gpolygon,polygonGate,flowFrame-method}
\alias{gpolygon,quadGate,character-method}
\alias{gpolygon,quadGate,filterResult-method}
\alias{gpolygon,quadGate,flowFrame-method}
\alias{gpolygon,rectangleGate,character-method}
\alias{gpolygon,rectangleGate,filterResult-method}
\alias{gpolygon,rectangleGate,flowFrame-method}

\title{ Drawing filter regions }

\description{ These methods extend the basic graphics
  \code{\link{polygon}} methods for drawing of
  \code{\link[flowCore:filter-class]{filter}} regions. They allow for multiple
  dispatch, since not all \code{\link[flowCore:filter-class]{filter}} types
  need to be evaluated for plotting, but this decision should be made
  internally.  }

\section{Methods}{
  
  \describe{

    \item{x = "filter", data = "missing"}{ General method for all
     objects inheriting from \code{\link[flowCore:filter-class]{filter}}. This
     is used as the default when no more explicit method is found. It
     tries to find the plotted parameters from the internal
     \code{flowViz.state} environment. This only works if the flow data
     has been plotted using the \code{plot} methods provided by this
     \code{flowViz} package. }

    \item{x = "filterResult", data = "missing"}{ General method for all
      \code{\link[flowCore:filterResult-class]{filterResult}} object. This basically
      extracts the \code{\link[flowCore:filter-class]{filter}} from the
      \code{\link[flowCore:filterResult-class]{filterResult}} and dispatches on
      that. }

    \item{x = "filterResult", data = "flowFrame"}{ For some
      \code{\link[flowCore:filter-class]{filter}} types we need the raw data
      to re-evaluate the filter. }

    \item{x = "curv1Filter", data = "ANY"}{ We either need a
      \code{\link[flowCore:filterResult-class]{filterResult}} or the raw data as a
      \code{\link[flowCore:flowFrame-class]{flowFrame}} for
      \code{\link[flowCore:curv1Filter-class]{curv1Filter}}s. }

    \item{x = "curv1Filter", data = "flowFrame"}{ see above }
    
    \item{x = "curv1Filter", data = "missing"}{ see above }
    
    \item{x = "curv1Filter", data = "multipleFilterResult"}{ see above }
    
    \item{x = "curv2Filter", data = "ANY"}{ We either need a
      \code{\link[flowCore:filterResult-class]{filterResult}} or the raw data as a
      \code{\link[flowCore:flowFrame-class]{flowFrame}} for
      \code{\link[flowCore:curv2Filter-class]{curv2Filter}}s.}

    \item{x = "curv2Filter", data = "flowFrame"}{ see above }

    \item{x = "curv2Filter", data = "multipleFilterResult"}{ see above }

    \item{x = "kmeansFilter", data = "ANY"}{ We don't know how to plot
      regions of a \code{\link[flowCore:kmeansFilter-class]{kmeansFilter}}, hence we
      warn. }

    \item{x = "norm2Filter", data = "ANY"}{ We either need a
      \code{\link[flowCore:filterResult-class]{filterResult}} or the raw data as a
      \code{\link[flowCore:flowFrame-class]{flowFrame}} for
      \code{\link[flowCore:norm2Filter-class]{norm2Filter}}s.}

    \item{x = "norm2Filter", data = "flowFrame"}{ see above }

    \item{x = "norm2Filter", data = "logicalFilterResult"}{ see above }

    \item{x = "polygonGate", data = "character"}{ We can plot a
      code{linkS4class{polygonGate}} directly from the gate definition. }

    \item{x = "polygonGate", data = "filterResult"}{ see above }

    \item{x = "polygonGate", data = "flowFrame"}{ see above }

    \item{x = "quadGate", data = "character"}{ We can plot a
      code{linkS4class{quadGate}} directly from the gate
      definition. }

    \item{x = "quadGate", data = "filterResult"}{ see above }

    \item{x = "quadGate", data = "flowFrame"}{ see above }

    \item{x = "rectangleGate", data = "character"}{ We can plot a
      code{linkS4class{rectangleGate}} directly from the gate
      definition. }

    \item{x = "rectangleGate", data = "filterResult"}{ see above }

    \item{x = "rectangleGate", data = "flowFrame"}{ see above }

  }
}


\details{ When plotting code{\link[flowCore:flowFrame-class]{flowFrame}}s using the
  \code{plot} method provided by \code{flowViz}, the plotted parameters
  are recorded, which makes it possible to correctly overlay the
  outlines of \code{\link[flowCore:filter-class]{filter}}s assuming that they are
  defined for the repective parameters. Warnings and error will be cast
  for the cases where the parameters are non-distinct or ambigious. }

\author{ F. Hahne }

\seealso{
  
  \code{\link[flowCore:filter-class]{filter}},
  \code{\link[flowCore:flowFrame-class]{flowFrame}}, \code{\link{glines}},
  \code{\link{gpoints}}
  
}


\keyword{methods}

