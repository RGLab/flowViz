\name{glpolygon-methods}
\docType{methods}
\alias{glpolygon}
\alias{glpolygon-methods}
\alias{glpolygon,curv1Filter,ANY-method}
\alias{glpolygon,curv1Filter,flowFrame-method}
\alias{glpolygon,curv1Filter,missing-method}
\alias{glpolygon,curv1Filter,multipleFilterResult-method}
\alias{glpolygon,curv2Filter,ANY-method}
\alias{glpolygon,curv2Filter,flowFrame-method}
\alias{glpolygon,curv2Filter,multipleFilterResult-method}
\alias{glpolygon,filter,missing-method}
\alias{glpolygon,filterResult,flowFrame-method}
\alias{glpolygon,filterResult,missing-method}
\alias{glpolygon,filterResult,ANY-method}
\alias{glpolygon,kmeansFilter,ANY-method}
\alias{glpolygon,norm2Filter,ANY-method}
\alias{glpolygon,norm2Filter,flowFrame-method}
\alias{glpolygon,norm2Filter,logicalFilterResult-method}
\alias{glpolygon,polygonGate,character-method}
\alias{glpolygon,polygonGate,filterResult-method}
\alias{glpolygon,polygonGate,flowFrame-method}
\alias{glpolygon,quadGate,character-method}
\alias{glpolygon,quadGate,filterResult-method}
\alias{glpolygon,quadGate,flowFrame-method}
\alias{glpolygon,rectangleGate,character-method}
\alias{glpolygon,rectangleGate,filterResult-method}
\alias{glpolygon,rectangleGate,flowFrame-method}


\title{ Drawing filter regions }

\description{ These methods extend the lattice \code{\link{lpolygon}}
  methods for drawing of \code{\linkS4class{filter}}
  regions. They allow for multiple dispatch, since not all
  \code{\linkS4class{filter}} types need to be evaluated for
  plotting, but this decision should be made internally.  }

\section{Methods}{
  
  \describe{

    \item{x = "filter", data = "missing"}{ General method for all
     objects inheriting from \code{\linkS4class{filter}}. This
     is used as the default when no more explicit method is found. It
     tries to find the plotted parameters from the internal
     \code{flowViz.state} environment. This only works if the flow data
     has been plotted using the \code{plot} methods provided by this
     \code{flowViz} package. }

    \item{x = "filterResult", data = "missing"}{ General method for all
      \code{\linkS4class{filterResult}} object. This basically
      extracts the \code{\linkS4class{filter}} from the
      \code{\linkS4class{filterResult}} and dispatches on
      that. }

    \item{x = "filterResult", data = "flowFrame"}{ For some
      \code{\linkS4class{filter}} types we need the raw data
      to re-evaluate the filter. }

    \item{x = "curv1Filter", data = "ANY"}{ We either need a
      \code{\linkS4class{filterResult}} or the raw data as a
      \code{\linkS4class{flowFrame}} for
      \code{\linkS4class{curv1Filter}}s. }

    \item{x = "curv1Filter", data = "flowFrame"}{ see above }
    
    \item{x = "curv1Filter", data = "missing"}{ see above }
    
    \item{x = "curv1Filter", data = "multipleFilterResult"}{ see above }
    
    \item{x = "curv2Filter", data = "ANY"}{ We either need a
      \code{\linkS4class{filterResult}} or the raw data as a
      \code{\linkS4class{flowFrame}} for
      \code{\linkS4class{curv2Filter}}s.}

    \item{x = "curv2Filter", data = "flowFrame"}{ see above }

    \item{x = "curv2Filter", data = "multipleFilterResult"}{ see above }

    \item{x = "kmeansFilter", data = "ANY"}{ We don't know how to plot
      regions of a \code{\linkS4class{kmeansFilter}}, hence we
      warn. }

    \item{x = "norm2Filter", data = "ANY"}{ We either need a
      \code{\linkS4class{filterResult}} or the raw data as a
      \code{\linkS4class{flowFrame}} for
      \code{\linkS4class{norm2Filter}}s.}

    \item{x = "norm2Filter", data = "flowFrame"}{ see above }

    \item{x = "norm2Filter", data = "logicalFilterResult"}{ see above }

    \item{x = "polygonGate", data = "character"}{ We can plot a
      code{linkS4class{polygonGate}} directly from the gate
      definition. }

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


\details{ When plotting code{\linkS4class{flowFrame}}s using the any of
  the lattice-type \code{plot} method provided by \code{flowViz}, the
  plotted parameters are recorded, which makes it possible to correctly
  overlay the outlines of \code{\linkS4class{filter}}s assuming that
  they are defined for the repective parameters. Warnings and error will
  be cast for the cases where the parameters are non-distinct or
  ambigious. These methods are meant to be used within lattice panel
  functions and are probably not of much use outside of those. }

\author{ F. Hahne }

\seealso{
  
  \code{\linkS4class{filter}},
  \code{\linkS4class{flowFrame}}, \code{\link{glpoints}}
  
}


\keyword{methods}

