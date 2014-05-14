\name{glpolygon-methods}
\docType{methods}
\alias{glpolygon}
\alias{glpolygon-methods}
\alias{glpolygon,curv1Filter,ANY-method}
\alias{glpolygon,complementFilter,ANY-method}
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
\alias{glpolygon,ellipsoidGate,character-method}
\alias{glpolygon,ellipsoidGate,filterResult-method}
\alias{glpolygon,ellipsoidGate,flowFrame-method}
\alias{glpolygon,subsetFilter,ANY-method}


\title{ Drawing filter regions }

\description{ These methods extend the lattice
  \code{\link[lattice:llines]{lpolygon}} 
  methods for drawing of \code{\link[flowCore:filter-class]{filter}}
  regions. They allow for multiple dispatch, since not all
  \code{\link[flowCore:filter-class]{filter}} types need to be evaluated for
  plotting, but this decision should be made internally.  }

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
      \code{\link[flowCore:filterResult-class]{filterResult}}
      object. This basically 
      extracts the \code{\link[flowCore:filter-class]{filter}} from the
      \code{\link[flowCore:filterResult-class]{filterResult}} and dispatches on
      that. }

    \item{x = "filterResult", data = "flowFrame"}{ For some
      \code{\link[flowCore:filter-class]{filter}} types we need the raw data
      to re-evaluate the filter. }

    \item{x = "curv1Filter", data = "ANY"}{ We either need a
      \code{\link[flowCore:filterResult-class]{filterResult}} or the raw
      data as a 
      \code{\link[flowCore:flowFrame-class]{flowFrame}} for
      \code{\link[flowStats:curv1Filter-class]{curv1Filter}}s. }

    \item{x = "curv1Filter", data = "flowFrame"}{ see above }
    
    \item{x = "curv1Filter", data = "missing"}{ see above }
    
    \item{x = "curv1Filter", data = "multipleFilterResult"}{ see above }
    
    \item{x = "curv2Filter", data = "ANY"}{ We either need a
      \code{\link[flowCore:filterResult-class]{filterResult}} or the raw
      data as a 
      \code{\link[flowCore:flowFrame-class]{flowFrame}} for
      \code{\link[flowStats:curv2Filter-class]{curv2Filter}}.}

    \item{x = "curv2Filter", data = "flowFrame"}{ see above }

    \item{x = "curv2Filter", data = "multipleFilterResult"}{ see above }

    \item{x = "kmeansFilter", data = "ANY"}{ We don't know how to plot
      regions of a \code{\link[flowCore:kmeansFilter]{kmeansFilter}}, hence we
      warn. }

    \item{x = "norm2Filter", data = "ANY"}{ We either need a
      \code{\link[flowCore:filterResult-class]{filterResult}} or the raw
      data as a 
      \code{\link[flowCore:flowFrame-class]{flowFrame}} for
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


\details{ When plotting \code{\link[flowCore:flowFrame-class]{flowFrames}}
  using the any of 
  the lattice-type \code{plot} method provided by \code{flowViz}, the
  plotted parameters are recorded, which makes it possible to correctly
  overlay the outlines of \code{\link[flowCore:filter-class]{filter}}
  assuming that 
  they are defined for the respective parameters. Warnings and error will
  be cast for the cases where the parameters are non-distinct or
  ambigious. These methods are meant to be used within lattice panel
  functions and are probably not of much use outside of those. }

\value{The methods will return the outlines of the gate region as
  polygon vertices.}

\author{ F. Hahne }

\seealso{
  
  \code{\link[flowCore:filter-class]{filter}},
  \code{\link[flowCore:flowFrame-class]{flowFrame}},
  \code{\link[flowViz:lpoints-methods]{glpoints}}
  
}


\keyword{methods}

