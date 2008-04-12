\name{gpoints-methods}
\docType{methods}
\alias{gpoints}
\alias{gpoints}
\alias{gpoints-methods}
\alias{gpoints,curv1Filter,flowFrame,character-method}
\alias{gpoints,curv2Filter,flowFrame,character-method}
\alias{gpoints,filter,flowFrame,missing-method}
\alias{gpoints,filter,missing,ANY-method}
\alias{gpoints,filterResult,flowFrame,character-method}
\alias{gpoints,kmeansFilter,flowFrame,character-method}
\alias{gpoints,norm2Filter,flowFrame,character-method}
\alias{gpoints,polygonGate,flowFrame,character-method}
\alias{gpoints,quadGate,flowFrame,character-method}
\alias{gpoints,rectangleGate,flowFrame,character-method}


\title{ Adding points within a gate to a plot }

\description{ These methods extend the basic graphics
  \code{\link{points}} methods for drawing of points contained within a
  \code{\linkS4class{filter}}. They allow for multiple
  dispatch, since not all \code{\linkS4class{filter}} types
  need to be evaluated for plotting, but this decision should be made
  internally. In any case, we need the raw data in the form of a
  \code{\linkS4class{flowFrame}}.}



\section{Methods}{
  
  \describe{

    \item{x = "filter", data = "flowFrame", channels = "missing"}{
      General method for all objects inheriting from
      \code{\linkS4class{filter}}. This is used as the default
      when no more explicit method is found. It tries to find the
      plotted parameters from the internal \code{flowViz.state}
      environment. This only works if the flow data has been plotted
      using the \code{plot} methods provided by this \code{flowViz}
      package. }

    \item{x = "filter", data = "missing", channels = "ANY"}{ This gives
      a useful error message when we don't get what we need. }

    \item{x = "filterResult", data = "flowFrame", channels =
      "character"}{ We can get all the information about a
      \code{\linkS4class{filter}} from its
      \code{\linkS4class{filterResult}} without the need to
      re-evaluate.}
    
   

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

    \item{x = "curv1Filter", data = "flowFrame", channels =
      "character"}{ We evaluate the
      \code{\linkS4class{filter}} on the
      \code{\linkS4class{flowFrame}} and plot the subset of
      selected points. By default, every subpopulation (if there are
      any) is colored differently.}

    \item{x = "curv2Filter", data = "flowFrame", channels = "character"}{
      see above }

    \item{x = "kmeansFilter", data = "flowFrame", channels =
      "character"}{ see above }

    \item{x = "norm2Filter", data = "flowFrame", channels =
      "character"}{ see above }

    \item{x = "polygonGate", data = "flowFrame", channels =
      "character"}{ see above }
    
    \item{x = "quadGate", data = "flowFrame", channels = "character"}{
      see above }
    
    \item{x = "rectangleGate", data = "flowFrame", channels =
      "character"}{ see above }

  }
}

\details{ When plotting code{\linkS4class{flowFrame}}s using the
  \code{plot} method provided by \code{flowViz}, the plotted parameters
  are recorded, which makes it possible to correctly overlay the points
  within \code{\linkS4class{filter}}s assuming that they are defined for
  the repective parameters. Warnings and error will be cast for the
  cases where the parameters are non-distinct or ambigious. }

\author{ F. Hahne }

\seealso{
  
  \code{\linkS4class{filter}}, \code{\linkS4class{flowFrame}},
  \code{\link{glines}}, \code{\link{gpolygon}}
  
  
}

\keyword{methods}

