\name{addName-methods}
\docType{methods}
\alias{addName}
\alias{addName-methods}
\alias{addName,curv1Filter,character-method}
\alias{addName,curv1Filter,logical-method}
\alias{addName,curv2Filter,character-method}
\alias{addName,curv2Filter,logical-method}
\alias{addName,ellipsoidGate,character-method}
\alias{addName,ellipsoidGate,logical-method}
\alias{addName,kmeansFilter,character-method}
\alias{addName,kmeansFilter,logical-method}
\alias{addName,polygonGate,character-method}
\alias{addName,polygonGate,logical-method}
\alias{addName,quadGate,character-method}
\alias{addName,quadGate,logical-method}
\alias{addName,quadGate,matrix-method}
\alias{addName,rectangleGate,character-method}
\alias{addName,rectangleGate,logical-method}

\title{Add gate names to a flowViz plot.}

\description{
  
 These methods add gate names to a \code{flowViz} plot, either derived
 from the population identifiers or as provided by the user. These
 methods are ment for internal use and are usually not called directly
 by the user.
 
}
\section{Methods}{
  \describe{

    \item{x = "curv1Filter", name = "character"}{ User-provided names. }
    
    \item{x = "curv1Filter", name = "logical"}{ Get names from the
      \code{\link[flowCore:filter-class]{filter}} or
      \code{\link[flowCore:filterResult-class]{filterResult}} object}
    
    \item{x = "curv2Filter", name = "character"}{ see above }
    
    \item{x = "curv2Filter", name = "logical"}{ see above }
  
    \item{x = "ellipsoidGate", name = "character"}{ see above }
    
    \item{x = "ellipsoidGate", name = "logical"}{ see above }
    
    \item{x = "kmeansFilter", name = "character"}{ see above }
    
    \item{x = "kmeansFilter", name = "logical"}{ see above }
    
    \item{x = "polygonGate", name = "character"}{ see above }
    
    \item{x = "polygonGate", name = "logical"}{ see above }
    
    \item{x = "quadGate", name = "character"}{ see above }
  
    \item{x = "quadGate", name = "logical"}{ see above }
    
    \item{x = "quadGate", name = "matrix"}{ see above }
    
    \item{x = "rectangleGate", name = "character"}{ see above }
    
    \item{x = "rectangleGate", name = "logical"}{ see above }

  }
}

\value{The methods are called for their side effects. No value is
  returned. }

\author{ F. Hahne }

\keyword{methods}

