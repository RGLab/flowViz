\name{flowViz.par.get}
\alias{flowViz.par.get}
\alias{flowViz.par.set}

\title{ Query and set session-wide graphical parameter defaults. }

\description{

  \code{flowViz.par.get} is the equivalent to
  \code{\link[lattice]{trellis.par.get}}. It queries the session wide
  defaults for all \code{lattice} and \code{flowViz} graphical
  parameters.

  \code{flowViz.par.set} is the equivalent to
  \code{\link[lattice:trellis.par.get]{trellis.par.set}}. It sets the
  same set of graphical parameters, either in the \code{flowViz} package
  or directly in \code{lattice}.
}

\usage{

flowViz.par.get(name = NULL)

flowViz.par.set(name, value, ..., theme, warn = TRUE, strict = FALSE, reset = FALSE)

}

\arguments{
  
  \item{name}{The name of a parameter category to set. }
  
  \item{value}{ A named list of values to set for category \code{name}
    or a list of such lists if \code{name} is missing.  }
  
  \item{\dots}{Further arguments that get passed on. }

  \item{theme}{ The theme to set. See
    \code{\link[lattice:trellis.par.get]{trellis.par.set}} for details. }

  \item{warn}{ This gets passed on directly to \code{trellis.par.set}. }

  \item{strict}{ This gets passed on directly to \code{trellis.par.set}. }
  
  \item{reset}{ \code{logical} scalar. When TRUE, drop the entire list of old graphical parameters and reset it with the supplied one. Default is FALSE, which updates the existing parameters. }
}

\details{
  
  Getting and setting graphical parameters in \code{flowViz} follows
  exactly the mechanism of the \code{lattice} package. For all purpose
  and intentions, \code{flowViz.par.get} and \code{flowViz.par.get} can
  be viewed as wrappers around their \code{lattice} counterparts
  \code{\link[lattice]{trellis.par.get}} and
  \code{\link[lattice:trellis.par.get]{trellis.par.set}} and you should
  consult their documentation for further details.

  We introduce four new categories of graphical parameters that are
  relevant for \code{flowViz} plots:
  
  \describe{
  \item{gate}{Controls the appearance of gate boundaries in
    \code{xyplots} (if \code{smooth=TRUE}) or of the points within a gate
    region (\code{smooth=FALSE}).  Available parameters are \code{col},
    \code{cex},\code{pch},\code{alpha},\code{lwd},\code{lty} and
    \code{fill}. }

  \item{gate.density}{Controls the appearance of gate boundaries in
    \code{densityplots}.  Available parameters are \code{col},
    \code{alpha},\code{lwd},\code{lty} and \code{fill}. }

  \item{flow.symbol}{Controls the appearance of 'regular' points in a
    \code{flowViz} plot.  Available parameters are \code{col},
    \code{cex},\code{pch},\code{alpha} and \code{fill}. }
  
  \item{gate.text}{Controls the appearance of the text used for gate
    names.  Available parameters are \code{col},
    \code{cex},\code{font},\code{alpha} and \code{lineheight}. }
  }
}

\value{
  
  \code{flowViz.par.get} returns a list of graphical parameter
  defaults, if \code{name} is not empty, only for this particular
  category. For an empty \code{name} argument, the function returns all
  parameter defaults, including the ones specified in the lattice
  package.

  \code{flowViz.par.set} is called for its side-effects of setting
  default parameters.
  
}

\note{

  Because parameter settings in \code{lattice} are device-dependent,
  \code{flowViz.par.get} will open a (default) device none is open at the
  time of the query.

}

\references{
  Deepayan Sarker, \emph{Lattice, Multivariate Data Visualization with R},
  Springer, New York, 2008
}

\author{ F. Hahne}


\seealso{

  \code{\link[lattice]{trellis.par.get}} and
  \code{\link[lattice:trellis.par.get]{trellis.par.set}}

}

\examples{

## Return all available parameters, including lattice ones
flowViz.par.get()

## Set the font for gate names 
flowViz.par.set("gate.text", list(font=2))

## Query only the gate.text category
flowViz.par.get("gate.text")

## Set a lattice parameter
plot.symbol <- trellis.par.get("plot.symbol")
flowViz.par.set("plot.symbol", list(col="red"))
trellis.par.get("plot.symbol")

## undo all settings
flowViz.par.set(list(plot.symbol=plot.symbol, gate.text=list(font=1)))

data(GvHD)
fs <- GvHD[1:2]

# using default ggplot2like theme
densityplot(~`FSC-H`, fs)
xyplot(`SSC-H`~`FSC-H`, fs, smooth = FALSE)
            
# reset it with default lattice theme            
flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)
densityplot(~`FSC-H`, fs)
xyplot(`SSC-H`~`FSC-H`, fs, smooth = FALSE)

}


