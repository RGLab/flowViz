## A nasty hack to get the current parameters that are plotted which need to
## be passed on to panel.xyplot.flowframe
panel.splom.flowframe <- function(x,
                                  frame,
                                  ...)
{
    cv <- current.viewport()$name
    dims <- as.numeric(strsplit(cv, ".", fixed=TRUE)[[1]][2:3])
    cn <- colnames(frame)
    channel.x.name <- cn[dims[1]]
    channel.y.name <- cn[dims[2]]
    ranges <- range(frame)
    xlim <- if(!length(grep("`", channel.x.name, fixed=TRUE))){
        tmp <- ranges[, channel.x.name]
        xd <- diff(tmp)/15
        tmp + c(-1,1)*xd
    }else NULL
    ylim <- if(!length(grep("`", channel.y.name, fixed=TRUE))){
        tmp <- ranges[, channel.y.name]
        yd <- diff(tmp)/15
        tmp + c(-1,1)*yd
    }else NULL
    
    panel.xyplot.flowframe.old(x=x, channel.x.name=channel.x.name,
                           channel.y.name=channel.y.name, frame=frame
				   			,xlim=xlim,ylim=ylim
                           ,...)
}


#' Method implementing Lattice scatter plot matrices for flow data.
#' 
#' 
#' This function create Trellis scatter plots matrices (splom) from flow
#' cytometry data.
#' 
#' 
#' The function draws a scatter plot matrix of the data for each flow parameter
#' in a \code{flowFrame}. For the most, one can think about this as a
#' rectangular arrangement of separate \code{\link[flowViz:xyplot]{xyplots}},
#' and most of that functionality is also available here. To be more precise,
#' the function repeatedly calls \code{panel.xyplot.flowframe} to do the actual
#' plotting. Please see its documentation for details.
#' 
#' @name splom
#' @aliases splom splom,flowFrame,missing-method panel.splom.flowframe
#' @docType methods
#' @param x A formula describing the structure of the plot and the variables to
#' be used in the display.
#' @param data A \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#' that serves as the source of data.
#' @param pscales This arguments is passed unchanged to the corresponding
#' methods in lattice, and is listed here only because it provides a different
#' default.  See documentation for the original methods for details.
#' @param time A character string giving the name of the data column recording
#' time. If not provided, we try to guess from the available parameters.
#' @param exclude.time Logical, specifying whether to exclude the time variable
#' from a scatter plot matrix. Defaults to \code{TRUE}.
#' @param names Logical specifying wether gate names should be added to the
#' plot. Currently, this feature is not supported for splom plots.
#' @param \dots More arguments, usually passed on to the underlying lattice
#' methods.
#' @author F. Hahne, D. Sarkar
#' @seealso
#' 
#' Not all standard lattice arguments will have the intended effect, but many
#' should.  For a fuller description of possible arguments and their effects,
#' consult documentation on lattice.
#' @keywords methods dplot
#' @examples
#' 
#' library(flowCore)
#' data(GvHD)
#' library(flowStats)
#' 
#' tf <- transformList(colnames(GvHD)[3:7], asinh)
#' dat <- tf %on% GvHD[[3]]
#' 
#' 
#' ## scatter plot matrix of individual flowFrames
#' lattice.options(panel.error=NULL)
#' splom(dat)
#' 
#' splom(dat[,1:3], smooth = FALSE)
#' 
#' 
#' ## displaying filters
#' rg <- rectangleGate("FSC-H"=c(200,400), "SSC-H"=c(300,700),
#' "FL1-H"=c(2,4), "FL2-A"=c(4,7))
#' splom(dat, filter=rg)
#' 
#' splom(dat, filter=rectangleGate("FSC-H"=c(400,800)))
#' 
#' splom(dat[,1:4], smooth = FALSE, filter=norm2Filter("FSC-H", "SSC-H", scale=1.5))
#' 
#' 
#' 
#' @export 
setMethod("splom",
          signature=signature(x="flowFrame", data="missing"),
          definition=function(x,
                              data, 
                              pscales,
                              time,
                              exclude.time=TRUE,
                              names=FALSE,
                              ...)
      {
          if(names){
              warning("Filter names are currently not supported for ",
                      "splom plots.", call.=FALSE)
              names=FALSE
          }
          if(missing(pscales))
             pscales <- lapply(range(x), function(x) list(limits=x, at=numeric(0)))
          gp <- list(...)[["par.settings"]]
          expr <- exprs(x)
          column.names <- colnames(expr)
          ## guess the time parameter and exclude if necessary
          expr <- exprs(x)
          if(missing(time))
              time <- flowCore:::findTimeChannel(expr)
          if(exclude.time && length(time))
              column.names <- column.names[column.names != time]
          splom(exprs(x[, column.names]),
                pscales=pscales, 
                panel=panel.splom.flowframe,
                frame=x,
                gp=gp,
                names=names,
                ...)
      })


