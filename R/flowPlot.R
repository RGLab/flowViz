#' Standard Plots for Flow Cytometry Data
#' 
#' A method that makes standard plots from a \code{flowFrame}. The user may
#' also provide various \code{filter} or \code{filterResult} arguments to
#' customize the plot.
#' 
#' The plot that is most commonly used in flow cytometry data analysis is
#' usuall called a "dot plot". In common statistical language, we would call
#' this a scatter plot. The basic idea is a 2-dimensional plot that shows the
#' location of every cell in regard to the measurements made on it, for
#' example, forward scatter vs side scatter. Most applications will, in
#' addition to the data, want to show information about one or more filters
#' (gates). Since there can be a very large number of cells in a sample, it is
#' common to show a smoothed version of the data that doesn't involve
#' registering every point on the graph.
#' 
#' @name flowPlot
#' @aliases flowPlot flowPlot,flowFrame-method
#' @docType methods
#' @param x An object of class \code{flowFrame} that contains the data to be
#' plotted.
#' @param child An optional argument of class \code{filterResult} that
#' specifies a subset of the data that are included in the \code{filterResult}
#' @param filter A \code{filter}, \code{filterResult} or
#' \code{filterResultList} object.
#' @param plotParameters A vector of charactors defining the x and y variables
#' in terms of columns in the data.
#' @param logx,logy Logical controlling wheterh the corresponding variables
#' will be log transfromed before passing to the panel function. Default to
#' \code{FALSE}.
#' @param parent An optional argument of class \code{filterResult} that
#' specifies a subset of the data that are inclueed in the \code{filterResult}.
#' @param colParent Specifying the color for \code{parent}. See \code{parent}
#' above.
#' @param colChild Specifies the color for \code{child}. See \code{chile}
#' above.
#' @param showFilter Logical, specifying whether to show the \code{filter}.
#' @param gate.fill Specifies the fill color of the gate. Default to
#' \code{transparent}.
#' @param gate.border Character or specifying the color of the gate border.
#' Default to \code{black}.
#' @param xlab,ylab Labels for data axes.
#' @param xlim,ylim Numeric vectors of length 2 specifying axis limits.
#' @param \dots More arguments, usually passed on to the underlying lattice
#' methods.
#' @author P. Haaland
#' @seealso \code{\link[flowCore:flowCore-package]{flowCore}}
#' @keywords methods
#' @examples
#' 
#' library(flowCore)
#' data(GvHD)
#' flowPlot(GvHD[["s5a01"]])
#' flowPlot(transform("SSC-H"=asinh,"FSC-H"=asinh) %on% GvHD[["s5a01"]])
#' 
#' 
#' @export 
setMethod("flowPlot",
          ## basic plot without a gate specified
          signature(x="flowFrame"),
          function(x, child, filter=NULL, plotParameters=c("FSC-H","SSC-H"),
          	logx=FALSE,logy=FALSE,
          	parent,colParent="Grey",colChild="Blue",
          	showFilter=TRUE,gate.fill="transparent",gate.border="black",
          	xlab,ylab,xlim,ylim,...){
            data <- exprs(x)
            ncells <- nrow(data)
            if(missing(xlab)){
              if(is.character(plotParameters)){
                xlab <-  plotParameters[1]
              }
              else {
                xlab <-  colnames(data)[plotParameters[1]]
              }
            }
            if(missing(ylab)){
              if(is.character(plotParameters)){
                ylab <-  plotParameters[2]
              }
              else {
                ylab <-  colnames(data)[plotParameters[2]]
              }
            }
            ## check for the parent. If it is a filterResult get the parentSet
            if(missing(parent)) {
   				## if the parent is missing assume that it is all cells
				selectParent = rep(TRUE,ncells)
 			}
			else {  	
				if(class(parent) != "logicalFilterResult") {
					stop("parent must be of class filterResult.")
				}
				else {
					if(is.null(parent@subSet)) {
						stop("The filterResult parent must have a subSet slot.")
					}
					else {
						selectParent = parent@subSet
					}  	
				}
			}
			## get the data and log if necessary
			xd =    data[selectParent,plotParameters[1]]
            yd =    data[selectParent,plotParameters[2]]   
            if(logx) {
            	xd[xd <=1] = 1
            	xd = log(xd,b=10)
            }
          	if(logy) {
          		yd[yd <=1] = 1
              	yd = log(yd,b=10)
              }

            ## the default limits are going to be 0,1 or the range of the data
            if(missing(xlim)) {
              	if(max(xd) <= 1.0) {
              		xlim = c(0,1)
              	}
              	else {
              		xlim <-  range(xd)
              	}
              }
              if(missing(ylim)) {
              	if(max(yd) <= 1.0) {
              		ylim = c(0,1)
              	}
              	else {
                	ylim <-  range(yd)
                }
              }

            if(missing(child)) {  
            
             
              smoothScatter(xd,yd,nrpoints=100,
                            xlim=xlim,
                            ylim=ylim,
                            xlab=xlab,
                            ylab=ylab,
                            transformation=function(x) x^(1/1),
                            colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)])),                 
                            ...)       
            }
            else {
              ## check for y and if it is a filterResult, get the subSet
              if(class(child) != "logicalFilterResult") {
                stop("child must be of class filterResult.")
              }
              else {
                if(is.null(child@subSet)) {
                  stop("The filterResult child must have a subSet slot.")
                }
                else {
                  selectChild <-  child@subSet
                }
              }
              ## check for the proper relationship between the parent and child
              ## if that is OK then we need to keep only the parts of child relevant to the parent
#              if(any(selectChild>selectParent)) {
#                stop("The child population must be contained in the parent population.")
#              }
#              else {
               	selectChild = selectChild[selectParent]
               	selectParent = rep(TRUE,length=length(selectChild))
#              }
              selectUnfiltered <- ((selectParent & !selectChild))

              plot(xd[selectChild],yd[selectChild],
                   xlab=xlab,
                   ylab=ylab,
                   xlim=xlim,
                   ylim=ylim,
                   type="n",
                   ...)
              ## filtered population
              points(xd[selectChild],yd[selectChild],col=colChild,...)
              ## unfiltered population
              points(xd[selectUnfiltered],yd[selectUnfiltered],col=colParent,...)
              if(showFilter){
              	childFilter = child@filterDetails[[1]]$filter
              	if(class(childFilter) == "subsetFilter") {
              		childFilter = childFilter@left
              	}
             	if(class(childFilter) == "rectangleGate") {
             		## a simple range gate, note that it may be possible that the graph
             		## isn't based on the parameters of the gate in which case the 
             		## gate won't show
             		if(length(childFilter@parameters)==1){
             			if(childFilter@parameters==plotParameters[1]) {
             				abline(v=childFilter@min,col=gate.border)
             				abline(v=childFilter@max,col=gate.border)
             			}
             			else if(childFilter@parameters==plotParameters[2]) {
             				abline(h=childFilter@min,col=gate.border)
             				abline(h=childFilter@max,col=gate.border)
             			}
             		}
             		## the classic rectangle
             		else if(length(childFilter@parameters)==2) {
             			## check to be sure the parameters actually match
             			if(all(childFilter@parameters %in% plotParameters)) {
             				rectends = c(childFilter@min[plotParameters[1]],childFilter@min[plotParameters[2]],childFilter@max[plotParameters[1]],childFilter@max[plotParameters[2]])
             				## draw the rectangles out to the boundaries if -Inf or Inf appears
             				if(rectends[1] == -Inf) {
             					rectends[1] = xlim[1]
             				}
             				if(rectends[2] == -Inf) {
             					rectends[2] = ylim[1]
             				}
             				if(rectends[3] == Inf) {
             					rectends[3] = xlim[2]
             				}
             				if(rectends[4] == -Inf) {
             					rectends[4] = ylim[2]
             				}
             				rect(rectends[1],rectends[2],rectends[3],rectends[4],col=gate.fill,border=gate.border,...)	
             			}
             		}
				}
# require(ellipse)
#				if(!is.null(attr(parentFilter@subSet,"center"))){
#				elps <- ellipse(x=attr(parentFilter@subSet,"cov"),centre=attr(parentFilter@subSet,"center"))
#				lines(elps)
#			}
#				if(!missing(ellCov)) {
#					elps <- ellipse(x=attr(parentFilter@subSet,"cov"),centre=attr(parentFilter@subSet,"center"))
#					lines(elps)
# 				}
					
 
              }
            }
          }
)
