#
## ==========================================================================
## Basic plot for fcsFrame object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("flowPlot",
          ## basic plot without a gate specified
          signature(x="flowFrame"),
          function(x, child, filter=NULL, plotParameters=c("FSC-H","SSC-H"),
          	logx=FALSE,logy=FALSE,
          	parent,colParent="Grey",colChild="Blue",
          	showFilter=TRUE,gate.fill="transparent",gate.border="black",
          	xlab,ylab,xlim,ylim,...){
          	require(geneplotter)
            data <- x@exprs
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
