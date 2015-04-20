##############################################################################
##                            flowSet methods                               ##
##############################################################################


prepanel.histogram.flowset.ex <- 
    function(x, frames, channel.x.name, xlim, type = "density" , margin = TRUE, breaks = 100, ...)
{

  if (length(nm <- as.character(x)) > 1)
    stop("must have only one flow frame per panel")
  
  if (length(nm) == 1){
    ranges <- range(frames[[nm]])
    if(missing(xlim)){
      #use default calculation
      xlim <- if(!length(grep("`", channel.x.name, fixed=TRUE))){
            tmp <- ranges[, channel.x.name]
            xd <- diff(tmp)/15
            tmp + c(-1,1)*xd
          }else NULL
    }
#    browser()
    if(type == "count"){
     #run hist to estimate the max count
     data <- as.vector(exprs(frames[[nm, channel.x.name]]))
     
     if(margin){
       r <- ranges[, channel.x.name]
       pl <- data<=r[1]
       pr <- data>=r[2]
       data <- data[!(pl | pr)]  
     }
     h <- hist(data, plot = FALSE, breaks = breaks)
     count.max <- max(h$counts)
     ylim = c(0, count.max)
    }else
        ylim = c(0,1)
    return(list(xlim=xlim, ylim = ylim))
  }else
    return(list())
}

panel.histogram.flowset.ex <- function(x,
    frames,
    filter = NULL,
    channel.x,
    channel.y
    ,stats = FALSE
    , overlay= NULL#a list of flowSet
    ,...)
{
  
  nm <- as.character(x)
  if (length(nm) < 1) return()
  
  filter <- .processFilter(filter, nm)
  
  if(!is.list(stats)){
    stats <- lapply(seq_along(nm), function(x) stats)
    names(stats) <- nm
    
  }
  
  overlay <- .process_overlay_flowSet(overlay, nm)
  
  panel.histogram.flowFrame(frame=frames[[nm]], filter=filter[[nm]], stats = stats[[nm]], overlay = overlay, ...)
}
panel.histogram.flowFrame <-
    function(frame,
        filter=NULL,
        margin=0.005,
        channel.x.name,
        channel.y.name,
        fill= "#619CFF",
        lty= 1,
        lwd= 1,
        alpha= 1,
        col= "transparent",
        gp
        ,stats = FALSE
        ,pos=0.5
        ,digits=2
        ,abs=FALSE
        ,checkName = TRUE
        , overlay = NULL
        , overlay.symbol = NULL
        , breaks = 100
        ,...
        )
{
          
  
          limits <- current.panel.limits()
          ylim <- limits[["ylim"]]
          validName <- !length(grep("\\(", channel.x.name))
          if(checkName)
            validName <- !(length(grep("\\(", channel.x.name)) ||
                  length(grep("\\(", channel.y.name)))
          else
            validName <- TRUE
          
          if(!validName)
            warning("Gate will not be plotted because channel names contain '(' character! Try to set checkName to FALSE to skip this check.")
          border <- col
          parm <- channel.x.name
          ptList<-plotType("gdensity", parm)
         
          xx <- exprs(frame)[,channel.x.name]
          r <- unlist(range(frame, channel.x.name))
          
          ## we ignore data that has piled up on the margins
          
          if(margin)#when margin is logical FALSE we skip marginal events filtering
          {
            
            pl <- xx<=r[1]
            pr <- xx>=r[2]
            xxt <- xx[!(pl | pr)]
            
            mbar(xx, list(pl, pr), r, 1, col, margin)
          }else
            xxt <- xx
          if(length(xxt)){

            
            panel.histogram(x=xxt, col=fill, border=NA, alpha=alpha, breaks = breaks, ...)
            
           
            if(!is.null(overlay))
            {
              overlayNames <- names(overlay)
              for(overLayName in overlayNames){
#                browser()
                thisOverlay <- overlay[[overLayName]]
                if(!is.null(thisOverlay)){
                  thisDat <- exprs(thisOverlay)
                  overlay.x <- thisDat[, channel.x.name]
                  
                  
                  this.overlay.symbol <- gp$overlay.symbol#default setting
                  user.overlay.symbol <- overlay.symbol[[overLayName]]#customized settings
                  #update the default with customized settings
                  if(!is.null(user.overlay.symbol))
                    this.overlay.symbol <- lattice:::updateList(this.overlay.symbol, overlay.symbol[[overLayName]])
                  
                  panel.histogram(overlay.x 
                      ,  border=NA
                      , col = this.overlay.symbol[["fill"]]
                      , alpha = this.overlay.symbol[["alpha"]]
                      , breaks = breaks
                      , ...
                  )  
                }
                
              }
              
            }                        
            ## add the filterResult if possible, we get them from the output of
            ## glpolygon (with plot=FALSE)
            if(!is.null(filter) && validName){    
              if(is(filter,"filters"))
              {
                fitGate <- FALSE #it is hard to plot multiple fitted gates
                
                mapply(filter,stats,FUN=function(curFilter,curStats){
                      
#                      browser()
                      names <- .getStats(curFilter,curStats, frame, digits, ...)
                      
                      
                      #this plot routine is only for 2-d scatter plot
                      #thus set plot as FALSE,just use it to get bounds
                      #add plot stats/names
                      bounds <- glpolygon(curFilter, frame,
                          channels=parm
                          ,ptList=ptList
                          ,verbose=FALSE
                          , plot=FALSE 
                          ,names=names
                          ,xlim=r
                          ,ylim=c(0,1)
                          ,strict=FALSE)
                      oo <- options(warn=-1)
                      on.exit(options(oo))
                      if(!is.na(bounds)){
                        ## iterate over gate regions
                        for(j in seq_along(bounds)){
                          tb <- bounds[[j]]
                          
                          
                            gpg<-gp$gate
                            panel.lines(x=c(tb[1],tb[1]),y=ylim
                                ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                            panel.lines(x=c(tb[2],tb[2]),y=ylim
                                ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                          }
                          
                      }
                      
                      options(oo)
                    })
              }else
              {
                  names <- .getStats(filter,stats, frame, digits, ...)
                  
        
                  #this plot routine is only for 2-d scatter plot
                  #thus set plot as FALSE,just use it to get bounds
                  #add plot stats/names
                  bounds <- glpolygon(filter, frame,
                      channels=parm
                      ,ptList=ptList
                      ,verbose=FALSE
                      , plot=FALSE 
                      ,names=names
                      ,xlim=r
                      ,ylim=c(0,1)
                      ,strict=FALSE)
                  oo <- options(warn=-1)
                  on.exit(options(oo))
                  if(!is.na(bounds)){
                    ## iterate over gate regions
                    for(j in seq_along(bounds)){
                      tb <- bounds[[j]]
                      
                      gpg<-gp$gate
                      
                      panel.lines(x=c(tb[1],tb[1]),y= ylim
                          ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                      panel.lines(x=c(tb[2],tb[2]),y= ylim
                          ,col=gpg$col,alpha=gpg$alpha, lwd=gpg$lwd,lty=gpg$lty)
                      
                      
                    }
                  }
                  
                  
                  
                  options(oo)
              }
            }
            
    
          }else{
            panel.lines(rl, rep(0,2), col="black")
          }
      
}

setMethod("histogram",
          signature(x = "formula", data = "flowSet"),
          function(x, data, ...){
            
            #construct lattice object
            thisTrellisObj <- .histogram.adapor(x, data, ...) 
              
            #pass frames slot
            thisData <- thisTrellisObj[["panel.args.common"]][["frames"]]
            thisTrellisObj[["panel.args.common"]][["frames"]] <- thisData@frames
            thisTrellisObj
    })

#' dispatch to different trellis object contructing function based on stack argument
#' @param stack \code{logical} indicating whether to stack `name` on y axis
.histogram.adapor <- function(x, data, stack = FALSE, ...){
  
      if(stack){
        stop("stacked histogram is not supported yet!")
#        thisObj <- .histogram.flowSet(x, data, ...) 
      }else
      {
    
        #add dummy y term in order to use .xyplot.flowSet to construct lattice object
        if(length(x[[2]]) == 1)
          xTerm <- x[[2]]
        else
          xTerm <- x[[2]][[2]]
        thisFormula <- eval(substitute(thisX~y, list(thisX = xTerm)))
        thisFormula[[3]] <- x[[2]]
#        browser()
        thisObj <- .xyplot.flowSet(thisFormula, data
            , panel = panel.histogram.flowset.ex
            , prepanel = prepanel.histogram.flowset.ex
            , ylab = ""
            , plotType = "densityplot"
            ,...)
        #append channel.name to work ncdfFlow::densityplot
        thisObj$panel.args.common$channel.name <- thisObj$panel.args.common$channel.x.name
      }
      thisObj
}


## ==========================================================================
## For the flowFrame method, we simply coerce to a flowSet and remove the
## axis annotation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("histogram",
          signature(x="formula", data="flowFrame"),
          function(x, data, overlay = NULL, ...){
            sn <- identifier(data)
            data <- as(data, "flowSet")
            sampleNames(data) <- sn
            
            overlay <- .process_flowFrame_overlay(overlay, sn)
                        
            histogram(x, data, overlap = 0, scales = list(y=list(draw=FALSE)), overlay = overlay, ...)              
          })

