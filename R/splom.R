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
    plotLims(xlim, ylim)
    panel.xyplot.flowframe(x=x, channel.x.name=channel.x.name,
                           channel.y.name=channel.y.name, frame=frame,
                           ...)
}


## A scatter plot matrix for flowFrames. We use the panel.xyplot.flowset
## function for the actual hard work
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


