

prepanel.xyplot.flowset <- 
    function(x, 
             frames, channel.x, channel.y,
             ...)
{
    if (length(nm <- as.character(x)) > 1)
        stop("must have only one flow frame per panel")
    if (length(nm) == 1)
    {
        xx <- evalInFlowFrame(channel.x, frames[[nm]])
        yy <- evalInFlowFrame(channel.y, frames[[nm]])
        list(xlim = range(xx, finite = TRUE),
             ylim = range(yy, finite = TRUE),
             dx = diff(xx), dy = diff(yy))
    }
    else list()
}

panel.xyplot.flowset <-
    function(x, 
             frames,
             channel.x, channel.y,
             channel.x.name, channel.y.name, 
             filter = NULL,
             filterResults = NULL,

             pch, smooth,
             ...)
{
    x <- as.character(x)
    if (length(x) > 1) stop("must have only one flow frame per panel")
    if (length(x) < 1) return()

    nm <- x
    xx <- evalInFlowFrame(channel.x, frames[[nm]])
    yy <- evalInFlowFrame(channel.y, frames[[nm]])

    ## this.filter.result represents result of applying
    ## filter, which is not necessarily
    ## filterResults[[nm]] (and we have no way of
    ## knowing if it is, so we will always need to
    ## recompute it).  However, it is only required if
    ## filter is specified, and then only for certain
    ## types of filters (specifically, those that have a
    ## 2D geometric representation), and even then only
    ## if the parameters match.  All these decisions are
    ## made inside the call to filterBoundary() later,
    ## but we define the 'this.filter.result' variable
    ## now because we don't want to compute it twice.

    this.filter.result <- NULL

    groups <- 
        if (!is.null(filterResults))
        {
            filterResults[[nm]]@subSet
        }
        else if (!is.null(filter))
        {
            this.filter.result <- filter(frames[[nm]], filter)
            this.filter.result@subSet
        }
        else NULL

    if (smooth) {
        panel.smoothScatter(xx, yy, ...)
    }
    else panel.xyplot(xx, yy, pch = pch,
                      groups = groups,
                      subscripts = seq_along(groups),
                      ...)

    if (!is.null(filter) && (is.list(displayFilter) || displayFilter))
    {
        display.pars <- if (is.list(displayFilter)) displayFilter else list(border = TRUE)
        filter.boundary <-
            filterBoundary(filter = filter,
                           parameters = c(channel.x.name, channel.y.name),
                           frame = frames[[nm]],
                           result = this.filter.result)
        do.call(panel.polygon,
                c(filter.boundary, display.pars))
    }
    
    ##                   if(displayFilter) {
    ##                       if(class(filter)!="rectangleGate") {stop("Only rectangleGate is supported for displayFilter")}			       
    ##                       hLine=c(attr(filter,"min")[gsub("`","",paste(channel.y))],attr(filter,"max")[gsub("`","",paste(channel.y))])
    ##                       vLine=c(attr(filter,"min")[gsub("`","",paste(channel.x))],attr(filter,"max")[gsub("`","",paste(channel.x))])
    ##                       if(!sum(is.na(hLine))) {panel.abline(h=hLine,col="red")}
    ##                       if(!sum(is.na(vLine))) {panel.abline(v=vLine,col="red")}
    ##                   }
    
}



setMethod("xyplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab, ylab,
                   as.table = TRUE,
                   prepanel = prepanel.xyplot.flowset,
                   panel = panel.xyplot.flowset,
                   pch = ".", smooth = TRUE,
                   filter = NULL,
                   filterResults = NULL,
                   displayFilter = TRUE,
                   ...)
      {
          pd <- pData(phenoData(data))
          uniq.name <- createUniqueColumnName(pd)
          ## ugly hack to suppress warnings about coercion introducing
          ## NAs (needs to be `undone' inside prepanel and panel
          ## functions):
          pd[[uniq.name]] <- factor(sampleNames(data)) 
          channel.y <- x[[2]]
          channel.x <- x[[3]]
          if (length(channel.x) == 3)
          {
              channel.x <- channel.x[[2]]
              x[[3]][[2]] <- as.name(uniq.name)
              x[[2]] <- NULL
          }
          else
          {
              x[[3]] <- as.name(uniq.name)
              x[[2]] <- NULL
          }
          channel.x.name <- expr2char(channel.x)
          channel.y.name <- expr2char(channel.y)
          channel.x <- as.expression(channel.x)
          channel.y <- as.expression(channel.y)

          if (missing(xlab)) xlab <- channel.x.name
          if (missing(ylab)) ylab <- channel.y.name

          densityplot(x, data = pd, 

                      prepanel = prepanel,
                      panel = panel,

                      frames = data@frames,
                      channel.x = channel.x,
                      channel.y = channel.y,
                      channel.x.name = channel.x.name,
                      channel.y.name = channel.y.name,
                      
                      filter = filter,
                      filterResults = filterResults,
                      as.table = as.table,

                      xlab = xlab,
                      ylab = ylab,
                      pch = pch, smooth = smooth,

                      ...)
      })






## xyplot for flowFrames, plotting columns against time

setMethod("xyplot",
          signature(x = "flowFrame", data = "missing"),
          function(x, data, xlab = time, ylab = "", time = "Time",
                   layout,
                   type = "l", smooth = FALSE,
                   ...)
      {
          expr <- exprs(x)
          if (!(time %in% colnames(expr)))
              stop("Name of the Time (X) variable must be specified as the 'time' argument")
          time.x <- expr[, time]
          fakedf <- data.frame(channel = setdiff(colnames(expr), time), time = 1)
          prepanel.xyplot.flowframe <- 
              function(x, y, time.x, expr, ...)
              {
                  xx <- time.x
                  yy <- expr[, as.character(y)]
                  list(xlim = range(xx, finite = TRUE),
                       ylim = range(yy, finite = TRUE),
                       dx = diff(xx), dy = diff(yy))
              }
          panel.xyplot.flowframe <- 
              function(x, y, time.x, expr, smooth = FALSE, ...)
              {
                  xx <- time.x
                  yy <- expr[, as.character(y)]
                  if (smooth) panel.smoothScatter(xx, yy, ...)
                  else panel.xyplot(xx, yy, ...)
              }
          if (missing(layout)) layout <- c(1, ncol(expr) - 1)
          xyplot(channel ~ time | channel, data = fakedf,
                 prepanel = prepanel.xyplot.flowframe,
                 panel = panel.xyplot.flowframe,
                 type = type, smooth = smooth,
                 time.x = time.x, expr = expr,
                 xlab = xlab, ylab = ylab,
                 layout = layout,
                 default.scales = list(y = list(relation = "free", rot = 0)),
                 ...)
      })


## xyplot with a formula.  We'll make this very simple; the drawback
## being that the expression matrix will be copied, the upshot being
## that all the fancy xyplot formula stuff will be valid.

setMethod("xyplot",
          signature(x = "formula", data = "flowFrame"),
          function(x, data,
                   smooth = TRUE, 
                   panel = if (smooth) panel.smoothScatter else panel.xyplot,
                   ...)
      {
          xyplot(x, data = as.data.frame(exprs(data)),
                 smooth = smooth, 
                 panel = panel,
                 ...)
      })


