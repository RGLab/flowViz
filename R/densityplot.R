##############################################################################
##                            flowSet methods                               ##
##############################################################################


## Dedicated prepanel function to set up dimensions
prepanel.densityplot.flowset <- 
    function(x, y, darg=list(n=50), frames, channel, overlap=0.3, ...)
{
    xl <- numeric(0)
    for(nm in as.character(x))
    {
        xx <- evalInFlowFrame(channel, frames[[nm]])
        xl <- c(xl, range(xx))
    }
    list(xlim = range(xl, finite=TRUE))
}


## Dedicated panel function to do the plotting and add gate boundaries
panel.densityplot.flowset <-
    function(x, y, darg=list(n=50), ref=FALSE, frames, channel,
             overlap = 0.3, filter=NULL,
             col = superpose.polygon$col,
             lty = superpose.polygon$lty,
             lwd = superpose.polygon$lwd,
             alpha = superpose.polygon$alpha,
             border = superpose.polygon$border,
             gpar=NULL, ...)
{
    superpose.line <- trellis.par.get("superpose.line")
    superpose.polygon <- trellis.par.get("superpose.polygon")
    reference.line <- trellis.par.get("reference.line")
    col.signif <- "red"
    
    ycode <- as.numeric(y)
    if (any(duplicated(ycode)))
        warning("Some combinations seem to have multiple samples.  \n  ",
                "Only one will be used.")
    nnm <- as.character(x)
    if(is(filter, "filter")){
        fl <- vector(length(nnm), mode="list")
        names(fl) <- nnm
        for(n in nnm)
            fl[[n]] <- filter
        filter <- fl
    }
    if(!is.null(filter) && !all(sapply(filter, is, "filter")))
        stop("'filter' must inherit from class 'filter' or a be named list of ",
             "such objects with names matching to frame names.", call.=FALSE)
    
    ny <- nlevels(y)
    col <- rep(col, length = ny)
    lty <- rep(lty, length = ny)
    lwd <- rep(lwd, length = ny)
    alpha <- rep(alpha, length = ny)
    border <- rep(border, length = ny)
    x <- as.character(x)
    height <- (1 + overlap)
    parm <- gsub("`", "", as.character(channel))
    plotType("gdensity", parm)
    if(is.null(gpar)){
        gpar <- flowViz.par()
        gpar$col <- "red"
    } else if(! "col" %in% names(gpar)){
        gpar$col <- "red"
    }
    class(gpar) <- "gpar"
    for (i in rev(seq_len(ny))){
        if (i %in% ycode)
        {
            nm <- x[match(i, ycode)]
            xx <- evalInFlowFrame(channel, frames[[nm]])
            h <- do.call(density, c(list(x=xx), darg))
            n <- length(h$x)
            max.d <- max(h$y)
            xl <- h$x[c(1, 1:n, n)]
            yl <- i + height * c(0, h$y, 0) / max.d
            panel.polygon(x=xl,y=yl, col=col[i], border=border[i],
                          lty=lty[i], lwd=lwd[i], alpha=alpha[i])
            ## add the filterResult if possible
            if(!is.null(filter[[nm]])){    
                bounds <- glpolygon(filter[[nm]], frames[[nm]],
                                    channels=parm,
                                    verbose=FALSE, plot=FALSE)
                for(j in seq_along(bounds)){
                    tb <- bounds[[j]]
                    if(ncol(tb) == 1 && colnames(tb) == parm){
                        sel <- xl >= min(tb) & xl <= max(tb)
                        if(any(sel)){
                            afun <- approxfun(xl, yl)
                            xr <- c(min(tb), seq(min(tb), max(tb), len=100),
                                    max(tb))
                            yr <- c(i, afun(xr[-c(1, length(xr))]), i)
                            grid.polygon(xr, yr, default.units = "native",
                                         gp=gpar)
                        }
                    }
                }
            }
            if (ref)
            {
                panel.abline(h=i,
                             col=reference.line$col,
                             lty=reference.line$lty,
                             lwd=reference.line$lwd,
                             alpha=reference.line$alpha)
            }
        }
    }
}


setMethod("densityplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, xlab,
                   as.table = TRUE, overlap = 0.3,
                   prepanel = prepanel.densityplot.flowset,
                   panel = panel.densityplot.flowset,
                   ...)
      {
          ocall <- sys.call(sys.parent())
          ccall <- match.call(expand.dots = FALSE)
          ccall <- manipulate.call(ocall, ccall)
          pd <- pData(phenoData(data))
          uniq.name <- createUniqueColumnName(pd)
          ## ugly hack to suppress warnings about coercion introducing
          ## NAs (needs to be `undone' inside prepanel and panel
          ## functions):
          pd[[uniq.name]] <- factor(sampleNames(data)) 
          channel <- x[[3]]  
          if (length(channel) == 3)
          {
              channel <- channel[[2]]
              x[[3]][[2]] <- as.name(uniq.name)
          }
          else x[[3]] <- as.name(uniq.name)
          channel.name <- expr2char(channel)
          channel <- as.expression(channel)
          if (missing(xlab)) xlab <- channel.name
          ccall$x <- x
          ccall$data <- pd
          ccall$prepanel <- prepanel
          ccall$panel <- panel
          ccall$frames <- data@frames
          ccall$channel <- channel
          ccall$as.table <- as.table
          ccall$overlap <- overlap
          ccall$xlab <- xlab
          ccall$horizontal <- TRUE
          ccall$lattice.options <-
              list(axis.padding =
                   list(factor = c(0.6, 1 + 2 * overlap)))
          ccall[[1]] <- quote(lattice::bwplot)
          ans <- eval.parent(ccall)
          ans$call <- ocall
          ans
      })
