##############################################################################
##                            flowSet methods                               ##
##############################################################################


## Dedicated prepanel function to set up dimensions
prepanel.densityplot.flowset <- 
    function(x, y, darg=list(n=50), frames, channel, channel.name,
             overlap=0.3, ...)
{
    xl <- range(eapply(frames, range, channel.name), finite=TRUE)
    list(xlim=xl + c(-1,1)*0.07*diff(xl))   
}


## Dedicated panel function to do the plotting and add gate boundaries
panel.densityplot.flowset <-
    function(x, y, darg=list(n=50), ref=FALSE, frames, channel,
             overlap = 0.3, channel.name, filter=NULL,
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
        gpar$lty <- "dotted"
        gpar$fill <- rgb(1,1,1,0.7)
    }
    else if(! "lty" %in% names(gpar))
        gpar$lty <- "dotted"
    else if(! "fill" %in% names(gpar))
        gpar$fill <- rgb(1,1,1,0.7)
    class(gpar) <- "gpar"
    for (i in rev(seq_len(ny))){
        if (i %in% ycode)
        {
            nm <- x[match(i, ycode)]
            xx <- evalInFlowFrame(channel, frames[[nm]])
            r <- unlist(range(frames[[nm]], channel.name))
            ## we ignore data that has piled up on the margins
            rl <- r + c(-1,1)*0.06*diff(r)
            pl <- xx<=r[1]
            pr <- xx>=r[2]
            xxt <- xx[!(pl | pr)]
            ## we need a smaller bandwidth than the default and keep it constant
            if(!("bw" %in% names(darg)))
                darg$bw <- diff(r)/50
            h <- do.call(density, c(list(x=xxt), darg))
            n <- length(h$x)
            max.d <- max(h$y)
            xl <- h$x[c(1, 1:n, n)]
            yl <- i + height * c(0, h$y, 0) / max.d
            panel.polygon(x=xl,y=yl, col=col[i], border=NA, alpha=alpha[i])
            ## we indicate piled up data by vertical lines (if > 1%)
            desat <- function(col, by=75)
                rgb(t(pmax(0, col2rgb(col)-by)),max=255)
            lx <- length(xx)
            if(sum(pl) > lx/100)
                panel.lines(rep(rl[1],2), c(i, i+sum(pl)/lx*height),
                            col=desat(col[i]), lwd=3)
            if(sum(pr) > lx/100)
                panel.lines(rep(rl[2],2), c(i, i+sum(pr)/lx*height),
                            col=desat(col[i]), lwd=3)
            ## add the filterResult if possible
            if(!is.null(filter[[nm]])){    
                bounds <- glpolygon(filter[[nm]], frames[[nm]],
                                    channels=parm,
                                    verbose=FALSE, plot=FALSE)
                oo <- options(warn=-1)
                on.exit(options(oo))
                if(!is.na(bounds)){
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
                options(oo)
            }
            panel.lines(x=xl,y=yl, col=border[i], lty=lty[i],lwd=lwd[i])
            panel.lines(rl, rep(i,2), col="black")
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
          pd[[uniq.name]] <- factor(sampleNames(data),
                                    levels=unique(sampleNames(data))) 
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
          ## That is super ugly!!! How do we get to the channel name
          ## from the formula???
          ccall$channel.name <- gsub("^.*\\(`|`\\).*$", "", channel.name)
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
