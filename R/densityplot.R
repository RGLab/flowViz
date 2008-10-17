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
    which.channel <- tail(which.packet(), 1)
    channel <- channel[[which.channel]]
    channel.name <- channel.name[which.channel]
    superpose.line <- trellis.par.get("superpose.line")
    superpose.polygon <- trellis.par.get("superpose.polygon")
    reference.line <- trellis.par.get("reference.line")
    col.signif <- "red"
    ycode <- as.numeric(y)
    validName <- !length(grep("\\(", channel.name))
    if (any(duplicated(ycode)))
        warning("Some combinations seem to have multiple samples.  \n  ",
                "Only one will be used.")
    nnm <- as.character(x)
    if(!is.null(filter)){
        if(!is.list(filter)){
            if(is(filter, "filter")){
                filter <- list(filter)
                names(filter) <- nnm
            }
        }else if(!is(filter, "filterResultList"))
            filter <- as(filter, "filterResultList")
        if(!(nnm %in% names(filter) || !is(filter[[nnm]] ,"filter"))){
            warning("'filter' must either be a filterResultList, a single\n",
                    "filter object or a named list of filter objects.",
                    call.=FALSE)
            filter <- NULL
        }
    }  
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
            if(!is.null(filter[[nm]]) && validName){    
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


analyzeDensityFormula <- function(x)
{
    ans <- list()
    if (length(x) == 2) {
        ans$left <- FALSE
        x <- x[[2]]
    }
    else if (length(x) == 3) {
        ans$left <- TRUE
        ans$left.symbol <- x[[2]]
        x <- x[[3]]
    }
    else stop("unexpected formula structure")
    if (length(x) == 1) {
        ans$conditioned <- FALSE
        ans$multiple.right <- FALSE
        ans$right.symbol <- x
    }
    else if (length(x) == 3) {
        switch(as.character(x[[1]]),
               "+" = {
                   ans$conditioned <- FALSE
                   ans$multiple.right <- TRUE
                   ans$right.symbol <- x
               },
               "|" = {
                   ans$conditioned <- TRUE
                   ans$cond.symbol <- x[[3]]
                   ans$right.symbol <- x[[2]]
                   ans$multiple.right <-
                       (length(x[[2]]) > 1 &&
                        as.character(x[[2]][[1]]) == "+")
               })
    }
    else stop("unexpected formula structure")
    if (ans$multiple.right) {
        x <- ans$right.symbol
        right.comps <- list()
        while ((length(x) > 1) && (as.character(x[[1]]) == "+")) {
            right.comps <- c(list(x[[3]]), right.comps)
            x <- x[[2]]
        }
        ans$right.comps <- rev(c(list(x), right.comps))
    }
    else ans$right.comps <- list(ans$right.symbol)
    ans
}

## str(analyzeDensityFormula(y ~ x1 + I(x2+x3) + x4 | a))


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
          pd[[uniq.name]] <-
              factor(sampleNames(data),
                     levels=unique(sampleNames(data))) 

          formula.struct <- analyzeDensityFormula(x)

          ## we want to add a column to pd for each channel, repeating
          ## pd as necessary.  We might want to skip this if there is
          ## only one channel, but for now we'll use it for
          ## conditioning even then.

          channel.name <-
              ## if (formula.struct$multiple.right) 
              sapply(formula.struct$right.comps, expr2char)
##               else
##                   channel.name <- expr2char(channel)
          pd <- rep(list(pd), length(channel.name))
          names(pd) <- channel.name
          pd <- do.call(lattice::make.groups, pd)
          ## FIXME: this won't work if pd already has a column named
          ## 'which'.  Should deal with that case somehow.

          ## Next task is to manipulate the formula.  The details of
          ## the transformation depends on whether there is a
          ## conditioning variable alread.
          ## y ~ channel ==> y ~ sample | which
          ## y ~ channel | var ==> y ~ sample | which + var
          
          new.x <- d1 ~ d2 | d3
          new.x[[2]] <- ## d1
              if (formula.struct$left) formula.struct$left.symbol
              else as.name("name")
          new.x[[3]][[2]] <- ## d2
              as.name(uniq.name)
          new.x[[3]][[3]] <- ## d3
              if (formula.struct$conditioned) {
                  ans <- (~.+.)[[2]]
                  ans[[3]] <- as.name("which")
                  ans[[2]] <- formula.struct$cond.symbol
                  ## probably not the ideal order, but I don't see how
                  ## to easily ake 'which' the first conditioning
                  ## variable (in case there is more than one
                  ## conditioning variable to begin with)
                  ans
              }
              else as.name("which")

##           channel <- x[[3]]  
##           if (length(channel) == 3)
##           {
##               channel <- channel[[2]]
##               x[[3]][[2]] <- as.name(uniq.name)
##           }
##           else x[[3]] <- as.name(uniq.name)

##           channel <- as.expression(channel)
##           if (missing(xlab)) xlab <- channel.name
##           ccall$x <- x
          if (missing(xlab)) xlab <- ""
          ccall$x <- new.x
          ccall$data <- pd
          ccall$prepanel <- prepanel
          ccall$panel <- panel
          ccall$frames <- data@frames
          ccall$channel <- formula.struct$right.comps ## channel
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
