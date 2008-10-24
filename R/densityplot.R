##############################################################################
##                            flowSet methods                               ##
##############################################################################


## Dedicated prepanel function to set up dimensions
prepanel.densityplot.flowset <- 
    function(x, y, darg=list(n=50), frames, 
             overlap=0.3, subscripts, ...,
             which.channel)
{
    channel.name <- unique(which.channel[subscripts])
    stopifnot(length(channel.name) == 1)
    xl <- range(eapply(frames, range, channel.name), finite=TRUE)
    list(xlim=xl + c(-1,1)*0.07*diff(xl))   
}


##FIXME: enable groups 
## Dedicated panel function to do the plotting and add gate boundaries
panel.densityplot.flowset <-
    function(x, y, darg=list(n=50), frames, channel,
             overlap = 0.3, channel.name, filter=NULL,
             fill=superpose.polygon$col,
             lty=superpose.polygon$lty,
             lwd=superpose.polygon$lwd,
             alpha=superpose.polygon$alpha,
             col=superpose.polygon$border,
             gpar, ...)
{
    which.channel <- tail(which.packet(), 1)
    lc <- length(channel)
    channel <- channel[[which.channel]]
    channel.name <- channel.name[which.channel]
    ycode <- as.numeric(y)
    validName <- !length(grep("\\(", channel.name))
    if (any(duplicated(ycode)))
        warning("Some combinations seem to have multiple samples.  \n  ",
                "Only one will be used.")
    nnm <- as.character(x)
    if(!is.null(filter)){
        if(lc > 1){
            if(!is.list(filter) && is(filter, "filterResultList") ||
               length(filter) != lc)
                stop("You are plotting several channels. Argument 'filter'/n",
                     "must be a list of the same length as number of channels.",
                     call.=FALSE)
            filter <- filter[[which.channel]]
        }
        if(!is.list(filter)){
            if(is(filter, "filter")){
                filter <- lapply(seq_along(nnm), function(x) filter)
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
    superpose.polygon <- trellis.par.get("superpose.polygon")
    border <- rep(col, length = ny)
    col <- rep(fill, length = ny)
    lty <- rep(lty, length = ny)
    lwd <- rep(lwd, length = ny)
    alpha <- rep(alpha, length = ny)
    x <- as.character(x)
    height <- (1 + overlap)
    parm <- gsub("`", "", as.character(channel))
    plotType("gdensity", parm)
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
            lx <- length(xx)
            if(sum(pl) > lx/100)
                panel.lines(rep(rl[1],2), c(i, i+sum(pl)/lx*height),
                            col=desat(col[i]), lwd=3)
            if(sum(pr) > lx/100)
                panel.lines(rep(rl[2],2), c(i, i+sum(pr)/lx*height),
                            col=desat(col[i]), lwd=3)
            ## add the filterResult if possible, we get them from the output of
            ## glpolygon (with plot=FALSE)
            if(!is.null(filter[[nm]]) && validName){    
                bounds <- glpolygon(filter[[nm]], frames[[nm]],
                                    channels=parm,
                                    verbose=FALSE, plot=FALSE)
                oo <- options(warn=-1)
                on.exit(options(oo))
                if(!is.na(bounds)){
                    ## iterate over gate regions
                    for(j in seq_along(bounds)){
                        tb <- bounds[[j]]
                        if(ncol(tb) == 1 && colnames(tb) == parm){
                            sel <- xl >= min(tb) & xl <= max(tb)
                            if(any(sel)){
                                afun <- approxfun(xl, yl)
                                xr <- c(min(tb), seq(min(tb), max(tb), len=100),
                                        max(tb))
                                yr <- c(i, afun(xr[-c(1, length(xr))]), i)
                                panel.polygon(xr, yr, border=gpar$col, col=gpar$fill,
                                              alpha=gpar$alpha, lwd=gpar$lwd,
                                              lty=gpar$lty)
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


analyzeDensityFormula <- function(x, dot.names)
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
        ans$right.comps <- c(list(x), right.comps)
    }
    else if (length(ans$right.symbol) == 1 && as.character(ans$right.symbol) == ".") {
        ans$multiple.right <- TRUE
        ans$right.comps <- lapply(dot.names, as.symbol)
    }
    else ans$right.comps <- list(ans$right.symbol)
    ans
}




setMethod("densityplot",
          signature(x = "formula", data = "flowSet"),
          function(x, data, channels, xlab,
                   as.table = TRUE, overlap = 0.3,
                   prepanel = prepanel.densityplot.flowset,
                   panel = panel.densityplot.flowset,
                   filter=NULL, scales=list(y=list(draw=F)), ...)
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
          if (missing(channels))
              channels <-
                  setdiff(colnames(data),
                          flowCore:::findTimeChannel(data))
          formula.struct <- analyzeDensityFormula(x, dot.names = channels)
          ## we want to add a column to pd for each channel, repeating
          ## pd as necessary.  We might want to skip this if there is
          ## only one channel, but for now we'll use it for
          ## conditioning even then.

          channel.name <-
              sapply(formula.struct$right.comps, expr2char)
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
          if (missing(xlab))
              xlab <- ""
          gp <- list(...)[["par.settings"]]
          gpar <- flowViz.par.get()
          if(!is.null(gp))
              gpar <- lattice:::updateList(gpar, gp)
          ccall$gpar <- gpar$gate.density
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
          ccall$subscripts <- TRUE
          ccall$default.scales <- list(x = list(relation = "free"))
          ccall$which.channel <-
              gsub("^.*\\(`|`\\).*$", "", as.character(pd$which))
          ccall$lattice.options <-
              list(axis.padding =
                   list(factor = c(0.6, 1 + 2 * overlap)))
          ccall[[1]] <- quote(lattice::bwplot)
          ans <- eval.parent(ccall)
          ans$call <- ocall
          ans
      })





## ==========================================================================
## Plot a view object. For everything but gates we have the modified data
## available, hence we can plot directly
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("densityplot",
          signature(x="formula", data="view"),
          function(x, data, ...) densityplot(x, Data(data), ...))


setMethod("densityplot",
          signature(x="view", data="missing"),
          function(x, data, channels, ...)
      {
          dat <- Data(x)
          if(is.null(dat))
              stop("Filter has not been applied to this view.\n",
                   "Run applyParentFilter(view, workflow).")
          if(missing(channels))
              channels <- if(is(x, "normalizeView"))
                  parameters(get(get(x@action)@normalization)) else
              colnames(dat)
          ## A curv1Filter overlay
          filter <- if(is(x, "normalizeView")){
              tmp <- lapply(channels, curv1Filter)
              names(tmp) <- channels
              tmp} else NULL
          ## We construct our own formula
          densityplot(~ ., data=dat, scales=list(y=list(draw=F)),
                      channels=channels, filter=filter)
      })


