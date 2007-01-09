
parallelCoordinates.rgl <-
    function(x,
             scaled = FALSE,
             alpha = 0.05,
             
             varnames = colnames(x),
             var.at = seq_len(ncol(x)),

             lim = c(min(x), max(x)),
             at = pretty(range(x, na.rm = TRUE)),
             lab = at,
             col = "orange",
             lwd = 0.1,
             add = FALSE,
             ...)
{

    ## ############################
    ##  INPUT
    ## ############################
    ## x       matrix of the data (rows are the observations,
    ##              columns are the variables)
    ## scaled  boolean; if true then scale the data to the range 0,1
    ##          centered about the minimum value
    ## group

    ## ############################
    ## OUTPUT
    ## ############################
    ## plots each row observation across all column variables

    if (!add)
        clear3d()
    
    if (scaled)
    {
        for (i in seq_len(ncol(x)))
        {
            rng <- range(x[, i], na.rm = TRUE)
            if (diff(rng) > 0)
                x[, i] <- (x[, i] - min(rng)) / diff(rng)
        }
    }
    
    zz <- rep(seq_len(nrow(x)) / 1000, each = ncol(x))
    yy <- rep(c(var.at, rev(var.at)), length = length(x))
    xx <- t(x) ## need ``row major'' order w.r.t. x

    alt.rows <- c(FALSE, TRUE)
    xx[, alt.rows] <- xx[rev(seq_len(nrow(xx))), alt.rows]

    savedP <- par3d(skipRedraw = TRUE)
    on.exit(par3d(savedP))

    cCol = material3d("color")

    rgl.linestrips(x = xx,
                   y = yy,
                   z = zz,
                   col = col,
                   lwd = lwd, 
                   alpha = alpha)

    material3d(color = cCol) 
    aspect3d(1,1,0.1)
    ## aspect3d(1,1,0.1)
    rgl.viewpoint(0, 0)

    if (!add)
    {
        rng <- range(xx)
        keep <- at >= rng[1] & at <= rng[2]
        axis3d("x", at = at[keep], lab = lab[keep])
        axis3d("y", at = var.at, lab = varnames)
        title3d(main = "Parallel Coordinates")
    }
}  

setGeneric("parallelCoord", function(object, ...)
    standardGeneric("parallelCoord"))

setMethod("parallelCoord", "flowFrame", function(object) {
    time = match("<Time>", names(object), nomatch=0)
    if( time > 0 )
       parallelCoordinates.rgl(exprs(object)[,-time])
    else
       parallelCoordinates.rgl(exprs(object))
  })

