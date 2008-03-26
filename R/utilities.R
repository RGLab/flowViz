

### Some utilities (used by the lattice methods, at least)


## needed for nice labels.  deparse() can potentially give results
## with length > 1
expr2char <- function(x) paste(deparse(x), collapse = "")


evalInFlowFrame <- function(expr, envir, enclos = baseenv())
{
    expr <- as.expression(expr)
    flowFrame2env <- function(ff)
    {
        ## function to convert a flowframe into an environment, so we
        ## can subsequently eval() things in it.  FIXME: defining this
        ## inside hoping for some scope advantage, which may not be
        ## real

        ffdata <- exprs(ff)
        e <- new.env()
        cn <- colnames(ff)
        for (i in seq_along(cn))
        {
            e[[ cn[i] ]] <- ffdata[, i]
        }
        e
    }

    ## FIXME: this copies things, which is potentially bad.  Options
    ## to explore are (1) do thing in C, which may turn out to be not
    ## too bad (2) do things on a limited number of rows at a time

    eval(expr, flowFrame2env(envir), enclos)
}


createUniqueColumnName <- function(x)
{
    make.unique(c(names(x), "sample"))[ncol(x) + 1]
}


## filter.object is an abstract filter (not filter result), and parameters is
## a character vector giving the names of the axes 



setGeneric("filterBoundary",
           function(filter.object, parameters, ...) standardGeneric("filterBoundary"))


## the default is to return nothing

setMethod("filterBoundary",
          signature(filter.object = "filter", parameters = "character"), 
          definition = function(filter.object, parameters, ...)
      {
          list(x = numeric(0), y = numeric(0))
      })


## for a subsetFilter (which is only a binary subset?), draw
## boundaries for components (if possible).  Based on a patch from
## Greg Warnes.  I'm not really sure whether this is general enough
## (it does work for a special case).

setMethod("filterBoundary", 
          signature(filter.object = "subsetFilter", parameters = "character"), 
          definition = function(filter.object, parameters, ...)
      {
          ## see class 'setOperationFilter'
          bnds <- lapply(filter@filters, filterBoundary, ...)
          bnds.x <- unlist(lapply(bnds, function(b) c(b$x, NA)))
          bnds.y <- unlist(lapply(bnds, function(b) c(b$y, NA)))
          list(x = bnds$x, y = bnds$y)
      })



## for norm2Filter, we need the filter result to determine the
## boundary.  If result is not NULL, use it.  Otherwise, apply filter
## to frame to get it.

setMethod("filterBoundary",
          signature(filter.object = "norm2Filter", parameters = "character"), 
          definition =
          function(filter.object, parameters,
                   frame, result = NULL, ...)
      {
          valid <-
              (length(parameters(filter.object)) == 2 &&
               length(parameters) == 2 &&
               setequal(parameters(filter.object), parameters))
          if (!valid)
              return(list(x = numeric(0), y = numeric(0)))
          if (is.null(result)) result <- filter(frame, filter.object)
          result.details <- filterDetails(result)
          if (length(result.details) != 1)
              stop("'result' represents more than one filter.\nThis should not have happened, please send a bug report")

          ## FIXME: the next section assumes details which may change
          ## (but don't know how else to access them)
          result.details <- result.details[[1]]
          if (length(result.details$filter@transformation) > 0) {
              warning("'result' appears to have been applied on transformed data.\nThese are not supported yet.")
              return(list(x = numeric(0), y = numeric(0)))
          }
          norm.center <- result.details$center[parameters]
          norm.cov <- result.details$cov[parameters, parameters]
          norm.radius <- result.details$radius
          ## now what?
          chol.cov <- t(chol(norm.cov))
          t <- seq(0, 2 * base::pi, length = 50)
          ans <- norm.center +
              (chol.cov %*% rbind(x = norm.radius * cos(t),
                                  y = norm.radius * sin(t)))
          ans <- as.data.frame(t(ans))
          names(ans) <- c("x", "y")
          ans
      })



## FIXME: I don't think the next two need to actually apply the filter
## at all...

setMethod("filterBoundary", 
          signature(filter.object = "rectangleGate", parameters = "character"), 
          definition =
          function(filter.object, parameters,
                   frame, result = NULL, ...)
      {
          valid <-
              (length(parameters(filter.object)) == 2 && ## although, 1 should also be OK
               length(parameters) == 2 &&
               setequal(parameters(filter.object), parameters))
          if (!valid)
              return (list(x = numeric(0), y = numeric(0)))
          if (is.null(result)) result <- filter(frame, filter.object)
          result.details <- filterDetails(result)
          if (length(result.details) != 1)
              stop("'result' represents more than one filter.\nThis should not have happened, please send a bug report")

          ## FIXME: the next section assumes details which may change (but don't know how else to access them)
          result.details <- result.details[[1]]
          rect.min <- result.details$filter@min[parameters]
          rect.max <- result.details$filter@max[parameters]
          ans <-
              list(x = c(rect.min[1], rect.max[1], rect.max[1], rect.min[1]),
                   y = c(rect.min[2], rect.min[2], rect.max[2], rect.max[2]))
          ans
      })

setMethod("filterBoundary", 
          signature(filter.object = "polygonGate", parameters = "character"), 
          definition =
          function(filter.object, parameters,
                   frame, result = NULL, ...)
      {
          valid <-
              (length(parameters(filter.object)) == 2 && ## although, 1 should also be OK
               length(parameters) == 2 &&
               setequal(parameters(filter.object), parameters))
          if (!valid)
              return (list(x = numeric(0), y = numeric(0)))
##           if (is.null(result)) result <- filter(frame, filter.object)
##           result.details <- filterDetails(result)
##           if (length(result.details) != 1)
##               stop("'result' represents more than one filter.\nThis should not have happened, please send a bug report")
		  
##           ## FIXME: the next section assumes details which may change (but don't know how else to access them)
##           result.details <- result.details[[1]][[1]]@boundaries
##           ans <- list(x =result.details[,parameters[[1]]], y = result.details[,parameters[[2]]])

          result.details <- filter.object@boundaries
          ans <- list(x = result.details[,parameters[[1]]],
                      y = result.details[,parameters[[2]]])
          ans
      })


setMethod("filterBoundary", 
          signature(filter.object = "curv2Filter", parameters = "character"), 
          definition =
          function(filter.object, parameters,
                   frame, result = NULL, ...)
      {
          valid <-
              (length(parameters(filter.object)) == 2 && 
               length(parameters) == 2 &&
               setequal(parameters(filter.object), parameters))
          if (!valid)
              return (list(x = numeric(0), y = numeric(0)))
          result.details <- attr(result@subSet, "polygons")
          ans <- list(x = unlist(sapply(result.details,
                      function(y) c(y$x, NA))),
                      y = unlist(sapply(result.details,
                      function(y) c(y$y, NA))))
          ans
      })



## mysterious match.call() stuff not understood all that well by the
## author (but bad things happen without them, trust me)...
          

manipulate.call <-
    function(ocall, ccall)
    ## ocall: call actually made by user
    ## ccall: result of match.call(expand.dots = FALSE)
    ## ccall will be modified and returned
{

    ## need to replace things in ccall$"..." with correspondingly
    ## named things in ocall, because language objects get replaced by
    ## ..1, ..2 etc which are never recovered (something to do with
    ## method dispatch)

    if (is.null(ccall$"...")) return(ccall) ## else
    dotnames <- names(ccall$"...")
    if ("" %in% dotnames)
    {
        stop("Unnamed arguments cannot be matched formal arguments")
    }
    ccall$"..." <- NULL
    ccall[dotnames] <- ocall[dotnames]
    ccall
}

