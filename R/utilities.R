

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

          ## FIXME: the next section assumes details which may change (but don't know how else to access them)
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
          ans <- norm.center + (chol.cov %*% rbind(x = norm.radius * cos(t), y = norm.radius * sin(t)))
          ans <- as.data.frame(t(ans))
          names(ans) <- c("x", "y")
          ans
      })





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






