
## things are getting too crowded in latticeMethods.R, so start
## splitting things out into different files.  This one is for
## parallel coordinate plots.


#' @export 
#' @rdname lattice-methods
setMethod("parallel",
          signature(x = "flowFrame", data = "missing"),
          function(x, data, 
                   reorder.by = function(x) var(x, na.rm = TRUE),
                   time = "Time", exclude.time = TRUE,
                   ...)
      {
          expr <- exprs(x)
          column.names <- colnames(expr)
          if (exclude.time) column.names <- column.names[column.names != time]
          if (!is.null(reorder.by))
          {
              column.order <- rev(order(apply(expr[, column.names], 2, reorder.by)))
              column.names <- column.names[column.order]
          }
          parallelplot(expr[, column.names],
                   ...)
      })



## setMethod("parallel",
##           signature(x = "formula", data = "flowSet"),
##           function(x, data, 
##                    reorder.by = function(x) var(x, na.rm = TRUE),
##                    time = "Time", exclude.time = TRUE,
##                    ...)
##       {
##           expr <- exprs(x)
##           column.names <- colnames(expr)
##           if (exclude.time) column.names <- column.names[column.names != time]
##           if (!is.null(reorder.by))
##           {
##               column.order <- rev(order(apply(expr[, column.names], 2, reorder.by)))
##               column.names <- column.names[column.order]
##           }
##           parallel(expr[, column.names], ...)
##       })


## Formula is like  ~ . | a + b. The '.' part is currently ignored.
#' @param filter flowCore filter
#' @export 
#' @rdname lattice-methods
setMethod("parallel",
          signature(x = "formula", data = "flowSet"),
          function(x, data,
                   time = "Time", exclude.time = TRUE,
                   filter = NULL,
                   xlab = NULL, ylab = NULL,
                   ...)
      {

          if(!is.null(filter)){
              if(!is.list(filter)){
                  if(is(filter, "filter"))
                      filter <- filter(data, filter)
              }else if(!is(filter, "filterResultList"))
                  filter <- as(filter, "filterResultList")
          }
          pd <- pData(phenoData(data))
          uniq.name <- createUniqueColumnName(pd)
          ## ugly hack to suppress warnings about coercion introducing NAs
          pd[[uniq.name]] <- factor(sampleNames(data))
          channels <- x[[2]]
          if (length(channels) == 3)
          {
              channels <- channels[[2]]
              x[[2]][[2]] <- as.name(uniq.name)
          }
          else x[[2]] <- as.name(uniq.name)

          ## channels is going to be interpreted as follows: If it is
          ## ".", use all channels except Time.  Otherwise, use
          ## f[,channels], where f is a flow frame


          column.names <- colnames(data)
          if (channels == as.symbol("."))
          {
              if (exclude.time)
                  column.names <- column.names[column.names != time]
          }
          else
              column.names <- eval(channels) ## e.g. c("FSC-H", "Time")

          prepanel.parallel.flowset <- 
              function(x, frames, column.names, ...)
              {
                  x <- as.character(x)
                  if (length(x) > 1) stop("must have only one flow frame per panel")
                  if (length(x) < 1) return()
                  nm <- x
                  z <- as.data.frame(exprs(frames[[nm]])[, column.names])
                  prepanel.default.parallel(z = z)
              }

          panel.parallel.flowset <-
              function(x, 
                       frames, channel, column.names, 
                       ...)
              {
                  x <- as.character(x)
                  if (length(x) > 1) stop("must have only one flow frame per panel")
                  if (length(x) < 1) return()

                  nm <- x
                  z <- as.data.frame(exprs(frames[[nm]])[, column.names])
                  if(!(nm %in% names(filter) || !is(filter[[nm]] ,"filterResult"))){
                      warning("'filter' must either be a filterResultList, a single\n",
                              "filter object or a named list of filter objects.",
                              call.=FALSE)
                      filter <- NULL
                  }
                  if (!is.null(filter))
                  {
                      this.filter.result <- filter[[nm]]
                      groups <- this.filter.result@subSet
                  }
                  else groups <- NULL
                  panel.parallel(1, 1, z = z,
                                 groups = groups,
                                 subscripts = seq_len(nrow(z)),
                                 ...)
              }

          densityplot(x, data = pd,
                      prepanel = prepanel.parallel.flowset,
                      panel = panel.parallel.flowset,
                      frames = data@frames,
                      column.names = column.names,
                      xlab = xlab, ylab = ylab,
                      default.scales =
                      list(x = list(at = c(0, 1), labels = c("Min", "Max")),
                           y =
                           list(alternating = FALSE, axs = "i", tck = 0,
                                at = seq_along(column.names),
                                labels = column.names)),
                      ...)
      })

