#' @export 
#' @rdname densityplot
setMethod("histogram",
    signature(x = "formula", data = "flowSet"),
    function(x, data, plotType, ...){
      densityplot(x, data, plotType = "histogram", ...)              
    })

#' @rdname densityplot
#' @aliases histogram,formula,flowFrame-method
setMethod("histogram",
          signature(x="formula", data="flowFrame"),
          function(x, data, ...){
            densityplot(x, data, plotType = "histogram", ...)              
          })

#override flowSet-version methods to pass data instead of data@frames
#' @aliases histogram,formula,ncdfFlowSet-method
#' @rdname densityplot
#' @export 
setMethod("histogram",
          signature(x = "formula", data = "ncdfFlowSet"),
          function(x, data, ...)
          {
            densityplot(x, data, plotType = "histogram", ...)
            
          })

#' @aliases histogram,formula,ncdfFlowList-method
#' @rdname densityplot
setMethod("histogram",
          signature(x = "formula", data = "ncdfFlowList"),
          function(x, data, ...)
          {
            selectMethod("histogram", signature = c("formula", "ncdfFlowSet"))(x, data, ...)
          })
      