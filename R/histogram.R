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

      
      