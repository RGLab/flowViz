
setMethod("histogram",
    signature(x = "formula", data = "flowSet"),
    function(x, data, plotType, ...){
      densityplot(x, data, plotType = "histogram", ...)              
    })
setMethod("histogram",
          signature(x="formula", data="flowFrame"),
          function(x, data, ...){
            densityplot(x, data, plotType = "histogram", ...)              
          })

      
      