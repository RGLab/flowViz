
 ##canLog is either length 1, or 
 densityPlot <- function(x, canLog=TRUE, ...) {
  if( !is(x, "flowFrame") )
    stop("only flowFrames")

  cN = colnames(x)
  ##drop Time if it is there
  wT =  match("Time", cN, nomatch=0) 
  cN = cN[-wT]

  numP = length(cN)

  Use = exprs(x)[,cN]
  if(length(canLog) == 1 ) 
    canLog=rep(canLog, numP)
  if( any(canLog) ) {
       doLog = canLog
       for(j in 1:numP) {
          if( canLog[j] )
             doLog[j] = needLog(Use[,j])
       }
  }
          
  densE = vector("list", length=numP)
  names(densE) = cN 
  for(i in 1:numP) 
      densE[[i]] = density(if(doLog[i]) log(Use[,i]) else Use[,i])

  opar = par(mfrow=c(2, ceiling(numP/2)))
  on.exit(par=opar)

  for(j in 1:numP) {
          if( doLog[j] )
              plot(densE[[j]], main=paste("log", cN[j]))
          else
              plot(densE[[j]], main=cN[j])
  }

 }

 needLog = function(x, nbins=10, frac=.1) {
     hN = ceiling(nbins/2)
     vals = table(cut(x, nbins))
     if(max(vals) == vals[1] && all(vals[hN:nbins]/max(vals) < frac ) ) 
         return(TRUE)
     else
       FALSE
  }

##some code from Deepayan, that can be used to get lattice versions

fakedf <- function(x)
{
    cnames <-
        if (is.null(colnames(x))) seq_len(ncol(x))
        else colnames(x)
    data.frame(column = seq_len(ncol(x)),
               colname = factor(cnames, levels = cnames))
}

#densityplot(~ column | colname, data = fakedf(x),
#            src.mat = x,
#            plot.points = FALSE,
#
#            prepanel = function(x, src.mat, ...) {
#                lattice:::prepanel.default.densityplot(src.mat[, x], ...)
#            },
#            panel = function(x, src.mat, ...) {
#                panel.densityplot(src.mat[, x], ...)
#            })


