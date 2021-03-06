

### Some utilities (used by the lattice methods, at least)


## needed for nice labels.  deparse() can potentially give results
## with length > 1
expr2char <- function(x) paste(deparse(x), collapse = "")

##return a formatted stats for display
.getStats <- function(curFilter,curStats,frame, digits, as.is = FALSE, ...){
  popNames<-identifier(curFilter)
  if(is.numeric(curStats)){#stats is explicitly provided
    show.stats <- TRUE
    stats <- curStats
  }else if(is.logical(curStats))#stats needs to be computed on the fly
  {
    
    if(curStats)
    {
      if (!is(curFilter, "filterResult")) 
        curFilter <- filter(frame, curFilter)
      curFres<-curFilter
      if(is(curFres, "multipleFilterResult")){
        ind <- lapply(seq_along(curFres), function(x) as(curFres[[x]], "logical"))
        count <- lapply(ind, sum)                         
        p.stats <- lapply(seq_along(ind), function(x) count[[x]]/length(ind[[x]]))
        stats <-p.stats
        names(stats) <- names(curFres)
      }else{
        ind <- as(curFres, "logical")
        count <- sum(ind)                         
        p.stats <- count/length(ind)
        stats <-p.stats
      }
      show.stats <- TRUE 
    }else
      show.stats <- FALSE
  }
  
  #parse names argument (gate name)
  names<-list(...)$names
  if(is.null(names))
    names<-FALSE

  #format stats when applicable
  if(show.stats)
  {
    if(is.numeric(unlist(stats))){
      if(as.is)
        stats <- as.character(stats)
      else{
        if(is.list(stats)){
          stats <- lapply(names(stats), function(x)
            paste(x,": ", format(stats[[x]]*100,digits=digits),"%",sep=""))
        }else{
          stats <- paste(format(stats*100,digits=digits),"%",sep="")   
        }
      }
    }
    #cat gate name if asked
    if(names)
    {
      stats <- paste(basename(popNames), stats, sep = "\n")
    }
    names(stats)<-rep(popNames, length(stats))
  }else
    stats <- names
  stats
}
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




## mysterious match.call() stuff not understood all that well by the
## author (but bad things happen without them, trust me)...        
manipulate.call <- function(ocall, ccall)
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



## Lighten up or darken colors
desat <- function(col, by=50)
{
    rgbcol <- col2rgb(col)
    up <- max(rgbcol-255)
    down <- min(rgbcol)
    if(by>0)
        rgbcol <- rgbcol + min(by, up)
    if(by<0)
        rgbcol <- rgbcol - min(abs(by), down)
    return(rgb(rgbcol[1,], rgbcol[2,], rgbcol[3,], maxColorValue=255))
}
