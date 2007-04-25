

### Some utilities (used by the lattice methods, at least)

evalInFlowFrame <- function(expr, envir, enclos = baseenv())
{

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

