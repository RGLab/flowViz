
#####################################################
#set the global default options for flowViz
########################################################
.flowVizEnv <- new.env()
assign("flowViz.options", list(), envir = .flowVizEnv)


.onLoad <- function(libname, pkgname) 
{
   
	flowViz.options(.defaultFlowVizOptions())
    
}
flowViz.getOption <- function(name)
{
	get("flowViz.options", envir = .flowVizEnv)[[name]]
}


flowViz.options <- function(...)
{
	new <- list(...)
	if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) new <- new[[1]]
	old <- .flowVizEnv$flowViz.options
	## any reason to prefer get("lattice.options", envir = .LatticeEnv)?
	
	## if no args supplied, returns full options list
	if (length(new) == 0) return(old)
	
	nm <- names(new)
	if (is.null(nm)) return(old[unlist(new)]) ## typically getting options, not setting
	isNamed <- nm != "" ## typically all named when setting, but could have mix
	if (any(!isNamed)) nm[!isNamed] <- unlist(new[!isNamed])
	
	## so now everything has non-"" names, but only the isNamed ones should be set
	## everything should be returned, however
	
	retVal <- old[nm]
	names(retVal) <- nm
	nm <- nm[isNamed]
	
	
	.flowVizEnv$flowViz.options <- lattice:::updateList(old, new[nm])
	
	## return changed entries invisibly
	invisible(retVal)
}


.defaultFlowVizOptions <- function()
{
	#modified lattice theme
	myTheme <- modifyList(standard.theme()
							,list(strip.background=list(col=rev(gray(seq(0.3,0.8,length=5))))
								,background=list(col="white")
								)
							)
#	show.settings(myTheme)	

    sillyName <- library
    sillyName(RColorBrewer)
    sillyName(grDevices)

	cR1<-rev(brewer.pal(11, "Spectral"))
	cR2<-IDPcolorRamp(21,t(col2hsv(c("blue","green","yellow","red"))),fr=c(0.7,0))
	list(argcolramp = colorRampPalette(cR1,bias=1)
		,argcolramp_bin = colorRampPalette(cR1,bias=1)
		,argcolramp_flowJo = colorRampPalette(cR2,bias=1)
		,flowVizTheme=myTheme
		)
}
	