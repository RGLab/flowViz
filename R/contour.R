setGeneric("contour",function(x,...) standardGeneric("contour"))


setMethod("contour",signature("flowFrame"),function(x,y=1:2,col="transparent",border=par("fg"),add=FALSE,
	xlab=names(x)[y[1]],ylab=names(x)[y[2]],nlevels=20,outliers=FALSE,bw=512,grid.size=c(65,65),lwd=1,lty=1,...) {
	if(length(y) != 2)
		stop("You must specify two dimensions!")
	if(is.character(y))
		y = match(y,colnames(x))
	z  = bkde2D(exprs(x)[,y],bw,grid.size)
	ll = contourLines(z$x1,z$x2,z$fhat,nlevels=nlevels)
	dens.range = range(sapply(ll,"[[","level"))

	if(!add) 
		plot(z$x1,z$x2,type="n",xlab=xlab,ylab=ylab,...)
	for(ct in ll) {
		polygon(ct$x,ct$y,border=border,col=col,lwd=lwd,lty=lty)
	}
	
	
	if(outliers) {
		
	}
	
	
	})
	
	
setMethod("contour",signature("flowSet"),function(x,y=1:2,col="transparent",border=par("fg"),lwd=1,lty=1,add=FALSE,...) {
	localPlot = function(...,outliers,bw,nlevels,grid.size) plot(...)
	if(is.character(y))
		y = match(y,colnames(x))
	if(!add) {
		mins = apply(fsApply(x,each_col,min),2,min)
		maxs = apply(fsApply(x,each_col,max),2,max)
		localPlot(c(mins[y[1]],maxs[y[1]]),c(mins[y[2]],maxs[y[2]]),type="n",xlab=colnames(x)[y[1]],ylab=colnames(x)[y[2]],...)		
	}
	col    = rep(col,length.out=length(x))
	border = rep(border,length.out=length(x))
	lwd    = rep(lwd,length.out=length(x))
	lty    = rep(lty,length.out=length(x))
	for(i in 1:length(x)) { contour(x[[i]],y,col=col[i],border=border[i],add=TRUE,...)}
})