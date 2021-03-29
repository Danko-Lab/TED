data_summary <- function(x) {
    m <- median(x)
    ymin <- as.numeric(quantile(x)[2])
    ymax <-as.numeric(quantile(x)[4])
    return(c(y=m,ymin=ymin,ymax=ymax))
 }



plot.theta.distribution <- function(theta.mat, 
								  	bulk.labels=NULL, 
								 	pdf.name=NULL, 
								 	cor.lim=T, 
								  	title="", 
								  	levels.1=NULL, 
								  	levels.2=NULL,
								  	plot.type=c("violin","box")){
	
	stopifnot(plot.type %in% c("violin","box"))
	
	if(is.null(levels.1)) levels.1 <- colnames(theta.mat)
	if(is.null(levels.2)) levels.2 <- unique(bulk.labels)
	
	if(!is.null(bulk.labels)){
		#exclude NA labels
		bulk.labels <- as.character(bulk.labels)
		theta.mat <- theta.mat[!is.na(bulk.labels),,drop=F]
		bulk.labels <- bulk.labels[!is.na(bulk.labels)]
			
		cell.type <- factor(rep(colnames(theta.mat), each=nrow(theta.mat)), levels=levels.1, ordered=T)
		bulk.labels <- factor(rep(bulk.labels, times=ncol(theta.mat)), levels=levels.2, ordered=T)
		theta.vec <- as.vector(as.matrix(theta.mat))
		theta.df <- cbind.data.frame(theta = theta.vec, cell.type = cell.type, bulk.labels = bulk.labels)
		p <- ggplot(theta.df, aes(x= cell.type, y= theta, fill= bulk.labels))	
	}
	else{
		cell.type <- factor(rep(colnames(theta.mat), each=nrow(theta.mat)), levels=levels.1, ordered=T)
		theta.vec <- as.vector(as.matrix(theta.mat))
		theta.df <- cbind.data.frame(theta = theta.vec, cell.type = cell.type)
		p <- ggplot(theta.df, aes(x= cell.type, y= theta, fill= cell.type))
	}
		
	if(plot.type=="violin") p <- p + geom_violin(trim=FALSE, scale="width") + 
									 stat_summary(fun.data= data_summary,geom="pointrange", color="white", position=position_dodge(0.9))
	if(plot.type=="box")  p <- p + geom_boxplot(position=position_dodge(0.9))
	
	if(is.logical(cor.lim)){
		if(cor.lim) p <- p + coord_cartesian(ylim = c(0, 1))
	} 	
	if(is.numeric(cor.lim)){
		p <- p + coord_cartesian(ylim = cor.lim)
	}
	
	p <- p + theme_classic() + ggtitle(title)

	#+ coord_fixed(1.5/1) 
	
	if(!is.null(pdf.name)) pdf(paste(pdf.name,".pdf",sep=""),pointsize=8,useDingbats=FALSE,width=20,height=5)
	print(p)
	if(!is.null(pdf.name)) dev.off()
	p
}







