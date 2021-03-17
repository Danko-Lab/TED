

#for a given number of clusters, cut the tree and get marker genes between clusters 
get.signature.genes.cls.num <- function(hc, dat.tmp, pval.cut, topN, cls.num, cell.type.labels){
	
	cls <- cutree(hc, cls.num)
	cell.type.cls.labels <- rep(NA,length(cell.type.labels))
	for(cls.i in unique(cls)){
		cell.type.cls.labels[cell.type.labels %in% names(cls)[cls== cls.i]] <- paste("cls", cls.i,sep="")
	}
	
	res <- findMarkers(x= t(dat.tmp), 
					   groups= cell.type.cls.labels, 
					   test.type= "t", 
					   pval.type= "all", 
					   direction="up",
					   sorted=F)
					   
	pval.cut.adj <- pval.cut * (cls.num-1)
	
	sig.gene.ids.cls.num <- unique(unlist(lapply(res,function(res.i) {
			res.i$lfc.min <- apply(res.i[,grepl("logFC",colnames(res.i)),drop=F],1,min)
			res.i <- res.i[res.i$p.value < pval.cut.adj & res.i$lfc.min>0,]
			if(!is.null(topN)) res.i <- res.i[order(res.i$lfc.min,decreasing=T)[1:min(topN,nrow(res.i))],]			
			rownames(res.i)
	} )))
	
	sig.gene.ids.cls.num
}


#return.cell.type.name a logical variable denoting if return the sigature genes separately for each cell type
get.signature.genes <- function(ref.dat,
								cell.type.labels,
								use.hclust,
								return.cell.type.name,
								psudeo.min= 1E-8, 
								pval.cut =0.01,
								topN=100,
								n.cores){
	
	#normalize ref.dat to prepare input for findMarker
	scran.sf <- computeSumFactors(x=t(ref.dat), clusters= cell.type.labels, BPPARAM = MulticoreParam(n.cores))
	dat.tmp <- ref.dat/scran.sf
	dat.tmp <- log2(dat.tmp+0.1) - log2(0.1)

	if(use.hclust){
		ref.collapsed <- collapse.exp.df(exp.df= ref.dat, sample.type.vec= cell.type.labels)
		ref.matched.norm <- norm.to.one(exp.df= ref.collapsed, psudeo.min= psudeo.min)
		ref.matched.norm.log <- log2(ref.matched.norm)

		ref.matched.norm.log.centered <- scale(ref.matched.norm.log,center=T,scale=F)
		cor.mat <- cor(t(ref.matched.norm.log.centered))

		hc <- hclust(as.dist(1-cor.mat),method="ward.D2")
		
		sig.gene.ids <- unique(unlist(mclapply(2:ncol(cor.mat),FUN= get.signature.genes.cls.num, 
												hc= hc, dat.tmp= dat.tmp, 
												pval.cut = pval.cut, topN= topN,
												cell.type.labels= cell.type.labels, 
												mc.cores=n.cores)))
	}
	else{
		res <- findMarkers(x= t(dat.tmp), 
					   groups= cell.type.labels, 
					   test.type= "t", 
					   pval.type= "all", 
					   direction="up",
					   sorted=F, 
					   BPPARAM = MulticoreParam(n.cores))
		
		sig.gene.id.list <- lapply(res,function(res.i) {
			res.i$lfc.min <- apply(res.i[,grepl("logFC",colnames(res.i)),drop=F],1,min)
			res.i <- res.i[res.i$p.value < pval.cut & res.i$lfc.min>0,]
			if(!is.null(topN)) res.i <- res.i[order(res.i$lfc.min,decreasing=T)[1:min(topN,nrow(res.i))],]			
			rownames(res.i)
		} )
		names(sig.gene.id.list) <- names(res)
		
		if(return.cell.type.name) return(sig.gene.id.list)
		else sig.gene.ids <- unique(unlist(sig.gene.id.list))
	}
	
	sig.gene.ids
}








