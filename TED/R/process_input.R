#remove genes with large technical batch effects (ribosomal, mitochondrial, or on sex chromosomes)
#genes curated from seqC, 10x, and gencode annotations
#only mouse and human genes are curated
#for users' convenience, I use human genocodeV22 for compatibility with TCGA
cleanup.genes <- function (ref.dat, 
						   species, 
						   gene.type,
						   input.type, 
						   exp.cells =1,
						   if.toupper=F){
	
	stopifnot(species %in% c("hs","mm"))
	stopifnot(prod(gene.type %in% c("RB","chrX","chrY","chrM"))==1)
	stopifnot(is.logical(if.toupper))
	
	if(! input.type %in% c("scRNA","GEP")) stop("Error: please specify the correct input.type!")
	if(input.type=="GEP"){
		exp.cells <- min(exp.cells,1)
		print("As the input is a collpased GEP, exp.cells is set to min(exp.cells,1)")
	}
	
	#load gene list
	if(species=="hs") gene.list <- read.table(system.file("extdata", "genelist.hs.txt", package="TED"),sep="\t",header=F,stringsAsFactors=F)
	if(species=="mm") gene.list <- read.table(system.file("extdata", "genelist.mm.txt", package="TED"),sep="\t",header=F,stringsAsFactors=F)
	
	gene.list <- gene.list[gene.list[,1] %in% gene.type,]
	
	#detect if EMSEMBLE ID (starts with ENS) or gene symbol is used
	if( sum(substr(colnames(ref.dat),1,3)=="ENS")> ncol(ref.dat)*0.8  ){
		# use 80% of colnames to avoid sometimes there are manually added gene name such as GFP-XXX.
		#use EMSEMBLE ID
		#strip the "." from ENSXXX.X
		print("EMSEMBLE IDs detected. Cleaning up genes based on EMSEMBLE IDs.")
		gene.ids <- unlist(lapply(colnames(ref.dat), function(gene.id) unlist(strsplit(gene.id,split="\\."))[1]))
		exclude.idx <- gene.ids %in% gene.list[,2]
	}
	else{
		#use gene symbols
		print("Gene symbols detected. Cleaning up genes based on gene symbols. Recommend to use EMSEMBLE IDs for more unique mapping.")
		if(if.toupper) gene.list[,3] <- toupper(gene.list[,3])
		exclude.idx <- colnames(ref.dat) %in% gene.list[,3]
	}

	cat("A total of ", sum(exclude.idx)," genes from", gene.type, " have been excluded","\n")
	ref.dat.filtered <- ref.dat[, ! exclude.idx] 
	
	if(exp.cells>0) {
		exclude.lowexp.idx <- apply(ref.dat.filtered>0,2,sum)>= exp.cells
		cat("A total of ", sum(!exclude.lowexp.idx)," lowly expressed genes have been excluded","\n")
		ref.dat.filtered <- ref.dat.filtered[, exclude.lowexp.idx]
	}
	else{
		cat("A total of 0 lowly expressed genes have been excluded","\n")
	}
	ref.dat.filtered
}


# this function aligns multiple expression matrices, stored as a list of matrices
align.exp.df<-function(exp.df.list, df.names=NULL){
	
	gene.id.list<-lapply(exp.df.list, colnames)
	gene.shared<-Reduce(intersect, gene.id.list)
	
	gene.id.list.matched <- lapply(1: length(exp.df.list), function(idx){
		exp.df <- exp.df.list[[idx]]
		gene.id <- gene.id.list[[idx]]
		exp.df[,match(gene.shared, gene.id),drop=F]
	}  )
	
	if(!is.null(df.names)) {
		stopifnot(length(df.names)==length(exp.df.list))
		names(gene.id.list.matched) <- df.names
	} 
	
	gene.id.list.matched
}


#sum up read count over the same cell (sub)type
collapse.exp.df<-function(exp.df, sample.type.vec){
	
	#remove NA in sample.type.vec
	non.na.idx <- !is.na(sample.type.vec)
	if(sum(!non.na.idx)>0) print("Warning: NA found in the cell (sub)type labels. These cells will be excluded!")
	sample.type.vec <- sample.type.vec[non.na.idx]
	exp.df <- exp.df[non.na.idx,]
	
	sample.type.uniq <-unique(sample.type.vec)
	
	exp.df.collapsed<-do.call(rbind,lapply(sample.type.uniq,function(sample.type) apply(exp.df[sample.type.vec ==sample.type, ,drop=F],2,sum)  ))
	
	rownames(exp.df.collapsed)<-sample.type.uniq
	
	exp.df.collapsed
}


#normalize expression vec, s.t. it sum up to one, with the lowest one= psudeo.min
norm.to.one<-function(exp.df, psudeo.min=1E-8){
	
	exp.df.norm <-  do.call(rbind, lapply(1:nrow(exp.df),function(row.idx) {
												vec <- exp.df[row.idx,]
												count.tot <- sum(vec)
												gene.num <- length(vec)
												psuedo.count <- psudeo.min * count.tot / (1- gene.num* psudeo.min)
												norm<- vec+ psuedo.count
												norm <- norm/sum(norm)
												norm[vec==0]<-psudeo.min
												norm
	} )) 

	rownames(exp.df.norm) <- rownames(exp.df)
	colnames(exp.df.norm) <- colnames(exp.df)
	exp.df.norm
}


process_GEP <- function (ref, mixture, psudeo.min, cell.type.labels){
	
	aligned.dat <- align.exp.df(exp.df.list=list(ref,mixture), df.names=NULL)
	
	#this is over subtypes (no need to collapse)
	ref.matched <- aligned.dat[[1]]
	mixture.matched <- aligned.dat[[2]]
	ref.matched.norm <- norm.to.one(exp.df=ref.matched, psudeo.min= psudeo.min)

	#collapse over cell types
	prior.matched <- collapse.exp.df(exp.df= ref.matched, sample.type.vec= cell.type.labels)
	prior.matched.norm <- norm.to.one(exp.df= prior.matched, psudeo.min= psudeo.min)
	
	return(list(ref.matched.norm= ref.matched.norm,
				prior.matched.norm = prior.matched.norm,
				mixture.matched= mixture.matched))
}


process_scRNA <- function (ref, mixture, psudeo.min, cell.type.labels, cell.subtype.labels){
	
	#collapse over subtypes to get reference profile
	ref.collapsed <- collapse.exp.df(exp.df=ref, sample.type.vec= cell.subtype.labels)
	
	#collapse over cell types to get prior profiles
	prior.collapsed <- collapse.exp.df(exp.df=ref, sample.type.vec= cell.type.labels)
	
	aligned.dat <- align.exp.df(exp.df.list=list(ref.collapsed, prior.collapsed, mixture), df.names=NULL)
	
	ref.matched <- aligned.dat[[1]]
	prior.matched <- aligned.dat[[2]]
	mixture.matched <- aligned.dat[[3]]
	
	ref.matched.norm <- norm.to.one(exp.df=ref.matched, psudeo.min= psudeo.min)
	prior.matched.norm <- norm.to.one(exp.df= prior.matched, psudeo.min= psudeo.min)
		
	return(list(ref.matched.norm= ref.matched.norm,
				prior.matched.norm = prior.matched.norm,
				mixture.matched= mixture.matched))
}








