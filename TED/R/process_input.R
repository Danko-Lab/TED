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


#sum up read count over the same cell type / patient
collapse.exp.df<-function(exp.df, sample.type.vec){
	
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








