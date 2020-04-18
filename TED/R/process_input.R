

#strip the extension number of gene id
strip.gene.id<-function(input.gene.ids){
	unlist(lapply(input.gene.ids, function(gene.id) unlist(strsplit(gene.id,split="\\."))[1]))
}


# this function aligns multiple expression matrices, stored as a list of matrices
align.exp.df<-function(exp.df.list, df.names=NULL){
	
	gene.id.list<-lapply(exp.df.list,function(exp.df)(rownames(exp.df)))
	
	gene.shared<-Reduce(intersect, gene.id.list)
	
	gene.id.list.matched <- lapply(1: length(exp.df.list), function(idx){
		exp.df <- exp.df.list[[idx]]
		gene.id <- gene.id.list[[idx]]
		exp.df[match(gene.shared, gene.id),,drop=F]
	}  )
	
	if(!is.null(df.names)) {
		stopifnot(length(df.names)==length(exp.df.list))
		names(gene.id.list.matched) <- df.names
	} 
	
	gene.id.list.matched
	
}



#sum up read count over the same cell type / patient
collpase.exp.df<-function(exp.df, sample.type.vec){
	
	sample.type.uniq <-unique(sample.type.vec)
	
	exp.df.collapsed<-do.call(cbind,lapply(sample.type.uniq,function(sample.type) apply(exp.df[, sample.type.vec ==sample.type],1,sum)  ))
	
	colnames(exp.df.collapsed)<-sample.type.uniq
	
	exp.df.collapsed
}


#normalize expression vec, s.t. it sum up to one, with the lowest one= psudeo.min

norm.to.one<-function(exp.df, psudeo.min=1E-8){
	
	apply(exp.df,2,function(vec) {
		
		count.tot <- sum(vec)
		gene.num <- length(vec)
		psuedo.count <- psudeo.min * count.tot / (1- gene.num* psudeo.min)
		norm<- vec+ psuedo.count
		norm <- norm/sum(norm)
		norm[vec==0]<-psudeo.min
		norm
	} ) 

}




process_GEP <- function (ref, mixture, psudeo.min){
	
	if( prod(grepl("\\.",colnames(ref))) ) colnames(ref) <- strip.gene.id(colnames(ref))
	if( prod(grepl("\\.",colnames(mixture))) ) colnames(mixture) <- strip.gene.id(colnames(mixture))
	
	aligned.dat <- align.exp.df(exp.df.list=list(t(ref),t(mixture)), df.names=NULL)
	
	ref.matched <- aligned.dat[[1]]
	mixture.matched <- aligned.dat[[2]]
	
	ref.matched.norm <- norm.to.one(exp.df=ref.matched, psudeo.min= psudeo.min)
	
	return(list(ref.matched.norm= t(ref.matched.norm),
				mixture.matched= t(mixture.matched)))
}



process_scRNA <- function (ref, mixture, psudeo.min, pheno.labels){
	
	if( prod(grepl("\\.",colnames(ref))) ) colnames(ref) <- strip.gene.id(colnames(ref))
	if( prod(grepl("\\.",colnames(mixture))) ) colnames(mixture) <- strip.gene.id(colnames(mixture))
	
	ref.collpased <- collpase.exp.df(exp.df=t(ref), sample.type.vec= pheno.labels)
	
	aligned.dat <- align.exp.df(exp.df.list=list(ref.collpased,t(mixture)), df.names=NULL)
	
	ref.matched <- aligned.dat[[1]]
	mixture.matched <- aligned.dat[[2]]
	
	ref.matched.norm <- norm.to.one(exp.df=ref.matched, psudeo.min= psudeo.min)
	
	return(list(ref.matched.norm= t(ref.matched.norm),
				mixture.matched= t(mixture.matched)))
}








