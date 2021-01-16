get.gibbs.idx <- function(chain.length, burn.in, thinning){
	all.idx <- 1: chain.length
	burned.idx <-  all.idx[-(1: burn.in)]
	thinned.idx <- burned.idx[seq(1,length(burned.idx), by= thinning)]
	thinned.idx
}



merge.gibbs.res <-function(gibbs.theta, gibbs.Znkg, gibbs.Zkg, type.to.subtype.mapping){
	
	print("merge subtypes")
	
	uniq.cell.types <- unique(type.to.subtype.mapping[,"cell.type"])
	
	#merge theta
	theta.merged <- do.call(cbind,lapply(uniq.cell.types, function(cell.type.i) {
		theta.cell.type.i <- gibbs.theta[, type.to.subtype.mapping[,"cell.type"]== cell.type.i,drop=F]
		theta.cell.type.i.merged <- apply(theta.cell.type.i,1,sum)
	}))
	
	colnames(theta.merged) <- uniq.cell.types
	rownames(theta.merged) <- rownames(gibbs.theta)

	percentage.tab<-apply(theta.merged,2,summary)	
	print(round(percentage.tab,3))

	#merge Znkg
	Znkg.merged <- array(0,
						  dim=c(dim(gibbs.Znkg)[1],length(uniq.cell.types),dim(gibbs.Znkg)[3]),
						  dimnames = list(dimnames(gibbs.Znkg)[[1]], 
	 			                         uniq.cell.types,
	 			                         dimnames(gibbs.Znkg)[[3]]))
	 			                   
	for(k in 1:length(uniq.cell.types)){
		Znkg.cell.type.k <- gibbs.Znkg[, type.to.subtype.mapping[,"cell.type"]== uniq.cell.types[k],, drop=F]
		Znkg.merged[,k,] <- apply(Znkg.cell.type.k,c(1,3),sum)
	}
	
	#merge Zkg
	Zkg.merged <- matrix(0,nrow=length(uniq.cell.types),ncol=ncol(gibbs.Zkg))
	Zkg.merged <- apply(Znkg.merged,c(2,3),sum)
	rownames(Zkg.merged) <- uniq.cell.types
	colnames(Zkg.merged) <- colnames(gibbs.Zkg)

	return(list(theta.merged= theta.merged, Znkg.merged= Znkg.merged, Zkg.merged= Zkg.merged))
}



run.gibbs <- function(input.phi, 
					  X,
					  type.to.subtype.mapping,
					  alpha,
					  thinned.idx,
					  n.cores){

	N<-nrow(X)
	G<-ncol(input.phi)
	K.tot <- nrow(input.phi)
		
	#instantiate starting point	
	theta.ini <- matrix(1/K.tot,nrow=N,ncol=K.tot)	
			
	gibbs.res <- draw.sample.gibbs(theta.ini = theta.ini,
				  			 		phi.hat = input.phi, 
				  			 		X = X, 
				  			 		alpha = alpha,
				  			 		thinned.idx = thinned.idx,
				  			 		conditional.idx = NULL,
				  			 		n.cores = n.cores,
				  			 		compute.posterior = F)		
	
	if(!is.null(type.to.subtype.mapping)){
		#perform merging
		gibbs.res.merged <- merge.gibbs.res(gibbs.theta = gibbs.res$gibbs.theta,
											  gibbs.Znkg = gibbs.res$Znkg,
											  gibbs.Zkg = gibbs.res$Zkg,
											  type.to.subtype.mapping = type.to.subtype.mapping)
		
		return(list(gibbs.theta = gibbs.res $ gibbs.theta,
				Znkg = gibbs.res$Znkg,
				theta.merged = gibbs.res.merged $ theta.merged, 
				Znkg.merged = gibbs.res.merged $ Znkg.merged,
				Zkg.merged = gibbs.res.merged $ Zkg.merged))
	}
	else{
		return(list(gibbs.theta = gibbs.res $ gibbs.theta,
				Znkg = gibbs.res$Znkg))
	}
		
}



run.gibbs.individualPhi <- function( phi.tum, 
									 phi.hat.env,
									 X,
									 alpha,
									 thinned.idx,
									 n.cores=40,
									 theta.ini=NULL){
	
	stopifnot(ncol(phi.tum)==ncol(phi.hat.env))	

	N<-nrow(X)
	G<-ncol(phi.tum)
	K.tot <- 1+ nrow(phi.hat.env)
	
	if(is.null(theta.ini)) theta.ini <- matrix(1/K.tot,nrow=N,ncol=K.tot)	
	
	stopifnot(nrow(phi.tum)==nrow(theta.ini))
	
	theta.final <- draw.sample.gibbs.individualPhi(theta.ini= theta.ini,
				  			 		phi.tum = phi.tum, 
				  			 		phi.hat.env = phi.hat.env,
				  			 		X=X, 
				  			 		alpha= alpha,
				  			 		thinned.idx= thinned.idx,
				  			 		n.cores= n.cores)		
	return(theta.final)
}


run.Ted.main <- function(input.phi,
						  input.phi.prior,
						  X,
						  tum.key,
						  type.to.subtype.mapping,
						  alpha,
						  sigma,
						  gibbs.control,
						  opt.control,
						  n.cores,
						  first.gibbs.only){

	thinned.idx <- get.gibbs.idx (chain.length = gibbs.control$chain.length, 
								  burn.in = gibbs.control$burn.in, 
								  thinning = gibbs.control$thinning)
	
	para <- list(X=X,
				 input.phi= input.phi,
				 input.phi.prior= input.phi.prior,
				 tum.key = tum.key,
				 type.to.subtype.mapping = type.to.subtype.mapping, 
				 alpha= alpha,
				 sigma= sigma,
				 gibbs.control= gibbs.control,
				 opt.control= opt.control,
				 n.cores= n.cores)
	
	print("run first sampling")
	first.gibbs.res <- run.gibbs (input.phi= input.phi, 
					  			  X=X,
					  			  type.to.subtype.mapping = type.to.subtype.mapping,
					  			  alpha= alpha,
					  			  thinned.idx = thinned.idx,
					  			  n.cores)
		
	if(first.gibbs.only) return(list(para= para, res= list(first.gibbs.res= first.gibbs.res)))
	
	print("pooling information across samples")
	env.prior.num <- -1 / (2* sigma ^2)
	env.prior.mat <- matrix(env.prior.num, nrow= nrow(input.phi.prior), ncol=ncol(input.phi.prior))
	
	Zkg <- first.gibbs.res $ Zkg.merged[match(rownames(input.phi.prior),rownames(first.gibbs.res $ Zkg.merged)),]
	batch.opt.res <- optimize.psi(input.phi = input.phi.prior,
					   			Zkg = Zkg,
					   			prior.mat = env.prior.mat,
					   			opt.control = opt.control,
					   			n.cores = n.cores)			   			
	phi.env.batch.corrected <- transform.phi.mat(input.phi = input.phi.prior, log.fold = batch.opt.res $opt.psi)

	print("run final sampling")
	if(!is.null(tum.key)){
		Zkg.tum <- first.gibbs.res$Znkg.merged[, tum.key,]	
		Zkg.tum.norm <- norm.to.one(Zkg.tum, min(input.phi))
			
		final.gibbs.theta <- run.gibbs.individualPhi ( phi.tum = Zkg.tum.norm, 
												 phi.hat.env= phi.env.batch.corrected, 
				   								 X=X,
				   								 alpha=1, 
				   								 thinned.idx = thinned.idx,
				   								 n.cores= n.cores)
	
		percentage.tab<-apply(final.gibbs.theta,2,summary)	
		print(round(percentage.tab,3))
	
		res <- list(first.gibbs.res= first.gibbs.res,
				Z.tum.first.gibbs = Zkg.tum,
				Zkg.tum.norm = Zkg.tum.norm,
				phi.env= phi.env.batch.corrected,
				final.gibbs.theta = final.gibbs.theta)
	}
	else{
		final.gibbs.theta <- run.gibbs (input.phi= phi.env.batch.corrected, 
					  			       X=X,
					  			       type.to.subtype.mapping = NULL,
					  			       alpha= alpha,
					  			       thinned.idx = thinned.idx,
					  			       n.cores) $gibbs.theta
		
		percentage.tab<-apply(final.gibbs.theta,2,summary)	
		print(round(percentage.tab,3))
		
		res <- list(first.gibbs.res= first.gibbs.res,
					phi.env= phi.env.batch.corrected,
					final.gibbs.theta = final.gibbs.theta)		
	}
	return(list(para= para, res= res))	
}


get.cormat <- function( Zkg.tum ){
	
	Zkg.tum.centered <- t(scale(t(Zkg.tum),center=T,scale=F))	
	cor.mat <- cor(Zkg.tum.centered) #pearson on centered vst
		
	cor.mat
}


run.Ted <- function(ref.dat, 
				X,
				cell.type.labels,
				cell.subtype.labels=NULL,
				tum.key=NULL,
				input.type,
				psudeo.min=1E-8, 
				alpha=1,
				sigma=2,
				outlier.cut=0.05,
				gibbs.control=list(chain.length=400,burn.in=200,thinning=2),
				opt.control=list(trace=0, maxit= 10000),
				n.cores=1,
				sig.gene=NULL,
				pdf.name=NULL,
				first.gibbs.only=F){
				    	
	#check input data format
	if( is.null(colnames(ref.dat))) stop("Error: please specify the gene names of ref.dat!")
	if( is.null(colnames(X))) stop("Error: please specify the gene names of X!")
	if (is.null(cell.type.labels)) stop("Error: please specify the cell.types.labels!")
	if (nrow(ref.dat) != length(cell.type.labels)) stop("Error: nrow(ref.dat) and length(cell.type.labels) need to be the same!")

	#check cell.type.labels and cell.subtype.labels
	if(is.null(cell.subtype.labels)) cell.subtype.labels <- cell.type.labels
	
	#creat mapping between cell type and phenotype (cell type is a superset of phenotype)
	type.to.subtype.mapping <- unique(cbind(cell.type=cell.type.labels, cell.subtype= cell.subtype.labels))
	if (max(table(type.to.subtype.mapping[,"cell.subtype"]))>1) stop("Error: one or more subtypes belong to multiple cell types!")
	if (length(unique(cell.type.labels)) > length(unique(cell.subtype.labels))) stop("Error: more cell types than subtypes!")


	#rm any genes with non numeric values (including NA values)
	print("removing non-numeric genes...")
	ref.dat <- ref.dat[,apply(ref.dat,2,function(ref.dat.gene.i) as.logical(prod(is.numeric (ref.dat.gene.i),is.finite (ref.dat.gene.i))))]
	X <- X[,apply(X,2,function(X.gene.i) as.logical (prod(is.numeric (X.gene.i), is.finite (X.gene.i)))),drop=F]
	
	#select features
	print("removing outlier genes...")
	X.norm <- apply(X,1,function(vec)vec/sum(vec))
	filter.idx <- apply(X.norm,1,max)<outlier.cut
	X<- X[, filter.idx, drop=F]
	num.genes.filtered <- sum(! filter.idx)
	cat("Number of outlier genes filtered=", num.genes.filtered,"\n")
	
	if(!is.null(sig.gene)) ref.dat <- ref.dat[,colnames(ref.dat) %in% sig.gene] 
	
	if(! input.type %in% c("scRNA","GEP")) stop("Error: please specify the correct input.type!")
	print("aligning reference and mixture...")
	
	if(input.type=="GEP"){
		stopifnot(nrow(ref.dat)==length(cell.type.labels))
		
		processed.dat <- process_GEP (ref= ref.dat,
									  mixture=X, 
									  psudeo.min= psudeo.min,
									  cell.type.labels= cell.type.labels)
	}
	if(input.type=="scRNA"){
		stopifnot(nrow(ref.dat)==length(cell.subtype.labels))
		
		processed.dat <- process_scRNA (ref= ref.dat,
										mixture=X, 
										psudeo.min= psudeo.min, 
										cell.type.labels= cell.type.labels, 
										cell.subtype.labels = cell.subtype.labels)
	}
		
	if(!is.null(tum.key)) {
		tum.idx <- which(rownames(processed.dat$prior.matched.norm)==tum.key)
		if(length(tum.idx)==0) stop("Error: tum.key is not matched to any rownames of the reference, please check the spelling!")
		processed.dat$prior.matched.norm <- processed.dat$prior.matched.norm[-tum.idx,]
	}
	else print("No tumor reference is speficied. Reference profiles are treated equally.")
		
	type.to.subtype.mapping <- type.to.subtype.mapping[match(rownames(processed.dat$ref.matched.norm),type.to.subtype.mapping[,"cell.subtype"]),]
	
	rm(ref.dat,X, cell.type.labels, cell.subtype.labels)
	gc()
			    	
	ted.res <- run.Ted.main (input.phi= processed.dat$ref.matched.norm, 
				  input.phi.prior = processed.dat$prior.matched.norm,
			 	  X= processed.dat$mixture.matched,
			 	  tum.key = tum.key,
			 	  type.to.subtype.mapping = type.to.subtype.mapping,
			 	  alpha=alpha,
			 	  sigma=sigma,
			 	  gibbs.control= gibbs.control,
			 	  opt.control= opt.control,
			 	  n.cores= n.cores,			 	  
			 	  first.gibbs.only= first.gibbs.only)

	#perform vst transformation and make heatmap for tumor 
	if(!is.null(tum.key)){
		Zkg.tum <- ted.res$res$Z.tum.first.gibbs
		#whether possible to compute vst transformation on tumor expression
		#if every gene has at least one zero, vst is not possible, then only export Zkg.tum.norm
		Zkg.tum.round <- t(round(Zkg.tum))
		
		if.vst <- sum(apply(Zkg.tum.round,1,min)==0)< nrow(Zkg.tum.round)
		
		if(if.vst & nrow(Zkg.tum.round)>1) {
			print("vst transformation is feasible")
			Zkg.tum.vst <- vst(Zkg.tum.round, nsub= min(nrow(Zkg.tum.round)/5,1000)) #adjust nsub, to avoid error when too few genes are used
			cor.mat <- get.cormat ( Zkg.tum= Zkg.tum.vst)
		}
		else {
			if(nrow(Zkg.tum.round)==1) print("only one mixture sample. vst transformation is NOT feasible")
			else print("every gene has at least one zero. vst transformation is NOT feasible")
			Zkg.tum.vst <- NULL
			cor.mat <- get.cormat ( Zkg.tum= Zkg.tum.norm )
		}
		
		ted.res$res$Zkg.tum.vst <- Zkg.tum.vst
		ted.res$res$cor.mat <- cor.mat
		
		if(!is.null(pdf.name) & nrow(Zkg.tum.round)>1) plot.heatmap(dat= cor.mat, pdf.name= pdf.name, cluster=T, self=T, show.value=F, metric="is.cor")
	}

	return(ted.res)
	
}



