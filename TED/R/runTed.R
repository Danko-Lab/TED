get.gibbs.idx <- function(chain.length, burn.in, thinning){
	all.idx <- 1: chain.length
	burned.idx <-  all.idx[-(1: burn.in)]
	thinned.idx <- burned.idx[seq(1,length(burned.idx), by= thinning)]
	thinned.idx
}



print.cell.percentage<-function(X, phi, gibbs.theta, tum.idx){
	
	print("inferred cell percentage")
	
	theta.tum <- gibbs.theta[,tum.idx,drop=F]
	theta.tum.merged <-   apply(theta.tum,1,sum)
	theta.env <- gibbs.theta[,-tum.idx,drop=F]
	
	theta.merged <- cbind(theta.tum.merged, theta.env)
	colnames(theta.merged) <- c("Tumor", rownames(phi)[-tum.idx])
	rownames(theta.merged) <- rownames(X)

	percentage.tab<-apply(theta.merged,2,summary)	
	print(round(percentage.tab,3))

	return(list(theta.merged = theta.merged, theta.tum= theta.tum, theta.env= theta.env))
}



run.gibbs <- function(input.phi, 
					  X,
					  alpha,
					  thinned.idx,
					  n.cores,
					  tum.idx){
	

	N<-nrow(X)
	G<-ncol(input.phi)
	K.tot <- nrow(input.phi)
		
	#instantiate starting point
	#print("instantiate variables")
	
	theta.ini <- matrix(1/K.tot,nrow=N,ncol=K.tot)	
			
	gibbs.res <- draw.sample.gibbs(theta.ini= theta.ini,
				  			 		phi.hat= input.phi, 
				  			 		X=X, 
				  			 		alpha= alpha,
				  			 		thinned.idx= thinned.idx,
				  			 		conditional.idx=NULL,
				  			 		n.cores= n.cores,
				  			 		compute.posterior=F)		
	if(!is.null(tum.idx)){
		theta.list <- print.cell.percentage(X, input.phi, gibbs.res$gibbs.theta, tum.idx)
	}
	else{
		theta.list <- list(theta.merged = gibbs.res$gibbs.theta, theta.tum= NULL)
	}
			
	return(list(theta.merged = theta.list$theta.merged, 
				theta.tum = theta.list$theta.tum, 
				Znkg = gibbs.res$Znkg,
				Zkg = gibbs.res$Zkg))
	
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


get.cormat <- function( Zkg.tum, sig.gene){
	
	Zkg.tum.centered <- t(scale(t(Zkg.tum),center=T,scale=F))	
	if(!is.null(sig.gene)) Zkg.tum.centered <- Zkg.tum.centered[rownames(Zkg.tum.centered) %in% sig.gene, ]
	cor.mat <- cor(Zkg.tum.centered) #pearson on centered vst
		
	cor.mat
}



run.Ted.main <- function(input.phi, 
				    X, 
				    alpha=1, 
				    sigma=2,
				    gibbs.control=list(chain.length=400,burn.in=200,thinning=2),
				    opt.control=list(trace=0, maxit= 100000),
				    file.name=NULL,
				    n.cores=1,
				    tum.idx,
				    sig.gene=NULL, 
				    pdf.name=NULL,
				    first.gibbs.only=F){

		
	if(!is.null(file.name)) sink(file=paste(file.name,".txt",sep=""))

	thinned.idx <- get.gibbs.idx (chain.length = gibbs.control$chain.length, 
								  burn.in = gibbs.control$burn.in, 
								  thinning = gibbs.control$thinning)
	
	para <- list(X=X, 
				 alpha= alpha,
				 sigma= sigma,
				 gibbs.control= gibbs.control,
				 opt.control= opt.control,
				 n.cores= n.cores,
				 sig.gene= sig.gene)
	
	
	print("run first sampling")
		
	first.gibbs.res <- run.gibbs (input.phi= input.phi, 
					  			  X=X,
					  			  alpha= alpha,
					  			  thinned.idx = thinned.idx,
					  			  n.cores,
					  			  tum.idx)
		
	if(first.gibbs.only) return(list(para= para, res= list(first.gibbs.res= first.gibbs.res)))
	
	print("correct batch effect")
	env.prior.num <- -1 / (2* sigma ^2)
	env.prior.mat <- matrix(env.prior.num, nrow= nrow(input.phi)-length(tum.idx), ncol=ncol(input.phi))
	
	if(!is.null(tum.idx)){
		
		batch.opt.res <- optimize.psi(input.phi = input.phi[-tum.idx,],
					   			Zkg = first.gibbs.res $Zkg[-tum.idx,],
					   			prior.mat = env.prior.mat,
					   			opt.control = opt.control,
					   			n.cores = n.cores)
					   			
		phi.env.batch.corrected <- transform.phi.mat(input.phi = input.phi[-tum.idx,], log.fold = batch.opt.res $opt.psi)

		Zkg.tum <- apply(first.gibbs.res$Znkg[, tum.idx,,drop=F],c(1,3),sum)
		colnames(Zkg.tum) <- colnames(X)
		rownames(Zkg.tum) <- rownames(X)
	
		Zkg.tum.norm <- norm.to.one(t(Zkg.tum), min(input.phi))
		
		#whether possible to compute vst transformation on tumor expression
		#if every gene has at least one zero, vst is not possible, then only export Zkg.tum.norm
		Zkg.tum.round <- t(round(Zkg.tum))
		if.vst <- sum(apply(Zkg.tum.round,1,min)==0)< nrow(Zkg.tum.round)
		
		if(if.vst) {
			print("vst transformation is feasible")
			Zkg.tum.vst <- vst(Zkg.tum.round)
			cor.mat <- get.cormat ( Zkg.tum= Zkg.tum.vst, sig.gene= sig.gene)
		}
		else {
			print("every gene has at least one zero. vst transformation is NOT feasible")
			Zkg.tum.vst <- NULL
			cor.mat <- get.cormat ( Zkg.tum= Zkg.tum.norm, sig.gene= sig.gene)
		}
		
		print("run final sampling")
		final.gibbs.theta <- run.gibbs.individualPhi ( phi.tum = t(Zkg.tum.norm), 
												 phi.hat.env= phi.env.batch.corrected, 
				   								 X=X,
				   								 alpha=1, 
				   								 thinned.idx = thinned.idx,
				   								 n.cores= n.cores)
	
		percentage.tab<-apply(final.gibbs.theta,2,summary)	
		print(round(percentage.tab,3))
	
		if(!is.null(pdf.name)) plot.heatmap(dat= cor.mat, pdf.name= pdf.name, cluster=T, self=T, show.value=F, metric="is.cor")
	
		if(!is.null(file.name)) sink()
		
		res <- list(first.gibbs.res= first.gibbs.res,
				Z.tum.first.gibbs = Zkg.tum,
				Zkg.tum.norm = Zkg.tum.norm,
				Zkg.tum.vst = Zkg.tum.vst,
				phi.env= phi.env.batch.corrected,
				final.gibbs.theta = final.gibbs.theta,
				cor.mat= cor.mat)
	
		return(list(para= para, res= res))	
	}
	
	else{
		batch.opt.res <- optimize.psi(input.phi = input.phi,
					   			Zkg = first.gibbs.res $Zkg,
					   			prior.mat = env.prior.mat,
					   			opt.control = opt.control,
					   			n.cores = n.cores)
					   			
		phi.env.batch.corrected <- transform.phi.mat(input.phi = input.phi, log.fold = batch.opt.res $opt.psi)
		
		print("run final sampling")

		final.gibbs.res <- run.gibbs (input.phi= phi.env.batch.corrected, 
					  			  X=X,
					  			  alpha= alpha,
					  			  thinned.idx = thinned.idx,
					  			  n.cores,
					  			  tum.idx)
	
		res <- list(first.gibbs.res= first.gibbs.res,
					phi.env= phi.env.batch.corrected,
					final.gibbs.theta = final.gibbs.res$theta.merged)

		return(list(para= para, res= res))	
		
	}
	
}




run.Ted <- function(ref.dat, 
				X,
				pheno.labels=NULL,
				tum.key=NULL,
				input.type,
				psudeo.min=1E-8, 
				alpha=1,
				sigma=2,
				outlier.cut=0.05,
				gibbs.control=list(chain.length=400,burn.in=200,thinning=2),
				opt.control=list(trace=0, maxit= 100000),
				file.name=NULL,
				n.cores=1,
				sig.gene=NULL,
				pdf.name=NULL,
				first.gibbs.only=F){
				    	
	#check input data format
	if( is.null(colnames(ref.dat))) stop("Error: please specify the gene names of ref.dat!")
	if( is.null(colnames(X))) stop("Error: please specify the gene names of X!")

	#rm any genes with non numeric values (including NA values)
	print("removing non-numeric genes...")
	ref.dat <- ref.dat[,apply(ref.dat,2,function(ref.dat.gene.i) as.logical(prod(is.numeric (ref.dat.gene.i),is.finite (ref.dat.gene.i))))]
	X <- X[,apply(X,2,function(X.gene.i) as.logical (prod(is.numeric (X.gene.i), is.finite (X.gene.i))))]
	
	
	print("removing outlier genes...")
	X.norm <- apply(X,1,function(vec)vec/sum(vec))
	filter.idx <- apply(X.norm,1,max)<outlier.cut
	X<- X[, filter.idx]
	num.genes.filtered <- sum(! filter.idx)
	cat("Number of outlier genes filtered=", num.genes.filtered,"\n")
	

	if(! input.type %in% c("scRNA","GEP")) stop("Error: please specify the correct input.type!")
	
	print("aligning reference and mixture...")
	
	if(input.type=="GEP"){
		processed.dat <- process_GEP (ref= ref.dat,
									  mixture=X, 
									  psudeo.min= psudeo.min)
	}
	if(input.type=="scRNA"){
		stopifnot(nrow(ref.dat)==length(pheno.labels))
		stopifnot(!is.null(tum.key))
		
		processed.dat <- process_scRNA (ref= ref.dat,
										mixture=X, 
										psudeo.min= psudeo.min, 
										pheno.labels= as.character(pheno.labels))
	}
	
	if(!is.null(tum.key)) {
		tum.idx <- which(grepl(tum.key,rownames(processed.dat$ref.matched.norm)))
		if(length(tum.idx)==0) stop("Error: tum.key is not matched to any rownames of the reference, please check the spelling!")
	}
	else {
		tum.idx <- NULL
		print("No tumor reference is speficied. Reference profiles are treated equally.")
	}
	
	
	rm(ref.dat,X, pheno.labels)
	gc()
			    	
	run.Ted.main (input.phi= processed.dat$ref.matched.norm, 
			 	  X= processed.dat$mixture.matched,
			 	  alpha=alpha,
			 	  sigma=sigma,
			 	  gibbs.control= gibbs.control,
			 	  opt.control= opt.control,
			 	  file.name= file.name,
			 	  n.cores= n.cores,
			 	  tum.idx= tum.idx,
			 	  sig.gene= sig.gene,
			 	  pdf.name= pdf.name,
			 	  first.gibbs.only= first.gibbs.only)

}
