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




likelihood.classfier <- function(exp.inferred, X.obs){
	
	likelihood.mat <- matrix(NA, nrow=nrow(X.obs), ncol=nrow(exp.inferred))
		
	for( i in 1:nrow(X.obs)){
		for (j in 1:nrow(exp.inferred)){
			exp.inferred.j <- exp.inferred[j,]
			# normalize such that they sum up to one!! not the case for tumor expression profiles
			exp.inferred.j <- exp.inferred.j  / sum(exp.inferred.j) 			
			likelihood.mat[i,j] <- sum( X.obs[i,] * log(exp.inferred.j))
		}
	}
	
	aligned.patient  <- apply(likelihood.mat,1,which.max)
	aligned.patient
}






get.accuracy <- function(gibbs.theta, opt.phi.hat, X.obs){
	
	stopifnot(ncol(gibbs.theta)==nrow(opt.phi.hat))

	exp.inferred<-  gibbs.theta %*% opt.phi.hat
	rownames(exp.inferred) <- rownames(X.obs)
	colnames(exp.inferred) <- colnames(opt.phi.hat)
		
	cor.mat<- cor(t(X.obs), t(exp.inferred),method="spearman")
	aligned.patient <- apply(cor.mat,1,which.max)
	accu.spearman <-  sum(aligned.patient==1:nrow(cor.mat)) / nrow(cor.mat)
	cat( accu.spearman, " " )
	
	aligned.patient <- likelihood.classfier (exp.inferred, X.obs)
	accu.naiveBayes <-  sum(aligned.patient==1:nrow(X.obs)) / nrow(X.obs)
	cat(accu.naiveBayes, " " )
	
	#print("naiveBayes classifier:")
	#print(table(aligned.patient==1:nrow(X.obs)))
	
	return(list(exp.inferred= exp.inferred, accu.spearman = accu.spearman , accu.naiveBayes= accu.naiveBayes ))
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
				Zkg = gibbs.res$Zkg,
				Znk = gibbs.res$Znk))
	
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




run.EM <- function(phi.tum,
				   phi.env, 
				   X,
				   alpha=1, sigma=2, 
				   thinned.idx,
				   theta.ini,
				   opt.control=list(trace=0, maxit= 100000),
				   n.cores=40,
				   EM.maxit=100,
				   EM.res=NULL,
				   X.tum=NULL,
				   print.accuracy=F,
				   compute.posterior=F){
	N<-nrow(X)
	G<-ncol(phi.tum)
	K.tum <- nrow(phi.tum)
	K.tot <- nrow(phi.tum) + nrow(phi.env)
	
	#if genewise sigma is used
	stopifnot(length(sigma)==1 | length(sigma)==G | length(sigma)== K.tum) # 3 ways to specify sigma
	
	if(is.null(nrow(sigma))) sigma <- matrix(sigma,nrow= K.tum,ncol=G, byrow=F)
		
	#instantiate starting point		
	
	if(!is.null(EM.res)) phi.tum <- EM.res$opt.phi.hat.tum
	opt.phi.hat <- rbind(phi.tum, phi.env)
	
	prior.mat <- -1 / (2* sigma ^2)
	
	log.posterior.vec <- Inf
	
	EM.cycle.time.start <- Sys.time()
	
	EM.cycle.idx<-1	
	print("starting EM cycles #...")
	while(EM.cycle.idx <= EM.maxit){ #EM convergence criterion
		cat(EM.cycle.idx, " ")
		#print("running gibbs sampling...")

		gibbs.res<-draw.sample.gibbs(theta.ini= theta.ini,
				  			 		phi.hat= opt.phi.hat, 
				  			 		X=X, 
				  			 		alpha= alpha,
				  			 		thinned.idx= thinned.idx,
				  			 		conditional.idx= (K.tum+1): K.tot,
				  			 		n.cores= n.cores,
				  			 		compute.posterior = compute.posterior)		
		
		
		
		#print("optimizing psi...")
		opt.res <- optimize.psi(input.phi = phi.tum,
					   			Zkg = gibbs.res$Zkg[1:K.tum,],
					   			prior.mat = prior.mat,
					   			opt.control = opt.control,
					   			n.cores = n.cores)
		
		#update the phi.hat
		opt.phi.hat.tum <- transform.phi.mat(input.phi =phi.tum, log.fold =opt.res$opt.psi)
		opt.phi.hat <- rbind(opt.phi.hat.tum, phi.env)
		
		if(compute.posterior){
			log.posterior <- compute.log.posterior (log.posterior.tum= opt.res$log.posterior,
								  					input.phi.env = phi.env,
								  					Zkg.env = gibbs.res$Zkg[-c(1:K.tum),],
								  					const = gibbs.res$const)
			log.posterior.increment <- log.posterior - log.posterior.vec[length(log.posterior.vec)] 
			log.posterior.vec <- c( log.posterior.vec, log.posterior )
			#cat("log.posterior=", opt.res$log.posterior,"\n")
			#cat("log.posterior.increment=", log.posterior.increment,"\n")
		}
				
		EM.cycle.idx<-EM.cycle.idx+1	
	}
	
	cat("\n")
	
	if(print.accuracy){
		print("accuracy at final EM:")			
		exp.all.inferred <- get.accuracy(gibbs.theta=gibbs.res$gibbs.theta, 
											  opt.phi.hat= opt.phi.hat, 
											  X.obs=X)
		
		exp.tum.inferred <- get.accuracy(gibbs.theta = gibbs.res$gibbs.theta[,1:K.tum], 
											  opt.phi.hat= opt.phi.hat[1:K.tum,], 
											  X.obs= X.tum)	
		cat("\n")
	}
	
	#print("EM reaches maximum cycle!")
	
	EM.cycle.time.end <- Sys.time()
	#cat("total time spent on EM", EM.cycle.time.end-EM.cycle.time.start,"minutes ", "\n" )
		
	return(list(K.tum = K.tum,
				theta.all = gibbs.res$gibbs.theta , 
				opt.phi.hat.tum = opt.phi.hat.tum,
				log.posterior = log.posterior.vec))
				# exp.all.inferred = exp.all.inferred$exp.inferred,
				# exp.tum.inferred = exp.tum.inferred$exp.inferred,
				# accu.sp.all = exp.all.inferred $accu.spearman,
				# accu.nb.all = exp.all.inferred $accu.naiveBayes,
				# accu.sp.tum = exp.tum.inferred $accu.spearman,
				# accu.nb.tum = exp.tum.inferred $accu.naiveBayes
					
}






norm.to.one<-function(exp.df, psudeo.min=1E-8){
	
	apply(exp.df,2,function(vec) {
		
		count.tot <- sum(vec)
		gene.num <- length(vec)
		psuedo.count <- psudeo.min * count.tot / (1- gene.num* psudeo.min)
		#print(psuedo.count)
		norm<- vec+ psuedo.count
		norm <- norm/sum(norm)
		norm[vec==0]<-psudeo.min
		norm
	} ) 

}




# get.cormat <- function( Zkg.tum.norm, sig.gene){
		
	# phi.hat.tum.centered <- t(scale(t(Zkg.tum.norm),center=T,scale=F))	
	# if(!is.null(sig.gene)) phi.hat.tum.centered <- phi.hat.tum.centered[rownames(phi.hat.tum.centered) %in% sig.gene, ]
	# cor.mat <- cor(phi.hat.tum.centered,method="spearman")
	
	# cor.mat
	
# }


get.cormat <- function( Zkg.tum, sig.gene){
	
	Zkg.tum.centered <- t(scale(t(Zkg.tum),center=T,scale=F))	
	if(!is.null(sig.gene)) Zkg.tum.centered <- Zkg.tum.centered[rownames(Zkg.tum.centered) %in% sig.gene, ]
	cor.mat <- cor(Zkg.tum.centered) #pearson on centered vst
		
	cor.mat
}



run.Ted <- function(input.phi, 
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

		Zkg.tum <- apply(first.gibbs.res$Znkg[, tum.idx,],c(1,3),sum)
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
	
		if(!is.null(file.name)) plot.heatmap(dat= cor.mat, pdf.name= pdf.name, cluster=T, self=T, show.value=F, metric="is.cor")
	
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





merge.tum.cls <- function(Zkg.tum.norm, cor.mat, K){
	
	dendrogram<- hclust(dist(cor.mat),method="ward.D2")
	cls <- cutree(dendrogram,k=K)
		
	phi.tum <- do.call(rbind, lapply(1:K, function(k) apply(Zkg.tum.norm[,cls==k,drop=F],1,mean) ))
	colnames(phi.tum) <- rownames(Zkg.tum.norm)
	rownames(phi.tum) <- paste("Tumor-", 1:K ,sep="")
	phi.tum
}




learn.embedding.withPhiTum <- function(ted.res,
							phi.tum,
							EM.maxit=50,
							alpha=NULL,
							sigma=NULL,
							file.name=NULL,
							gibbs.control=NULL,
							opt.control=NULL,
							n.cores=NULL,
							print.accuracy=F,
							compute.posterior=F){
	
	Zkg.tum.norm <- ted.res$res$Zkg.tum.norm
	phi.env <- ted.res$res$phi.env
	theta.merged <- ted.res $res$ final.gibbs.theta
	
	X <- ted.res$para$X
	
	if(is.null(alpha)) alpha <- ted.res$para$alpha 
	if(is.null(sigma)) sigma <- ted.res$para$sigma
	if(is.null(gibbs.control)) gibbs.control <- ted.res$para$gibbs.control
	if(is.null(opt.control)) opt.control <- ted.res$para$opt.control
	if(is.null(n.cores)) n.cores <- ted.res$para$n.cores			    
	
	thinned.idx <- get.gibbs.idx (chain.length = gibbs.control$chain.length, 
								  burn.in = gibbs.control$burn.in, 
								  thinning = gibbs.control$thinning)
	
		
	current.K <- nrow(phi.tum)	
	
	cat("K=", current.K ,"\n")
				
	theta.ini<- cbind(do.call(cbind,lapply(1: current.K, function(k) theta.merged[,"Tumor"] / current.K)), theta.merged[,-1])
			
	embed.EM.res <- run.EM(phi.tum = phi.tum, 
								   phi.env = phi.env,
								   X=X,
								   alpha= alpha,
								   sigma= sigma,
								   thinned.idx= thinned.idx,
								   opt.control= opt.control,
								   n.cores= n.cores,
								   EM.maxit= EM.maxit,
								   theta.ini= theta.ini,
								   X.tum= t(Zkg.tum.norm),
								   print.accuracy= print.accuracy,
								   compute.posterior= compute.posterior)
				
	return(embed.EM.res)

}






learn.embedding.Kcls <- function(ted.res,
							K.vec,
							EM.maxit=50,
							alpha=NULL,
							sigma=NULL,
							file.name=NULL,
							gibbs.control=NULL,
							opt.control=NULL,
							n.cores=NULL,
							print.accuracy=F,
							compute.posterior=F){
	
	Zkg.tum.norm <- ted.res$res$Zkg.tum.norm
	phi.env <- ted.res$res$phi.env
	theta.merged <- ted.res $res$ final.gibbs.theta
	
	X <- ted.res$para$X
	
	if(is.null(alpha)) alpha <- ted.res$para$alpha 
	if(is.null(sigma)) sigma <- ted.res$para$sigma
	if(is.null(gibbs.control)) gibbs.control <- ted.res$para$gibbs.control
	if(is.null(opt.control)) opt.control <- ted.res$para$opt.control
	if(is.null(n.cores)) n.cores <- ted.res$para$n.cores			    
	
	thinned.idx <- get.gibbs.idx (chain.length = gibbs.control$chain.length, 
								  burn.in = gibbs.control$burn.in, 
								  thinning = gibbs.control$thinning)
		
	#run clustering of Zkg.tum.norm
	embed.EM.res.list <- list()
	
	for(current.K in K.vec){
			cat("current.K=", current.K ,"\n")
			phi.tum <- merge.tum.cls (Zkg.tum.norm, cor.mat= ted.res$res$cor.mat, current.K)
				
			theta.ini<- cbind(do.call(cbind,lapply(1: current.K, function(k) theta.merged[,"Tumor"] / current.K)), theta.merged[,-1])
			
			embed.EM.res <- run.EM(phi.tum = phi.tum, 
								   phi.env = phi.env,
								   X=X,
								   alpha= alpha,
								   sigma= sigma,
								   thinned.idx= thinned.idx,
								   opt.control= opt.control,
								   n.cores= n.cores,
								   EM.maxit= EM.maxit,
								   theta.ini= theta.ini,
								   X.tum= t(Zkg.tum.norm),
								   print.accuracy= print.accuracy,
								   compute.posterior= compute.posterior)
			
			embed.EM.res.list[[length(embed.EM.res.list)+1]] <- embed.EM.res
			
	}
	
	return(embed.EM.res.list)

}








