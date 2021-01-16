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



merge.tum.cls <- function(Zkg.tum.norm, cor.mat, K){
	
	dendrogram<- hclust(dist(cor.mat),method="ward.D2")
	cls <- cutree(dendrogram,k=K)
		
	phi.tum <- do.call(rbind, lapply(1:K, function(k) apply(Zkg.tum.norm[cls==k,,drop=F],2,mean) ))
	colnames(phi.tum) <- colnames(Zkg.tum.norm)
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








