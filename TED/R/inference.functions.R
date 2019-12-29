#utility function


#input.phi is a vector of length G
#log.fold is a vector of length G 

transform.phi <- function (input.phi, log.fold){
		
	stab.constant <- max(log.fold) 		
  	log.fold.stab <- log.fold - stab.constant
  	pert.phi <- input.phi * exp(log.fold.stab)
  	pert.phi <- pert.phi/sum(pert.phi)
  	
  	pert.phi
  	
}


#input.phi is a matrix of length K*G
#log.fold is a vector of length G 

transform.phi.mat <- function (input.phi, log.fold){
		
  	pert.phi.mat<-matrix(NA, nrow=nrow(input.phi), ncol=ncol(input.phi))
  	rownames(pert.phi.mat) <- rownames(input.phi)
  	colnames(pert.phi.mat) <- colnames(input.phi)
  	
  	for(i in 1:nrow(pert.phi.mat)) pert.phi.mat[i,] <- transform.phi (input.phi[i,], log.fold[i,])
  	
  	pert.phi.mat
}





#E step functions
sample.n <- function(X.i, theta.ini.i, phi.hat, alpha, gibbs.idx, conditional.idx, compute.posterior){
		
	G<- ncol(phi.hat)
	K.tot<-nrow(phi.hat)
	
	gibbs.theta.i <- theta.ini.i
	
	gibbs.theta.i.conditional <- theta.ini.i[conditional.idx]
	gibbs.theta.i.free.scale.factor <- 1-sum(gibbs.theta.i.conditional)
	free.idx <- 1: (length(gibbs.theta.i) - length(gibbs.theta.i.conditional))
	
	gibbs.Z.i <- array(NA,c(K.tot,G))
	gibbs.theta.i.sum <- rep(0, length(theta.ini.i))
	Zkg.i <- array(0,c(K.tot,G))
	multinom.coef <- 0
	z.logtheta <- 0
	theta.dirichlet.alpha <- 0
	
	for(idx in 1:max(gibbs.idx)){
	
		#sample Z for patient i
		prob.mat <- phi.hat * gibbs.theta.i #multiply by each column
		
		for (g in 1:G) gibbs.Z.i[,g] <- rmultinom(n=1, size=X.i[g], prob= prob.mat[,g])
						
		# sample theta for patient i
		gibbs.Z.i.k <- apply(gibbs.Z.i,1,sum) #total count for each cell type
		gibbs.theta.i <- rdirichlet(1, gibbs.Z.i.k[free.idx] + alpha)[1,]
		gibbs.theta.i <- c(gibbs.theta.i * gibbs.theta.i.free.scale.factor, gibbs.theta.i.conditional)
				
		if(idx %in% gibbs.idx) {
			Zkg.i <- Zkg.i + gibbs.Z.i
			gibbs.theta.i.sum <- gibbs.theta.i.sum + gibbs.theta.i
			
			if(compute.posterior){
				multinom.coef <-  multinom.coef  - sum(lfactorial(gibbs.Z.i))
				z.logtheta <- z.logtheta + (gibbs.Z.i.k[gibbs.theta.i!=0] %*% log(gibbs.theta.i[gibbs.theta.i!=0])) [1]
				theta.dirichlet.alpha <- theta.dirichlet.alpha + sum( (alpha-1) * log(gibbs.theta.i[gibbs.theta.i!=0]) )
			}
		}
	}
	
	gibbs.samples.num <- length(gibbs.idx)
	
	if(compute.posterior) const <- (multinom.coef + z.logtheta + theta.dirichlet.alpha) / gibbs.samples.num
	else const <- NULL
	
	return(list(Zkg.i = Zkg.i / gibbs.samples.num, 
				gibbs.theta.i = gibbs.theta.i.sum / sum(gibbs.theta.i.sum), 
				const= const))
}
	

draw.sample.gibbs <-function(theta.ini, phi.hat, X, alpha, thinned.idx, conditional.idx, n.cores, compute.posterior){
	
	N<-nrow(X)
	G<-ncol(phi.hat)
	K.tot <-  nrow(phi.hat) 
	

	Zkg.theta <- mclapply(1:N, FUN=function(i){
								   #cat(i," ")
								   sample.n(X.i=X[i,],
								   theta.ini.i = theta.ini[i,], 
								   phi.hat = phi.hat,
								   alpha = alpha,
								   gibbs.idx= thinned.idx,
								   conditional.idx= conditional.idx,
								   compute.posterior= compute.posterior)}, mc.cores=n.cores)
		 
	 gibbs.theta <- array(NA,c(N, K.tot))
	 rownames(gibbs.theta) <- rownames(X)
	 colnames(gibbs.theta) <- rownames(phi.hat)
	 
	 Zkg <- array(0,c(K.tot,G))
	 rownames(Zkg) <- rownames(phi.hat)
	 colnames(Zkg) <- colnames(phi.hat)
	 
	 Znkg <- array(NA,c(N, K.tot, G))
	 
	 const<-0
	 
	 for (n in 1:N) {
	 	Znkg[n,,] <- Zkg.theta[[n]]$Zkg.i
	 	Zkg <- Zkg + Zkg.theta[[n]]$Zkg.i
	 	gibbs.theta[n,] <- Zkg.theta[[n]]$gibbs.theta.i
	 	const <- const + Zkg.theta[[n]]$const
	 }

	return(list(gibbs.theta= gibbs.theta, 
				Znkg = Znkg,
				Zkg = Zkg,
				const = const))

}






draw.sample.gibbs.individualPhi <-function(theta.ini, phi.tum, phi.hat.env, X, alpha, thinned.idx, n.cores){
	
	N <- nrow(X)
	K.tot <-  nrow(phi.hat.env) +1 
		
	Zkg.theta <- mclapply(1:N, FUN=function(i){
								   sample.n(X.i=X[i,],
								   			theta.ini.i = theta.ini[i,], 
								   			phi.hat = rbind(phi.tum[i,], phi.hat.env),
								   			alpha = alpha,
								   			gibbs.idx= thinned.idx,
								   			conditional.idx= NULL,
								   			compute.posterior=F)}, mc.cores=n.cores)
								   					 
	 gibbs.theta <- array(NA,c(N, K.tot))
	 rownames(gibbs.theta) <- rownames(X)
	 colnames(gibbs.theta) <- c("Tumor", rownames(phi.hat.env))
	 	 
	 for (n in 1:N) gibbs.theta[n,] <- Zkg.theta[[n]]$gibbs.theta.i

	return(gibbs.theta)

}












#M step, optimize psi


log.posterior.psi <- function (parameters.psi.k,
								 input.phi,
								 input.phi.log,
								 Zkg.k, Zk.k,
								 prior.vec){
   
  stablizing.constant <- max(parameters.psi.k) 		
  psi.stab <- parameters.psi.k-stablizing.constant
  scale.const <- sum(input.phi * exp(psi.stab)) 		
  pert.phi.log.k <- input.phi.log + psi.stab -log(scale.const)
 
  # should find time to merge the stablizing operations in Rcgmin
 
  log.likelihood <-  (Zkg.k %*% pert.phi.log.k) [1]

  log.prior <- sum(prior.vec * parameters.psi.k ^2) #elementalwise prod
  
  return(-(log.likelihood + log.prior))
}


log.posterior.psi.grad <- function (parameters.psi.k, 
							   		input.phi,
							   		input.phi.log,
							   		Zkg.k, Zk.k,
							   		prior.vec){
   	  
  stablizing.constant <- max(parameters.psi.k) 		
  psi.stab <- parameters.psi.k-stablizing.constant
  pert.phi <- input.phi * exp(psi.stab)
  pert.phi <- pert.phi/sum(pert.phi)

  gradient <- Zkg.k - (Zk.k * pert.phi )  

  log.prior.grad <- 2* prior.vec* parameters.psi.k
  
  log.posterior.grad <-  -(gradient + log.prior.grad)
  
  return(log.posterior.grad)
}






optimize.psi<-function(input.phi,
					   Zkg,
					   prior.mat,
					   opt.control,
					   n.cores){
					   	
	# now consider sigma is a lenght(K)*G matrix, specifying gene and cell type-specific prior distribution. 
	
	#modify this line later, as final version may work with genewise sigma	
	G<-ncol(input.phi)
				 
	Zk <- apply(Zkg,1,sum)

	#optimize psi envir for each patient
		  			  					  		
	opt.res <- mclapply(1:nrow(input.phi), function(idx){
	  			  		Rcgminu(par= rep(0,G),
	  								fn= log.posterior.psi,
	  								gr= log.posterior.psi.grad,
	  			  					control= opt.control, 
	  			  					input.phi = input.phi[idx,],
	  			  					input.phi.log = log(input.phi[idx,]),
	  			  					Zkg.k=Zkg[idx,], Zk.k= Zk[idx],
	  			  					prior.vec = prior.mat[idx,])
	  			  		},mc.cores= n.cores)


	opt.psi <- do.call(rbind,lapply(opt.res, function(res) res$par))
	log.posterior.tum<-unlist(lapply(opt.res, function(res) res$value))
		
	return(list(opt.psi= opt.psi, log.posterior= log.posterior.tum))				   	

}









compute.log.posterior <- function(log.posterior.tum,
								  input.phi.env,
								  Zkg.env,
								  const){
								  	
	#add back the minus likelihood term of environment cells
	log.posterior.env <- - sum(log(input.phi.env) * Zkg.env) #elemental wise prod
	
	#add back the const , note log.posterior is the minus of the actual log of posterior
	log.posterior <- sum(log.posterior.tum, log.posterior.env) - const
	
	log.posterior					  	
}





