
logsumexp <- function (x) {
	y = max(x)
	y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}


#diagonal covariance
log.posterior.logitNormal <- function (parameters.psi.k,
								 mu.k,
								 precision.k,
								 Zkg.k, Zk.k,
								 nonzero.idx){
  
  psi_k_centered <- parameters.psi.k-mu.k
  prior <- 0.5* drop(psi_k_centered^2 %*% precision.k) #sign reversed
  likelihood <- Zk.k*logsumexp(parameters.psi.k) - drop(Zkg.k[nonzero.idx] %*% parameters.psi.k[nonzero.idx]) #sign reversed
  ll = prior +  likelihood 
  #print(cor(parameters.psi.k, psi.array[2,1,]))      
  return(ll)
}


log.posterior.logitNormal.grad <- function (parameters.psi.k,
								 	mu.k,
								 	precision.k,
								 	Zkg.k, Zk.k,
								 	nonzero.idx){
   	  
  
  prior.grad <- precision.k * (parameters.psi.k-mu.k) #sign reversed
  
  
  likelihood.grad <-  Zk.k * softmax(parameters.psi.k) - Zkg.k #sign reversed
  ll.grad = prior.grad + likelihood.grad
  
  return(ll.grad)
}


#if.resacle=T, then take the log(softmax(x)) 

get_MAP <-function( X, mu.vec, prec.vec, psi.ini, 
					   opt.control=list(trace=0, maxit= 10000), n.cores=30,if.resacle=T){

	M <- nrow(X)
	cat("psi. current sample ID:")
	opt.psi <- do.call(rbind,mclapply(1:M, FUN=function(m){
							cat(m," ");
							opt.psi.i <-  optim(par= psi.ini[m,], #initialize using the gibbs.psi.i from previous opt
												fn= log.posterior.logitNormal,
												gr= log.posterior.logitNormal.grad,
												method="BFGS",
												control= opt.control,
												mu.k = mu.vec,
												precision.k = prec.vec,
												Zkg.k = X[m,], Zk.k = sum(X[m,]),
												nonzero.idx = X[m,] >0,
												lower = -30, upper = 10)$par
							if(if.resacle) opt.psi.i <- opt.psi.i-logsumexp(opt.psi.i)
							return(opt.psi.i)
						}, mc.cores=n.cores))
	cat("\n") 

	return(opt.psi)
}


#filter.by.cor a logical vector, denoting if the prior mean and variance is estimated by a subset of bulk samples whose correlation with the piror phi is higher than the correlation between the sum and prior phi 
MAP_PhiEnv <- function(ted.obj, 
					   filter.by.cor=T,
					   opt.control=list(trace=0, maxit= 100000), 
					   n.cores){
	
	theta.merged <- ted.obj $res$ first.gibbs.res$theta.merged

	psudeo.min <- ted.obj$ para$psudeo.min
	
	input.phi <- ted.obj $ para$ input.phi.prior
	tum.key <- ted.obj $ para$ tum.key
	
	Znkg.merged <- ted.obj $ res$ first.gibbs.res$Znkg.merged
	Znkg.merged.env <- Znkg.merged[,!dimnames(Znkg.merged)[[2]] %in% tum.key,]
	
	phi.env.n <- array(NA,dim=dim(Znkg.merged.env))
	dimnames(phi.env.n) <- dimnames(Znkg.merged.env)
	
	for (k in 1:dim(Znkg.merged.env)[2]){
		
		cell.type.name.k <- dimnames(Znkg.merged.env)[[2]][k]
		cat("compute cell type",k,":", cell.type.name.k,"\n")
		
		Zng.k <-  Znkg.merged.env[,k,] #N by G matrix
		Zng.k.norm <- norm.to.one(exp.df=Zng.k, psudeo.min= psudeo.min)
		
		# only fit on gene with non-zero expression. For all zero genes, no need to include in the "regress out" step anyways.
		express.genes.id <- apply(Zng.k,2,max)>0
		Zng.k_nonzero <- Zng.k[, express.genes.id]
		Zng.k.norm_nonzero <- Zng.k.norm[, express.genes.id]
		Zng.k.norm_nonzero.log <-  log(Zng.k.norm_nonzero)
		
		#compute prior mean
		prior.sum <-  apply(Zng.k_nonzero,2,sum) # prior.sum= posterior sum from BayesPrism gibbs sampling
		prior.mu <- log(prior.sum/sum(prior.sum)) # estimate mu by pulling all reads, so that samples with higher cell type fractions are more heavily weighted
		
		if(filter.by.cor){
			#only use subset of samples to compute pirior.
			#useful to exclude samples with essentially non-existent cells
			#selection criteria: sample with higher spearman correlation with the posterior sum than the prior (enriching for samples with higher theta)					#if no samples availble, then use the top 20% (unlikely to happen)
			
			input.phi.k <- input.phi[rownames(input.phi)== cell.type.name.k, express.genes.id]
			
			cor.with.posterior.sum <- drop(cor(t(Zng.k_nonzero), prior.sum,method="spearman"))
			cor.with.prior <- drop(cor(t(Zng.k_nonzero), input.phi.k,method="spearman"))
			
			selected.sample.idx <- cor.with.posterior.sum > cor.with.prior
			
			if(sum(selected.sample.idx)<3){
				cat("warning: cell type",
					dimnames(Znkg.merged.env)[[2]][k],
					":fewer than 3 samples have higher correlation with posterior sum than the prior. Use top 20% \n")
				theta.k <- theta.merged[,colnames(theta.merged)== cell.type.name.k]
				selected.sample.idx <- theta.k > quantile(theta.k,0.8)
			} 
			
			#re-exclude zero genes
			express.genes.id <- apply(Zng.k[selected.sample.idx,],2,max)>0
			Zng.k_nonzero <- Zng.k[, express.genes.id]
			Zng.k.norm_nonzero <- Zng.k.norm[, express.genes.id]
			Zng.k.norm_nonzero.log <-  log(Zng.k.norm_nonzero)
		
			#re-compute prior mean
			prior.sum <-  apply(Zng.k_nonzero[selected.sample.idx,],2,sum) # prior.sum= posterior sum from BayesPrism gibbs sampling
			prior.mu <- log(prior.sum/sum(prior.sum))
			#compute prior variance (with psuedo-count, and then take log)
			prior.prec <- 1 / apply(Zng.k.norm_nonzero.log[selected.sample.idx,],2,var)
		}
		else{
			#compute prior variance (with psuedo-count, and then take log)
			prior.prec <- 1 / apply(Zng.k.norm_nonzero.log,2,var)
		}

		#optimize the log posterior
		
		# psi.ini <- matrix(NA,nrow=nrow(Zng.k_nonzero),ncol=ncol(Zng.k_nonzero))
		# for(n in 1:nrow(psi.ini)) psi.ini[n,] <- prior.mu
		# decided to initilaize with MLE estimator (not much different than initilaze with psi.ini)
		phi.env.n[,k, express.genes.id] <- get_MAP ( X= Zng.k_nonzero, 
													 mu.vec= prior.mu,
													 prec.vec= prior.prec,
													 psi.ini= Zng.k.norm_nonzero.log,
													 n.cores= n.cores,
													 opt.control= opt.control)		
	}
	
	ted.obj$res$phi.env.n <- phi.env.n
	ted.obj
}





