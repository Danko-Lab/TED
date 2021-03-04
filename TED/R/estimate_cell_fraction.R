# parameters are ordered as cell mean (length=j), alpha(length=i), beta(length=i) and cell var(length=i*j)
# cell var is ordered by row, i.e. i,j={1,1},{1,2},{1,3}...
# , where mu_batch.i = alpha_batch.i + beta_batch.i * mu_batch.base
# X.sum.mat and X2.sum.mat with dim=c(i,j) are sum over X_i,j and X_i,j^2 for batch i and celltype j
# count.tab, with dim=c(i,j), is the total number of cells of batch i and celltype j

log.likelihood <- function (parameters,
							X.sum.mat,
							X2.sum.mat,
						    count.tab,
							base.batch.idx){
  I <- nrow(count.tab)
  J <- ncol(count.tab)
  
  mu <- parameters[1:J]
  alpha <- parameters[(J+1):(J+I)]

  sigma2 <- (X2.sum.mat / count.tab) - (X.sum.mat/count.tab)^2

  log.lik <- 0 # this is minus logLik
  for(i in 1:I){
  	for (j  in 1:J){
  		n.i.j <- count.tab[i,j]
  		if(n.i.j>0){  													 	
  			X.sum.i.j <- X.sum.mat[i,j]
  			X2.sum.i.j <- X2.sum.mat[i,j]
  			mu.transformed.i.j <- alpha[i] + mu[j]
  			sigma2.i.j <- X2.sum.i.j/n.i.j - 2* mu.transformed.i.j * X.sum.i.j/n.i.j + mu.transformed.i.j^2 # always set sigma to MLE

  			#compute minus logLik
  			log.lik <- log.lik + 
  					   0.5* n.i.j * log(sigma2.i.j) + 
  					   (X2.sum.i.j - 2* X.sum.i.j * mu.transformed.i.j + (mu.transformed.i.j^2) * n.i.j) / (2 * sigma2.i.j)

  		}
  	}	
  }
  log.lik
}


log.likelihood.grad <- function (parameters,
								 X.sum.mat,
								 X2.sum.mat,
								 count.tab,
								 base.batch.idx){
  I <- nrow(count.tab)
  J <- ncol(count.tab)
  
  mu <- parameters[1:J]
  alpha <- parameters[(J+1):(J+I)]

  grad.mu <- rep(0,J)
  grad.alpha <- rep(0,I)
  	  
  for(i in 1:I){
 	for (j  in 1:J){
 		n.i.j <- count.tab[i,j]
 		if(n.i.j > 0){
  			X.sum.i.j <- X.sum.mat[i,j]
  			X2.sum.i.j <- X2.sum.mat[i,j]
  			mu.transformed.i.j <- alpha[i] + mu[j]
  			sigma2.i.j <- X2.sum.i.j/n.i.j - 2* mu.transformed.i.j * X.sum.i.j/n.i.j + mu.transformed.i.j^2 # always set sigma to MLE
  		
  			shared.term.i.j <- (- X.sum.i.j + n.i.j* mu.transformed.i.j) / sigma2.i.j
  		
 			grad.mu[j] <- grad.mu[j] +  shared.term.i.j 		
 			grad.alpha[i] <- grad.alpha[i] + shared.term.i.j
 		}
 	}
 	if(i == base.batch.idx) grad.alpha[i] <- 0
  }
  	  	
  return(c(grad.mu, grad.alpha))
}


initial.param <- function(X.sum.mat,
						  X2.sum.mat,
						  count.tab,
						  base.batch.idx){
	
	I <- nrow(count.tab)
  	J <- ncol(count.tab)
  	
	mu.all <- apply(X.sum.mat,2,sum)/apply(count.tab,2,sum)
	mu.base <- X.sum.mat[base.batch.idx,] / count.tab[base.batch.idx,]
	mu.ini <- as.numeric(lm( mu.base - mu.all ~ 1 )$coefficients[1]) + mu.all 
		
	alpha.ini <- c()
	for(i in 1:I){
		if(i== base.batch.idx) {
			alpha.ini <- c(alpha.ini,0)
		}
		else{			
			mu.i <- X.sum.mat[i,] / count.tab[i,]
			lm.i <- lm(mu.i - mu.base ~ 1)
			alpha.i <- as.numeric(lm.i$coefficients[1])
						
			alpha.ini <- c(alpha.ini, alpha.i )			
		}
	}
	
	return(as.numeric(c(mu.ini, alpha.ini)))
}


estimate_sf <- function(ref.dat, 
						cell.type.labels, 
						batch.labels, 
						opt.control=list(trace=1, maxit= 200)){
	
	cell.type.labels <- as.character(cell.type.labels)
	batch.labels <- as.character(batch.labels)
	
	#get total and log2(total) library size for each cell
	cell.tot <- apply(ref.dat,1,sum)
	tot.log <- log2(cell.tot)
	
	#generate batch by cell type count matrix
	count.tab <- as.matrix(table(cbind.data.frame(batch.labels,cell.type.labels)))
	
	#select batch with the most complete cell types (if tie, then order by the total number of cells) as the base batch (scaling factor=1)
	batch.tot.cell.type.count <- apply(count.tab>0,1,sum)
	batch.tot.cell.count <- apply(count.tab,1,sum)
	base.batch.idx <- order(batch.tot.cell.type.count, batch.tot.cell.count, decreasing=T)[1]
		
	I <- nrow(count.tab)
  	J <- ncol(count.tab)
  	
  	X.sum.mat <- matrix(0,nrow=nrow(count.tab),ncol=ncol(count.tab))					  
	X2.sum.mat <- matrix(0,nrow=nrow(count.tab),ncol=ncol(count.tab))					  

  	for(i in 1:I){
 		for (j in 1:J){
 			if(count.tab[i,j] < 2){
 				#clean up entry with fewer than 2 cells 
 				count.tab[i,j] <- 0
 			}
 			else{
				tot.log.i.j <- tot.log[batch.labels ==rownames(count.tab)[i] & 
  									cell.type.labels== colnames(count.tab)[j]]
  				X.sum.mat[i,j] <- sum(tot.log.i.j)					  
				X2.sum.mat[i,j] <- sum(tot.log.i.j^2)
  			}
		}
	}
	
	#initialize parameters
	ini.param <- initial.param (X.sum.mat= X.sum.mat,
						  		X2.sum.mat=X2.sum.mat,
						  		count.tab= count.tab,
						  		base.batch.idx= base.batch.idx)
	
	#minimize minus log likelihood
	opt.param <-  optim(par= ini.param,
						fn= log.likelihood,
						gr= log.likelihood.grad,
						method="BFGS",
						control= opt.control,
						X.sum.mat= X.sum.mat,
						X2.sum.mat=X2.sum.mat,
						count.tab= count.tab,
						base.batch.idx= base.batch.idx)$par
	

	mu.log <- opt.param[1:ncol(count.tab)]
	mu.log <- mu.log - median(mu.log)
	mu <- 2^ mu.log
	
	names(mu) <- colnames(count.tab)
	
	mu
}


convert.theta.mat <- function(theta.mat,
						   	  cell.sf){
	cell.sf.matched <- as.numeric(cell.sf [match(colnames(theta.mat),names(cell.sf))])
	
	t(apply(theta.mat,1,function(theta.mat.i){
		theta.mat.i.converted <- theta.mat.i/cell.sf.matched
		theta.mat.i.converted <- theta.mat.i.converted/sum(theta.mat.i.converted)
		theta.mat.i.converted}
		))
}


convert.cell.fraction <- function(ted.obj,
								  cell.sf){
	first.gibbs.theta <- ted.obj$res$first.gibbs.res$theta.merged
	final.gibbs.theta <- ted.obj$res$final.gibbs.theta
	
	first.gibbs.cell.fraction <- convert.theta.mat (first.gibbs.theta,cell.sf)
	print("first.gibbs.cell.fraction:")
	percentage.tab<-apply(first.gibbs.cell.fraction,2,summary)	
	print(round(percentage.tab,3))
	ted.obj$res$first.gibbs.res$first.gibbs.cell.fraction <- first.gibbs.cell.fraction

	if(!is.null(final.gibbs.theta)){
		final.gibbs.cell.fraction <- convert.theta.mat (final.gibbs.theta,cell.sf)
		print("final.gibbs.cell.fraction:")
		percentage.tab<-apply(final.gibbs.cell.fraction,2,summary)	
		print(round(percentage.tab,3))
		ted.obj$res$final.gibbs.cell.fraction <- final.gibbs.cell.fraction
	}
		
	ted.obj
}
								  
								  
								  
       
