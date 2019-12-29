library(TED)

#run tcga
load("tcga.gbm.example.rdata")

#run deconvolution

tcga.ted <-  run.Ted  (input.phi = t(ref.norm), 
				    X=t(tcga.tumor.pc.NOchrY), 
				    alpha=1E-8, 
				    sigma=2,
				    gibbs.control=list(chain.length=150,burn.in=50,thinning=2),
				    opt.control=list(trace=0, maxit= 100000),
				    n.cores=45,
				    tum.idx=1:60,
				    pdf.name="tcga.tumor")


#run embedding learning

tcga.ebd.res.k4 <- learn.embedding.Kcls  (ted.res = tcga.ted,
									 	  K.vec = 4,
									 	  EM.maxit=50,
									 	  n.cores =25)

