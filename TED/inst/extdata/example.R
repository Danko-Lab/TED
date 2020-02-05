library(TED)
load(system.file("extdata", "tcga.gbm.example.rdata", package="TED"))


tcga.ted <-  run.Ted  (ref.dat = t(ref.norm), 
				    X=t(tcga.tumor.pc.NOchrY),
				    tum.key="tumor",
				    input.type="GEP",
				    n.cores=45,
				    pdf.name="tcga.tumor")



tcga.ebd.res.k4 <- learn.embedding.Kcls  (ted.res = tcga.ted,
									 	  K.vec = 4,
									 	  EM.maxit=50,
									 	  n.cores =25)

