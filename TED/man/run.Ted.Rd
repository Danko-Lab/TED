\name{run.Ted}
\alias{run.Ted}
\title{
TED deconvolution module
}
\description{
Run Bayesian deconvolution to estimate cell type composition and gene expression. Optionally the function can compute the correlation matrix for estimated tumor expression in each sample and plot hirachial clustering and heatmap based on the correlation matrix.
}
\usage{
run.Ted(input.phi, 
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
      first.gibbs.only=F)
}
\arguments{
  \item{input.phi}{matrix for the reference expression profiles of scRNA-seq. Each row is a cell type; each column is a gene.}
  \item{X}{matrix for the bulk RNA-seq profile. Each row is a sample; each column is a gene.}
  \item{alpha}{One positive numerical parameter or a numeircal vector of length equal nrow(input.phi), denoting the dirichlet hyper-parameter. Default=1, which represents a uniform prior over the simplex of theta. For sparser priors, use 0<alpha<1. Note that alpha usually does not affect the results, due to the dominating likelihood term resulted from the high sequencing depth of bulk RNA-seq.}
  \item{sigma}{One positive numerical parameter or a numeircal vector of length equal number of genes (for gene-specific prior), denoting the prior of the standard deviation of log fold change between the true expression and the reference.Default=2, which represents a weak gene-wise prior. User may provide their own sigma based on prior knowledge, such as differential expression analysis. }
  \item{gibbs.control}{A list of parameters controling the Gibbs sampling. Default chain.length=400, burn.in=200, thinning=2.}
  \item{opt.control}{A list of parameters controling the optimization by Rcgmin, Default trace=0, maxit= 100000. }
  \item{file.name}{A file name if wants to write log to files. Default=NULL (print to scree). }
  \item{ncores}{Number indicating CPU cores in parallel computing. Default=1}
  \item{tum.idx}{A numerical vector indicating which rows of input.phi represent tumor expression profiles.}
  \item{sig.gene}{A character vector indicating the naems of the genes to be used to compute correaltion matrix. Default=NULL (use all genes). }
  \item{pdf.name}{The name of the pdf file containing the correlation heatmap of tumor gene expression to be ploted. Default=NULL (do not plot).}
  \item{first.gibbs.only}{A logical parameter denoting if to only run the first gibbs sampling, i.e. the initial estimates of theta and Z. Default: FALSE}
}
\details{
}
\value{
A list of result is returned in this function. 
  \item{para}{All input data and parameters.}
  \item{res}{All output of TED. }
  \item{res$first.gibbs.res$theta.merged}{Initial estimates of theta0}
  \item{res$first.gibbs.res$Znkg}{Mean of posterior of gene expression in each cell type and each patient.}
  \item{res$first.gibbs.res$Z.tum.first.gibbs}{Mean of posterior of gene expression of tumor in each patient.}
  \item{res$first.gibbs.res$Zkg.tum.norm}{Depth normalized Z.tum.first.gibbs (A pseudo count is added, such that the zero-valued genes have the same value as the min(phi.input)). Refered to as the psi.tum in the TED paper)}
  \item{res$first.gibbs.res$Zkg.tum.vst}{Variance stablized transformed value of Z.tum.first.gibbs}
  \item{res$first.gibbs.res$phi.env}{Batch effect corrected expression profiles of stromal cells (refered to as the psi.str in the TED paper)}
  \item{res$ final.gibbs.theta}{Updated theta after batch correction and tumor expression estimates.}
  \item{res$ cor.mat}{The correlation matrix of the estimated tumor expression profiles across bulk RNA-seq samples.}
}
\examples{
}
\keyword{ TED }
