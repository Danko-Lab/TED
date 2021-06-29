\name{run.Ted}
\alias{run.Ted}
\title{
Bayesian deconvolution module
}
\description{
Run Bayesian deconvolution to estimate cell type composition and gene expression. Optionally the function can compute the correlation matrix for estimated tumor expression in each sample and plot hirachial clustering and heatmap based on the correlation matrix.
}
\usage{
run.Ted(ref.dat, 
      X,
      cell.type.labels,
      cell.subtype.labels=NULL,
      tum.key=NULL,
      input.type= c("scRNA","GEP"),
      pseudo.min=1E-8,
      alpha=1, 
      sigma=2,
      outlier.cut=0.01,
      outlier.fraction=0.1,
      ref.min.exp=1,
      gibbs.control=list(chain.length=1000,burn.in=500,thinning=2), 
      opt.control=list(trace=0, maxit= 100000),
      file.name=NULL, 
      n.cores=1,
      n.cores.2g=NULL,
      pdf.name=NULL,
      first.gibbs.only=F)
}
\arguments{
  \item{ref.dat}{a matrix containing the raw count of scRNA-seq OR collpased gene expression profile scRNA-seq. Each row is a cell (input.type="scRNA") or a cell type (input.type="GEP"); each column is a gene. rownames of ref.dat are cell IDs (input.type="scRNA") or cell type names(input.type="GEP"); colnames of ref.dat are gene IDs.}
  \item{X}{matrix for the bulk RNA-seq profile. Each row is a sample; each column is a gene. rownames of X are sample IDs of mixtures; colnames of X are gene IDs. Note that run.TED will only use the shared genes between ref.dat and X.}
  \item{cell.type.labels}{a character or factor vector indicating the cell type of each row of the raw count matrix of scRNA-seq or gene expression profile (GEP). The length has to be equal to nrow(ref.dat)}
  \item{cell.subtype.labels}{a character or factor vector indicating the cell subtype of each row of the raw count matrix of scRNA-seq or gene expression profile (GEP). The length needs be equal to nrow(ref.dat). Default is NULL, which uses the same value of cell.type.labels. Note that TED computes the posterior sum over the subtypes to get the total fraction / expression of each cell type. This allows a more fine-grained definition of cell types / cell states.}
  \item{tum.key}{The character in cell.type.labels denoting the tumor cells, e.g. "tumor" or "malignant".}
  \item{input.type}{A character string. ="scRNA" if using the raw count matrix of scRNA-seq as input; ="GEP" if using a collpased gene expression profile, preferably summed raw count(un-normalized), as input.}
  \item{pseudo.min}{A numeric value indicating the minimum (non-zero) value of phi. Default=1E-8.}
  \item{alpha}{One positive numerical parameter or a numeircal vector of length equal nrow(input.phi), denoting the dirichlet hyper-parameter. Default=1, which represents a uniform prior over the simplex of theta. For sparser priors, use 0<alpha<1. Note that alpha usually does not affect the results, due to the dominating likelihood term resulted from the high sequencing depth of bulk RNA-seq.}
  \item{sigma}{One positive numerical parameter or a numeircal vector of length equal number of genes (for gene-specific prior), denoting the prior of the standard deviation of log fold change between the true expression and the reference.Default=2, which represents a weak gene-wise prior. User may provide their own sigma based on prior knowledge, such as differential expression analysis. }
  \item{outlier.cut & outlier.fraction}{Filter genes in X whose expression fraction is greater than outlier.cut (Default=0.01. previous version used 0.05) in more than outlier.fraction (Default=0.1) of bulk data. Typically for dataset with reasonalble quality contol in mapping, very few genes will be filtered. Removal of outlier genes will ensure that the inference will not be dominated by outliers, which sometimes may be resulted from poor QC in mapping.}  
  \item{ref.min.exp}{Keep genes in ref.dat whose maximum expression value across all reference cell types is greater than ref.min.exp. Default=1, i.e. A gene will be included in the deconvolution only if it has at least 2 reads in at least one cell type. Setting ref.min.exp= -1 will include all genes. Filtering out the lowly expressed/unexpressed gene in the reference will make the inference more robust to outliers in the bulk, and will speed up the inference. If using normalized ref.dat, users need to set it to the desired value.}  
  \item{gibbs.control}{A list of parameters controling the Gibbs sampling. Default chain.length=1000, burn.in=500, thinning=2. A list of parameters controling the Gibbs sampling. Default chain.length=1000, burn.in=500, thinning=2. Previous version default is chain.length=400, burn.in=200, thinning=2. Default chain length has been increased to accomondate spatial transcriptomic data which usually has lower depth than conventional bulk data, and hence may need longer chain to reach the stationary distribution.}
  \item{opt.control}{A list of parameters controling the optimization by Rcgmin, Default trace=0, maxit= 100000. }
  \item{file.name}{A file name if wants to write log to files. Default=NULL (print to scree). }
  \item{ncores}{Number of CPU threads used for parallel computing. Default=1}
  \item{n.cores.2g}{Number of CPU threads used for parallel computing for the final Gibbs sampling. Default=NULL (same as ncores). Recommended to set to a number smaller than n.cores, if deconvolving large number of mixtures, such as Visium data, or the number of cell types is large, to avoid memory overflow.}
  \item{pdf.name}{The name of the pdf file containing the correlation heatmap of tumor gene expression to be ploted. Default=NULL (do not plot).}
  \item{first.gibbs.only}{A logical parameter denoting if to only run the first gibbs sampling, i.e. the initial estimates of theta and Z. Default: FALSE}
}
\details{
}
\value{
A list of result is returned in this function. 
  \item{para}{All input data and parameters.}
  \item{res}{All output of TED. }
  \item{res$first.gibbs.res$gibbs.theta}{Initial estimates of fraction for all cell subtypes in each bulk sample.}
  \item{res$first.gibbs.res$Znkg}{Initial estimates of the mean of posterior read count for each cell subtypes  in each bulk sample.}
  \item{res$first.gibbs.res$theta.merged}{Initial estimates of fraction summed across cell types in each bulk sample.}
  \item{res$first.gibbs.res$Znkg.merged}{Initial estimates of the mean of posterior read count summed across cell types in each bulk sample.}
  \item{res$Zkg.tum}{Mean of posterior of gene expression of tumor in each patient.}
  \item{res$Zkg.tum.norm}{Depth normalized Zkg.tum (A pseudo count is added, such that the zero-valued genes have the same value as the min(phi.input)). Refered to as the psi.tum in the TED paper)}
  \item{res$Zkg.tum.vst}{Variance stablized transformed value of Zkg.tum. If vst transformation is not feasible, return NULL. }
  \item{res$phi.env}{Batch effect corrected expression profiles of stromal cells (refered to as the psi.str in the TED paper)}
  \item{res$final.gibbs.theta}{Updated theta after batch correction and tumor expression estimates.}
  \item{res$cor.mat}{The correlation matrix of the estimated tumor expression profiles across bulk RNA-seq samples.}
}
\examples{
}
\keyword{ TED }
