\name{cleanup.genes}
\alias{cleanup.genes}
\title{
Utility function to remove highly expressed outlier genes that are sensitive to batch effects from ref.dat

}
\description{

# This function automatically detects whether ENSEMBLE IDs or gene symbols are used. Recommend to use EMSEMBLE IDs for more unique mapping.
# This function can remove genes from four categories : 1) ribosomal genes, 2) mitochondrial genes, 3) genes in chrM and chrY. 
# It works only for human and mouse genes, by comparing to a list of curated genes from GENCODE and common annotations used by seqC and cell ranger scRNA-seq pipelines. 
# For users' convenience of deconvolving TCGA, the genocodeV22 is used for construction of human database. 
# Mouse GENCODE version used is gencode.vM26. 
# Ribosomal and mitochondrial genes are often highly expressed, providing little signal for cell type deconvolution and are sensitive to batch effects. It is always advisable to remove them. 
# When deconvolving samples from different sexes, genes from chrY and/or chrX are recommend to be removed. Removal of these genes are better done using the exact annotation that generated the scRNA-seq data to avoid the differences between gene annotations. 
# Only apply this function on scRNA-seq reference, but NOT on the bulk mixture. As information from all genes is needed to identify outliers in bulk by run.Ted.

}
\usage{
cleanup.genes(ref.dat, 
              species=c("mm","hs"), 
              gene.type=c("RB","chrM","chrX","chrY"), 
              input.type=c("scRNA","GEP"), 
              exp.cells=1,
              if.toupper=F)
}
\arguments{
  \item{ref.dat}{A matrix containing the raw count of scRNA-seq OR collpased gene expression profile scRNA-seq. Each row is a cell or a cell type; each column is a gene. rownames of ref.dat are cell IDs or cell type names; colnames of ref.dat are gene IDs/names.}
  \item{species}{A character string denoting whether sample is mouse ("mm") or human ("hs").}
  \item{gene.type}{A character vector denoting the category(s) of genes to be removed. "RB"=ribosomal genes; "chrM"=mitochondrial genes; "chrX" and "chrY" denoting genes on the sex chromosome. see description for details.}
  \item{input.type}{A character string. ="scRNA" if using the raw count matrix of scRNA-seq as input; ="GEP" if using a collpased gene expression profile, preferably summed raw count(un-normalized), as input. If ="GEP", exp.cells is set to min(exp.cells,1). }
  \item{exp.cells}{A numeric variable denoting to keep genes with the minimum number of cells (input.type="scRNA") or cell types (input.type="GEP") that show non-zero expression. Default=1, meaning to exclude genes has zero expression across all cells. A typical range is 1 to 10, depending on the number of cells in each cluster. Excluding lowly expressed genes will also speed up deconvolution.}
  \item{if.toupper}{A logical variable denoting if the gene names in the ref.dat has been converted to all uppercases (for seqC output). Default is FALSE.}
  }
\details{
}
\value{
  ref.dat with genes removed from the columns
}
\examples{
}
\keyword{ TED }
