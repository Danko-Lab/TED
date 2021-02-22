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
cleanup.genes(ref.dat, species=c("mm","hs"), gene.type=c("RB","chrM","chrX","chrY"), exp.cells=1)
}
\arguments{
  \item{ref.dat}{a matrix containing the raw count of scRNA-seq OR collpased gene expression profile scRNA-seq. Each row is a cell or a cell type; each column is a gene. rownames of ref.dat are cell IDs or cell type names; colnames of ref.dat are gene IDs/names.}
  \item{species}{a character string denoting whether sample is mouse ("mm") or human ("hs").}
  \item{gene.type}{a character vector denoting the category(s) of genes to be removed. "RB"=ribosomal genes; "chrM"=mitochondrial genes; "chrX" and "chrY" denoting genes on the sex chromosome. see description for details.}
  \item{exp.cells}{a numeric variable denoting to keep genes with the minimum number of cells that show non-zero expression. Default=1, meaning to exclude genes has zero expression across all cells. A range to choose is typically from 1-10, depending on the number of cells in each cluster. Excluding lowly expressed genes will also speed up deconvolution.}
  }
\details{
}
\value{
  ref.dat with genes removed from the columns
}
\examples{
}
\keyword{ TED }
