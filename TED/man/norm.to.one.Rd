\name{norm.to.one}
\alias{norm.to.one}
\title{
Utility function to prepare the input.phi

}
\description{
Normalize the scRNA-seq reference matrix such that the sum of expression of each cell type=1, while the min of expression = the specified small number.
}
\usage{
norm.to.one (exp.df, pseudo.min=1E-8)
}
\arguments{
  \item{exp.df}{matrix for the raw expression profiles of scRNA-seq. Each column is a cell type; each row is a gene.}
  \item{pseudo.min}{The desired minimum value after normalization.}
}
\details{
}
\value{
The normalized expression matrix ready to be used as the input.phi for TED.
}
\examples{
}
\keyword{ TED }
