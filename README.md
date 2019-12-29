TED
========

Tumor micro-Environment Deconvolution (TED): A Fully Bayesian Inference of Tumor Microenvironment composition and gene expression

TED is comprised of the deconvolution modules and the embedding learning module. The deconvolution module leverages cell type-specific expression profiles from scRNA-seq and implements a fully Bayesian inference to jointly estimate the posterior distribution of cell type composition and cell type-specific gene expression from bulk RNA-seq expression of tumor samples. The embedding learning module uses Expectation-maximization (EM) to approximate the tumor expression using a linear combination of tumor pathways while conditional on the inferred expression and fraction of non-tumor cells estimated by the deconvolution module. 


Cite TED:
--------
Bayesian Inference of Cell Composition and Gene Expression Reveals Tumor-Microenvironment Interactions

Tinyi Chu and Charles Danko

--------



Workflow of tfTarget
--------
<img src="img/img1.png">


Requires
--------

* R packages:
	
	DESeq2, parallel, MCMCpack, gplots
	
Installation
--------

* If all dependent packages and commands have been installed, please use the following codes to install/update the package in R terminal. 

```````
library("devtools");
install_github("Danko-Lab/TED/TED")
```````


Usage
----------
library(TED)

use ?function_name for more details

#utility function:
norm.to.one

#R functions:
run.Ted, learn.embedding.withPhiTum, learn.embedding.Kcls

	
Output
----------
use ?function_name for more details




Documents
----------

* R vignette:
 (Coming soon)

* R manual:
 (Coming soon)
