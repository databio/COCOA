# PCRSA: Principal Component Region Set Analysis
[![Build Status](https://travis-ci.org/databio/PCRSA.svg?branch=master)](https://travis-ci.org/databio/PCRSA)
A method for annotating principal components of DNA methylation data with region sets. This software is still in beta version so there may be some changes to user interface and function behaviour (as well as documentation improvements) but the core utilities of the package are functional.

PCA of DNA methylation data can be hard to interpret. Our method annotates principal components (PCs) with region sets. 
A region set is a set of genomic regions that share a biological annotation, for instance transcription factor binding regions, regions with a certain histone modification, or chromatin accessibility regions. 
PCRSA can identify meaningful sources of variation for each principal component, giving insight into the biological significance of the PC and into the variation among samples.

## Installing PCRSA
PCRSA may be installed from Github:

devtools::install_github("databio/PCRSA")

or locally after downloading/cloning the source code:

install.packages("path/to/PCRSA/directory", repos=NULL, type="source")

## Learning How to Use PCRSA
A vignette is included with the package that shows how to use the main PCRSA functions [(here)](./vignettes/IntroToPCRSA.Rmd). An additional vignette that shows how to use PCRSA with a region set database as you normally would in an analysis is [(here)](./vignettes/PCRSA_Workflow.Rmd).
