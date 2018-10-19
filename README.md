# COCOA: Coordinate Covariation Analysis
[![Build Status](https://travis-ci.org/databio/PCRSA.svg?branch=master)](https://travis-ci.org/databio/PCRSA)
**COCOA is a method for understanding variation among samples**. COCOA can be used with data that includes genomic coordinates such as DNA methylation. 
To describe the method on a high level, COCOA uses a database of "region sets" and principal component analysis (PCA) of the data 
to identify sources of variation among samples. A region set is a set of genomic regions that share a biological annotation, 
for instance transcription factor (TF) binding regions, histone modification regions, or open chromatin regions. 
In contrast to some other common techniques, COCOA is unsupervised, meaning that samples do not have to be divided into groups 
such as case/control or healthy/disease, although COCOA works in those situations as well. Also, COCOA focuses on continuous variation 
between samples instead of having cutoffs. Because of this, COCOA can be used as a complementary method alongside "differential" methods 
that find discrete differences between groups of samples and it can also be used in situations where there are no groups.  
COCOA can identify biologically meaningful sources of variation between samples and increase understanding of 
variation in the data. 

So far, the package has been validated on DNA methylation data but we are planning to expand the package to work with genomic range-based data (eg ATAC-seq) as well. The current implementation could theoretically work with any data that has single genomic coordinates, each with an associated value (eg DNA methylation, genetic data).

## Installing COCOA
COCOA may be installed from Github:

devtools::install_github("databio/COCOA")

or locally after downloading/cloning the source code:

install.packages("path/to/COCOA/directory", repos=NULL, type="source")

## Learning How to Use COCOA
A vignette is included with the package that shows [how to use the main COCOA functions](http://code.databio.org/COCOA/articles/IntroToCOCOA.html) and walks you through an example application. An additional vignette shows [how to use COCOA with a region set database](http://code.databio.org/COCOA/articles/COCOA_Workflow.html) as you normally would in an analysis.
