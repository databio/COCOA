# COCOA: Coordinate Covariation Analysis
[![Build Status](https://travis-ci.com/databio/COCOA.svg?branch=master)](https://travis-ci.com/databio/COCOA)
**COCOA is a method for understanding epigenetic variation among samples.** COCOA can be used with epigenetic data that includes genomic coordinates and an epigenetic signal, such as DNA methylation and chromatin accessibility data. To describe the method on a high level, COCOA quantifies inter-sample variation with either a supervised or unsupervised technique then uses a database of "region sets" to annotate the variation among samples. A region set is a set of genomic regions that share a biological annotation, for instance transcription factor (TF) binding regions, histone modification regions, or open chromatin regions. COCOA can identify region sets that are associated with epigenetic variation between samples and increase understanding of variation in your data.


## Installing COCOA
To install from Bioconductor (recommended):

```{r, eval=FALSE, message=FALSE, warning=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("COCOA")
```

COCOA may also be installed from Github:

```{r, eval=FALSE, message=FALSE, warning=FALSE}
devtools::install_github("databio/COCOA")
```

or locally after downloading/cloning the source code:

```{r, eval=FALSE, message=FALSE, warning=FALSE}
install.packages("path/to/COCOA/directory", repos=NULL, type="source")
```

## Learning How to Use COCOA
A vignette is included with the package that shows an [overview of COCOA](http://code.databio.org/COCOA/articles/IntroToCOCOA.html) and walks you through multiple analysis scenarios with code.
