# load libraries and prepare environment for other scripts

library(LOLA)
library(simpleCache)
library(data.table)
library(GenomicRanges)
library(caret)
library(RGenomeUtils)
library(gridExtra) #marrangeGrob for colorClusterPlots()
# some of the environmental variables from aml/.../00-init.R will need to be reset
source(paste0(Sys.getenv("CODE"), "aml_e3999/src/00-init.R" )) 
# the AML init script will set a different plots directory
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/plots/"))

source(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/PRA.R"))
source(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/visualization.R"))


brcaMetadata = fread(paste0(Sys.getenv("CODE"), 
                         "PCARegionAnalysis/metadata/brca_metadata.csv"))
# only keep patients who have definitive status for ER and PGR
brcaMetadata = brcaMetadata[brcaMetadata$ER_status %in% 
                                                     c("Positive", "Negative"), ]
brcaMetadata = brcaMetadata[brcaMetadata$PGR_status %in% 
                                c("Positive", "Negative"), ]
