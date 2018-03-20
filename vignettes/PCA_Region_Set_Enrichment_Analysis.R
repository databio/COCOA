
library(LOLA)

# 


# load LOLA database
lolaPath = system.file("extdata", "hg19", package="LOLA")
regionSetDB = loadRegionDB(lolaPath)

GRList


# prepare DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))
simpleCache("combinedBRCAMethyl_noXY")
pca = prcomp(combinedBRCAMethyl_noXY[["methylProp"]], center = TRUE)
coordinates = combinedBRCAMethyl_noXY[["coordinates"]]
pcWeights = as.data.table(pca$rotation)


# run PC region set enrichment analysis
pcRegionSetEnrichment(loadingMat=pcWeights, coordinateDT = coordinates, 
                      GRList, 
                      PCsToAnnotate = c("PC1", "PC2"))



# confirm that results make sense
# PCs that are enriched for ATACseq from a certain cell type should be 
# able to separate that cell type from the others
# visualize


# could also do the analysis for hormone receptors in breast cancer:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4754852/
load("allBRCAexpression.RData")
# figuring out which patients are ER+ (inexact estimate based on gene expression data)
ER1Ind = grep(pattern = "ENSG00000091831", x = myExprDT[, "Gene"],
     ignore.case = TRUE) # ESR1
hist(as.numeric(myExprDT[ER1Ind, 2:ncol(myExprDT)]),
     breaks = seq(0,300,1))
sum(as.numeric(myExprDT[ER1Ind, 2:ncol(myExprDT)]) < 1) / (ncol(myExprDT)-1)
# progesterone receptoer positive
PGRInd = grep(pattern = "ENSG00000082175", x = myExprDT[, "Gene"],
              ignore.case = TRUE) # PGR
hist(as.numeric(myExprDT[PGRInd, 2:ncol(myExprDT)]),
     breaks = seq(0,400,1))
sum(as.numeric(myExprDT[PGRInd, 2:ncol(myExprDT)]) < 2) / (ncol(myExprDT)-1)
