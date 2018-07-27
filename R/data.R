#' A matrix with loadings
#' 
#' This object contains loadings
#' for PCA of DNA methylation data.
#' DNA methlyation data is Illumina 450k 
#' microarray data from breast cancer
#' patients from The Cancer Genome Atlas
#' (TCGA-BRCA, https://portal.gdc.cancer.gov/).
#' Each row corresponds to one cytosine and
#' the coordinates for these cytosines are
#' in the object brcaCoords1, (data("brcaCoords1")). 
#' Only cytosines on chr1 are included to keep
#' the example data small.
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaLoadings1
#' @usage data(brcaLoadings1)
#' @format A matrix object
NULL

#' A data.table, data.frame object with coordinates for cytosines
#' from chr1 included in the PCA. 
#' 
#' Corresponds to 
#' the rows of brcaLoadings1 and brcaMethylData1.
#' DNA methlyation data is Illumina 450k 
#' microarray data from breast cancer
#' patients from The Cancer Genome Atlas
#' (TCGA-BRCA, https://portal.gdc.cancer.gov/).
#' Only cytosines on chr1 are included to keep
#' the example data small.
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaCoord1
#' @usage data(brcaCoord1)
#' @format A data.table, data.frame object
NULL

#' A matrix with DNA methylation levels 
#' from chromosome 1 for four patients.
#' 
#' This object contains methylation levels (0 to 1)
#' for cytosines in chromosome 1 that were covered by
#' the DNA methylation microarray (Illumina 450k microarray).
#' Each row corresponds to one cytosine and
#' the coordinates for these cytosines are
#' in the object brcaCoords1, (data("brcaCoords1")). 
#' Only cytosines on chr1 are included to keep
#' the example data small. Columns are patients,
#' with TCGA patient identifiers as column names. The patients
#' were selected based on their PC1 score from
#' the PCA of DNA methylation on all chromosomes. The
#' patients with the two highest PC1 scores and 
#' the two lowest PC1 scores are included 
#' (see data("brcaPCScores") for the actual scores).
#' DNA methlyation data is Illumina 450k 
#' microarray data from breast cancer
#' patients from The Cancer Genome Atlas
#' (TCGA-BRCA, https://portal.gdc.cancer.gov/).
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaMethylData1
#' @usage data(brcaMethylData1)
#' @format A matrix object
NULL

#' A matrix with principal component scores for 
#' PCs 1-4 for four breast cancer patients.
#' 
#' This object contains PC scores for four patients
#' for PCs 1-4. Columns are PCs. Rows are patients,
#' with TCGA patient identifiers as row names. The patients
#' were selected based on their PC1 score from
#' the PCA of DNA methylation on all chromosomes. The
#' patients with the two highest PC1 scores and 
#' the two lowest PC1 scores are included.
#' DNA methlyation data is Illumina 450k 
#' microarray data from breast cancer
#' patients from The Cancer Genome Atlas
#' (TCGA-BRCA, https://portal.gdc.cancer.gov/).
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaPCScores
#' @usage data(brcaPCScores)
#' @format A matrix object
NULL

#' Estrogen receptor alpha binding regions.
#' 
#' Binding regions for estrogen receptor alpha (ESR1).
#' Only includes regions in chr1 to keep
#' the example data small.
#'
#' @docType data
#' @keywords datasets
#' @name esr1_chr1
#' @usage data(esr1_chr1)
#' @format A GRanges object
NULL

#' Gata3 binding regions.
#' 
#' Binding regions for gata3.
#' Only includes regions in chr1 to keep
#' the example data small.
#'
#' @docType data
#' @keywords datasets
#' @name gata3_chr1
#' @usage data(gata3_chr1)
#' @format A GRanges object
NULL

#' Nrf1 binding regions.
#' 
#' Binding regions for Nrf1.
#' Only includes regions in chr1 to keep
#' the example data small.
#'
#' @docType data
#' @keywords datasets
#' @name nrf1_chr1
#' @usage data(nrf1_chr1)
#' @format A GRanges object
NULL

#' Atf3 binding regions.
#' 
#' Binding regions for Atf3.
#' Only includes regions in chr1 to keep
#' the example data small.
#'
#' @docType data
#' @keywords datasets
#' @name atf3_chr1
#' @usage data(atf3_chr1)
#' @format A GRanges object
NULL

# # script for generating package data
# # restricting data to reduce how much memory the package takes up
# # first run part of 02-brca_PCRSA to get brcaMList and filteredMData
# 
# library(simpleCache)
# library(GenomicRanges)
# 
# setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))
# 
# # loading values
# simpleCache("allMPCA_657", assignToVariable = "mPCA")
# loadingMat <- mPCA$rotation
# 
# # coordinates
# coordinateDT <- brcaMList[["coordinates"]]
# 
# # region sets
# # some of the top ranked region sets for PC1 (in top 5)
# esr1 <- GRList[rsName == "Human_MCF-7_ESR1_E2-45min_Brown.bed"][[1]]
# gata3 <- GRList[rsName == "wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak"][[1]]
# # some of the bottom ranked region sets for PC1 (in bottom 1% of >2200 region sets)
# nrf1 <- GRList[rsName == "wgEncodeAwgTfbsSydhHepg2Nrf1IggrabUniPk.narrowPeak"][[1]]
# atf3 <- GRList[rsName == "wgEncodeAwgTfbsSydhK562Atf3UniPk.narrowPeak"][[1]]
# 
# 
# # restrict to chromosome 22 (reduce file sizes)
# chr22Ind <- coordinateDT$chr == "chr22"
# coord22 <- coordinateDT[chr22Ind,]
# # subset loadings, also restrict PCs
# load22 <- loadingMat[chr22Ind, paste0("PC", 1:4)]
# 
# grList <- GRangesList(er, gata, nrf1, atf3)
# ol22 <- lapply(X = grList, function(x) findOverlaps(query = x, subject = MIRA:::dtToGr(coord22)))
# 
# 
# # restrict to chromosome 1 (reduce file sizes)
# chr1Ind <- coordinateDT$chr == "chr1"
# brcaCoord1 <- coordinateDT[chr1Ind,]
# # subset loadings, also restrict PCs
# brcaLoadings1 <- loadingMat[chr1Ind, paste0("PC", 1:4)]
# 
# # restrict region sets to chr1
# grList <- GRangesList(er, gata3, nrf1, atf3)
# ol1 <- lapply(X = grList, function(x) findOverlaps(query = x, subject = MIRA:::dtToGr(coord1)))
# 
# # to restrict to chr1
# esr1_chr1 <- esr1[ seqnames(esr1) == "chr1"]
# gata3_chr1 <- gata3[ seqnames(gata3) == "chr1"]
# nrf1_chr1 <- nrf1[ seqnames(nrf1) == "chr1"]
# atf3_chr1 <- atf3[ seqnames(atf3) == "chr1"]
# 
# # save("", file = "coord1.RData") # reduce ~4x from in-memory size
# save("brcaLoadings1", file = "brcaLoadings1.RData", compress = "xz") # reduce ~7x
# save("brcaCoord1", file = "brcaCoord1.RData", compress = "xz")
# 
# save("esr1_chr1", file = "esr1_chr1.RData", compress = "xz")
# save("gata3_chr1", file = "gata3_chr1.RData", compress = "xz")
# save("nrf1_chr1", file = "nrf1_chr1.RData", compress = "xz")
# save("atf3_chr1", file = "atf3_chr1.RData", compress = "xz")
# 
# ### Get DNA methylation level data, filteredMData from another script
# filteredMData <- cbind(brcaMList[["coordinates"]], filteredMData)
# chr1MData <- filteredMData[filteredMData$chr == "chr1", ]
# 
# # selecting patients to include based on PC1, two with low PC scores, two with 
# # high PC scores
# pc1Order <- sort(mPCA$x[, "PC1"], decreasing = FALSE)
# lowScore <- names(pc1Order)[1:2]
# highScore <- names(pc1Order)[(length(pc1Order) - 1):length(pc1Order)]
# 
# brcaMethylData1 <- chr1MData[, c(lowScore, highScore), with=FALSE]
# brcaMethylData1 <- as.matrix(brcaMethylData1)
# save("brcaMethylData1", file = "brcaMethylData1.RData", compress = "xz")
# 
# ## also getting PC scores for PC 1-4 for these patients
# brcaPCScores <- mPCA$x[c(lowScore, highScore), paste0("PC", 1:4)]
# save("brcaPCScores", file="brcaPCScores.RData", compress="xz")
