#' A matrix with loadings
#' 
#' This object contains loadings
#' for PCA of DNA methylation data.
#' DNA methlyation data is from breast cancer
#' patients from The Cancer Genome Atlas
#' (TCGA-BRCA, https://portal.gdc.cancer.gov/).
#' Each row corresponds to one cytosine and
#' the coordinates for these cytosines are
#' in the object brcaCoords1. Only
#' cytosines on chr1 are included to keep
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
#' the rows of brcaLoadings1.
#' DNA methlyation data is from breast cancer
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
# 
# library(simpleCache)
# library(GenomicRanges)
# 
# setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))
# 
# # loading values
# simpleCache("allMPCA_657", assignToVariable = "mPCA")
# loadingMat = mPCA$rotation
# 
# # coordinates
# coordinateDT = brcaMList[["coordinates"]]
# 
# # region sets
# # some of the top ranked region sets for PC1 (in top 5)
# esr1 = GRList[rsName == "Human_MCF-7_ESR1_E2-45min_Brown.bed"][[1]]
# gata3 = GRList[rsName == "wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak"][[1]]
# # some of the bottom ranked region sets for PC1 (in bottom 1% of >2200 region sets)
# nrf1 = GRList[rsName == "wgEncodeAwgTfbsSydhHepg2Nrf1IggrabUniPk.narrowPeak"][[1]]
# atf3 = GRList[rsName == "wgEncodeAwgTfbsSydhK562Atf3UniPk.narrowPeak"][[1]]
# 
# 
# # restrict to chromosome 22 (reduce file sizes)
# chr22Ind = coordinateDT$chr == "chr22"
# coord22 = coordinateDT[chr22Ind,]
# # subset loadings, also restrict PCs
# load22 = loadingMat[chr22Ind, paste0("PC", 1:4)]
# 
# grList = GRangesList(er, gata, nrf1, atf3)
# ol22 = lapply(X = grList, function(x) findOverlaps(query = x, subject = MIRA:::dtToGr(coord22)))
# 
# 
# # restrict to chromosome 1 (reduce file sizes)
# chr1Ind = coordinateDT$chr == "chr1"
# brcaCoord1 = coordinateDT[chr1Ind,]
# # subset loadings, also restrict PCs
# brcaLoadings1 = loadingMat[chr1Ind, paste0("PC", 1:4)]
# 
# # restrict region sets to chr1
# grList = GRangesList(er, gata3, nrf1, atf3)
# ol1 = lapply(X = grList, function(x) findOverlaps(query = x, subject = MIRA:::dtToGr(coord1)))
# 
# # to restrict to chr1
# esr1_chr1 = esr1[ seqnames(esr1) == "chr1"]
# gata3_chr1 = gata3[ seqnames(gata3) == "chr1"]
# nrf1_chr1 = nrf1[ seqnames(nrf1) == "chr1"]
# atf3_chr1 = atf3[ seqnames(atf3) == "chr1"]
# 
# # save("", file = "coord1.RData") # reduce ~4x from in-memory size
# save("brcaLoadings1", file = "brcaLoadings1.RData", compress = "xz") # reduce ~7x
# save("brcaCoord1", file = "brcaCoord1.RData", compress = "xz")
# 
# save("esr1_chr1", file = "esr1_chr1.RData", compress = "xz")
# save("gata3_chr1", file = "gata3_chr1.RData", compress = "xz")
# save("nrf1_chr1", file = "nrf1_chr1.RData", compress = "xz")
# save("atf3_chr1", file = "atf3_chr1.RData", compress = "xz")


