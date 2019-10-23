# A matrix with loadings
# 
# This object contains loadings
# for PCA of DNA methylation data.
# DNA methlyation data is Illumina 450k 
# microarray data from breast cancer
# patients from The Cancer Genome Atlas
# (TCGA-BRCA, https://portal.gdc.cancer.gov/).
# Each row corresponds to one cytosine and
# the coordinates for these cytosines are
# in the object brcaMCoord1, (data("brcaMCoord1"),
# hg38 genome). 
# Only cytosines on chr1 are included to keep
# the example data small.
#
#
# @docType data
# @keywords datasets
# @name brcaLoadings1
# @usage data(brcaLoadings1)
# @format A matrix object


#' A data.frame object with coordinates for cytosines
#' from chr1 included in the PCA. 
#' 
#' Corresponds to 
#' the rows of brcaLoadings1 and brcaMethylData1.
#' DNA methlyation data is Illumina 450k 
#' microarray data from breast cancer
#' patients from The Cancer Genome Atlas
#' (TCGA-BRCA, https://portal.gdc.cancer.gov/).
#' Coordinates correspond to the hg38 genome version.
#' Only cytosines on chr1 are included to keep
#' the example data small.
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaMCoord1
#' @usage data(brcaMCoord1)
#' @format A data.frame object
NULL

#' A matrix with DNA methylation levels 
#' from chromosome 1 for four patients.
#' 
#' This object contains methylation levels (0 to 1)
#' for cytosines in chromosome 1 that were covered by
#' the DNA methylation microarray (Illumina 450k microarray).
#' Each row corresponds to one cytosine and
#' the coordinates for these cytosines are
#' in the object brcaMCoord1, (data("brcaMCoord1"), hg38 genome). 
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
#' (TCGA-BRCA, https://portal.gdc.cancer.gov/), 
#' hg38 genome.
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaPCScores
#' @usage data(brcaPCScores)
#' @format A matrix object
NULL

#' A data.frame with principal component scores for 
#' PCs 1-4 for 657 breast cancer patients as well
#' as a column with estrogen receptor status.
#' 
#' This object contains PC scores for 657 patients
#' for PCs 1-4. Columns are PCs as well
#' as a column with estrogen receptor status. Rows are patients,
#' with TCGA patient identifiers as row names. Patients
#' were selected from all BRCA patients in TCGA based on having complete
#' metadata information for estrogen receptor status 
#' and progesterone receptor status as well as 
#' having 450k microarray data. 
#' PCA was done on the Illumina 450k
#' DNA methlyation  
#' microarray data (TCGA-BRCA, https://portal.gdc.cancer.gov/),
#' hg38 genome.
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaPCScores657
#' @usage data(brcaPCScores657)
#' @format A data.frame object
NULL

#' Estrogen receptor alpha binding regions.
#' 
#' Binding regions for estrogen receptor alpha (ESR1).
#' hg 38 genome version. Only includes regions in chr1 to keep
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
#' hg38 genome version.
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
#' hg38 genome version.
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
#' hg38 genome version.
#' Only includes regions in chr1 to keep
#' the example data small.
#'
#' @docType data
#' @keywords datasets
#' @name atf3_chr1
#' @usage data(atf3_chr1)
#' @format A GRanges object
NULL

#' Example COCOA Results (made up)
#' 
#' A data.frame with example COCOA results.
#' 5 region sets with names given by rsScores$rsName.
#' Each region set has a score for each PC. Scores
#' for real region sets would normally be orders of 
#' magnitude smaller.
#'
#' @docType data
#' @keywords datasets
#' @name rsScores
#' @usage data(rsScores)
#' @format A data.frame object
NULL

#' A data.frame object with coordinates for BRCA ATAC-seq peak regions
#' from chr1. 
#' 
#' Corresponds to 
#' the rows of brcaATACData1.
#' The ATAC-seq data is from breast cancer
#' patients from The Cancer Genome Atlas
#' (TCGA-BRCA, Corces et. al, 2018, doi: 10.1126/science.aav1898,
#' https://atacseq.xenahubs.net/download/brca/brca_peak_Log2Counts_dedup).
#' Coordinates correspond to the hg38 genome version.
#' Only regions on chr1 are included to keep
#' the example data small.
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaATACCoord1
#' @usage data(brcaATACCoord1)
#' @format A data.frame object 
NULL

#' A matrix with ATAC-seq counts in peak regions 
#' from chromosome 1 for four patients.
#' 
#' Each row corresponds to one region and
#' the coordinates for these regions are
#' in the object brcaATACCoord1, (data("brcaATACCoord1"), hg38 genome). 
#' Only regions on chr1 are included to keep
#' the example data small. Columns are patients,
#' with TCGA patient identifiers as column names. 
# The patients with the two highest PC1 scores and
# the two lowest PC1 scores are included
# (see data("brcaPCScores") for the actual scores).
#' ATAC-seq data is from breast cancer
#' patients from The Cancer Genome Atlas
#' (TCGA-BRCA, Corces et. al, 2018, doi: 10.1126/science.aav1898,
#' https://atacseq.xenahubs.net/download/brca/brca_peak_Log2Counts_dedup).
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaATACData1
#' @usage data(brcaATACData1)
#' @format A matrix object
NULL

#' A data.frame with patient metadata for breast 
#' cancer patients.
#' 
#' Has metadata for patients for which DNA methylation
#' or chromatin accessibility data was included as package data.
#' Rows are patients,
#' with TCGA patient identifiers as row names and the column "subject_ID". 
#' Also includes columns: ER_status, ER_percent, age_diagnosis, days_to_death,
#' and days_to_last_follow_up.
#' Metadata is from The Cancer Genome Atlas
#' (TCGA-BRCA, https://portal.gdc.cancer.gov/).
#'
#'
#' @docType data
#' @keywords datasets
#' @name brcaMetadata
#' @usage data(brcaMetadata)
#' @format A data.frame object
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
# esr1 <- GRList[rsName == 
#                    "Human_MCF-7_ESR1_E2-45min_Brown.bed"][[1]]
# gata3 <- GRList[rsName == 
#                     "wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak"][[1]]
# # some of the bottom ranked region sets 
# # for PC1 (in bottom 1% of >2200 region sets)
# nrf1 <- GRList[rsName == 
#                    "wgEncodeAwgTfbsSydhHepg2Nrf1IggrabUniPk.narrowPeak"][[1]]
# atf3 <- GRList[rsName == 
#                    "wgEncodeAwgTfbsSydhK562Atf3UniPk.narrowPeak"][[1]]
# 
# 
# # restrict to chromosome 22 (reduce file sizes)
# chr22Ind <- coordinateDT$chr == "chr22"
# coord22 <- coordinateDT[chr22Ind,]
# # subset loadings, also restrict PCs
# load22 <- loadingMat[chr22Ind, paste0("PC", 1:4)]
# 
# grList <- GRangesList(er, gata, nrf1, atf3)
# ol22 <- lapply(X = grList, 
#                function(x) findOverlaps(query = x, 
#                                         subject = MIRA:::dtToGr(coord22)))
# 
# 
# # restrict to chromosome 1 (reduce file sizes)
# chr1Ind <- coordinateDT$chr == "chr1"
# brcaMCoord1 <- coordinateDT[chr1Ind,]
# brcaMCoord1 <- as.data.frame(brcaMCoord1)
# # subset loadings, also restrict PCs
# brcaLoadings1 <- loadingMat[chr1Ind, paste0("PC", 1:4)]
# 
# # restrict region sets to chr1
# grList <- GRangesList(er, gata3, nrf1, atf3)
# ol1 <- lapply(X = grList, 
#               function(x) findOverlaps(query = x, 
#                                        subject = MIRA:::dtToGr(coord1)))
# 
# # to restrict to chr1
# esr1_chr1 <- esr1[ seqnames(esr1) == "chr1"]
# gata3_chr1 <- gata3[ seqnames(gata3) == "chr1"]
# nrf1_chr1 <- nrf1[ seqnames(nrf1) == "chr1"]
# atf3_chr1 <- atf3[ seqnames(atf3) == "chr1"]

# # add sequence info (seqInfo)
# exRSNames <- c("esr1_chr1", "gata3_chr1", "nrf1_chr1", "atf3_chr1")
# 
# # for loop is to avoid having to retype expression for each region set
# for (i in seq_along(exRSNames)) {
#     
#     # restrict seqinfo to only chromosomes that are present (chr1)
#     assignString = paste0("seqlevels(", exRSNames[i], ") <- as.character(unique(seqnames(", exRSNames[i], ")))")
#     eval(parse(text=assignString))
#     
#     assignString = paste0("isCircular(", exRSNames[i], ") <- rep(FALSE, length(seqinfo(", exRSNames[i], ")))")
#     eval(parse(text=assignString))
#     
#     # assign reference genome
#     assignString = paste0("genome(", exRSNames[i], ") <- rep(\"hg38\", length(seqinfo(", exRSNames[i], ")))")
#     eval(parse(text=assignString))
# 
#     # chr1 has length of 248,956,422 based on below link
#     # https://www.ncbi.nlm.nih.gov/grc/human/data
#     assignString = paste0("seqlengths(", exRSNames[i], ") <- 248956422")
#     eval(parse(text=assignString))
# }

# 
# # save("", file = "coord1.RData") # reduce ~4x from in-memory size
# save("brcaLoadings1", 
#      file = "brcaLoadings1.RData", 
#      compress = "xz") # reduce ~7x
# save("brcaMCoord1", file = "brcaMCoord1.RData", compress = "xz")
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
# # selecting patients to include based on PC1, two with low PC scores, 
# # and two with high PC scores
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
##################
# # signalMat, signalCoord, patientMetadata
# source(paste0(Sys.getenv("CODE"), "COCOA_paper/src/00-init.R"))
# loadBRCADNAm(loadingMat = FALSE, pcScores = FALSE)
# # making second round of BRCA DNAm data
# 
# signalCoord <- COCOA:::dtToGr(signalCoord)
# data("brcaMethylData1")
# 
# data("esr1_chr1")
# data("gata3_chr1")
# data("atf3_chr1")
# data("nrf1_chr1")
# 
# allTfs <- c(esr1_chr1, gata3_chr1, atf3_chr1, nrf1_chr1)
# 
# # expand regions so meta region profile can be created
# allTfs <- resize(allTfs, width=12000, fix="center")
# 
# # regions that are covered by any TF
# allTfs <- reduce(allTfs)
# # length(unique(queryHits(findOverlaps(query = brcaATACCoord1, subject = allTfs))))
# overlapsRegionSets <- rep(FALSE, length(signalCoord))
# overlapsRegionSets[unique(queryHits(findOverlaps(query = signalCoord,
#                                                  subject = allTfs)))] <- TRUE
# 
# 
# brcaMCoord1 <- signalCoord[overlapsRegionSets]
# 
# # select samples, try to get a somewhat even distribution across PC's 1 and 4
# data("brcaPCScores657")
# 
# set.seed(1234)
# PC1Vec <- brcaPCScores657[, c("PC1", "PC4")]
# PC1Dist <- as.matrix(dist(PC1Vec))
# neighbor3Dist <- apply(X = PC1Dist, 2, function(x) mean(sort(x)[2:4]))
# 
# sampleInd <- sample(1:nrow(brcaPCScores657), size = 300, replace = FALSE, prob = neighbor3Dist)
# plot(brcaPCScores657[sampleInd, c("PC1", "PC4")])
# sampleNames <- row.names(brcaPCScores657)[sampleInd]
# sampleNamesDNAM <- sampleNames
# 
# brcaMethylData1 <- signalMat[overlapsRegionSets, sampleNames]
# save(brcaMethylData1, file="brcaMethylData1.RData", compress="xz")
# save(brcaMCoord1, file="brcaMCoord1.RData", compress="xz")

############ data for rsRankingIndex
# setwd("/data")
# PC1 <- c(1, 2, 3, 5, 7)
# PC2 <- c(5, 1, 3, 2, 4)
# rsName <- c("rs1", "rs2", "rs3", "rs4", "rs5")
# rsScores <- data.frame(PC1, PC2, rsName)
# save(rsScores, file = "rsScores.RData", compress = "xz")
# rsRankingIndex(rsScores = rsScores, PCsToAnnotate = c("PC1", "PC2"))

############ brcaPCScores657
# # # first load patient metadata and PCA for the 657 patients
# # brcaPCScores657 <- as.data.frame(allMPCA$x[, c("PC1", "PC2", "PC3", "PC4")])
# # patientMetadata <- as.data.frame(patientMetadata)
# # row.names(patientMetadata) <- patientMetadata$subject_ID
# # brcaPCScores657 <- cbind(brcaPCScores657, ER_Status = patientMetadata[row.names(brcaPCScores657), "ER_status"])
# # save(brcaPCScores657, file="brcaPCScores657.RData", compress = "xz")
# 
########### BRCA ATAC-seq data
# # first load BRCA ATAC-seq data from Corces et al. (signalMat, signalCoord)
# loadBRCAatac(signalMat = TRUE, signalCoord = TRUE,
#              loadingMat = TRUE, pcScores = TRUE)
# 
# # subset to only ATAC regions that are on chr1 and also
# # overlap with the package's existing region sets
# chr1Ind <- as.logical(seqnames(signalCoord) == "chr1")
# 
# data("esr1_chr1")
# data("gata3_chr1")
# data("atf3_chr1")
# data("nrf1_chr1")
# allTfs <- c(esr1_chr1, gata3_chr1, atf3_chr1, nrf1_chr1)
# # expand regions so meta region profile can be created
# allTfs <- resize(allTfs, width=12000, fix="center")
# allTfs <- reduce(allTfs)
# # length(unique(queryHits(findOverlaps(query = brcaATACCoord1, subject = allTfs))))
# overlapsRegionSets <- rep(FALSE, length(signalCoord))
# overlapsRegionSets[unique(queryHits(findOverlaps(query = signalCoord,
#                                                  subject = allTfs)))] <- TRUE
# keepInd <- chr1Ind & overlapsRegionSets
# 
# brcaATACCoord1 <- signalCoord[keepInd]
# brcaATACData1 <- signalMat[keepInd, ]
# 
# save("brcaATACCoord1", file = "brcaATACCoord1.RData", compress = "xz")
# 
# sampleNames <- colnames(brcaATACData1)
# # pca of all the data
# a = prcomp(t(signalMat))
# plot(a$x[, 1:2])
# patientMetadata <- as.data.frame(patientMetadata)
# row.names(patientMetadata) <- patientMetadata$subject_ID
# annoScores <- cbind(data.frame(a$x),
#                     ER_status=as.factor(patientMetadata[sampleNames, "ER_status"]))
# ggplot(data = annoScores, mapping = aes(x=PC1, y=PC2)) + geom_point(aes(color=ER_status))
# sampleInd = ((annoScores$PC1 > 12) | (annoScores$PC1 < -15)) & !is.na(annoScores$ER_status)
# 
# # take out samples without ER status
# brcaATACData1 <- brcaATACData1[, !is.na(annoScores$ER_status)]
# sampleNames <- colnames(brcaATACData1)
# 
# # trying to get a somewhat even distribution across PC1
# set.seed(1234)
# PC1Vec <- a$x[sampleNames, "PC1"]
# PC1Dist <- as.matrix(dist(PC1Vec))
# neighbor3Dist <- apply(X = PC1Dist, 2, function(x) mean(sort(x)[2:4]))
# sampleInd <- sample(1:ncol(brcaATACData1), size = 37, replace = FALSE, prob = neighbor3Dist)
# 
# sampleNames <- sampleNames[sampleInd]
# sampleNamesATAC <- sampleNames
# a = prcomp(t(brcaATACData1[, sampleInd]))
# plot(a$x[, 1:2])
# patientMetadata <- as.data.frame(patientMetadata)
# row.names(patientMetadata) <- patientMetadata$subject_ID
# annoScores <- cbind(data.frame(a$x),
#                     ER_status=as.factor(patientMetadata[sampleNames, "ER_status"]))
# ggplot(data = annoScores, mapping = aes(x=PC1, y=PC2)) + geom_point(aes(color=ER_status))
# 
# brcaATACData1 <- brcaATACData1[, sampleInd]
# 
# # only keep 3 sig figs
# brcaATACData1 <- round(brcaATACData1, 2)
# 
# save("brcaATACData1", file = "brcaATACData1.RData", compress = "xz")
# 
# # we no longer need to include the loadings or PC scores because we can
# # do PCA to generate them in the vignette
# # brcaATACLoadings1 <- loadingMat[keepInd, paste0("PC", 1:4)]
# # save("brcaATACLoadings1", file = "brcaATACLoadings1.RData", compress = "xz")
# # brcaATACPCScores <- pcScores[c(), paste0("PC", 1:4)]
# # save("brcaATACPCScores", file = "brcaATACPCScores", compress = "xz")
# # 
# ###########################################
# # get brcaMetadata
# sampleNames <- unique(c(sampleNamesDNAM, sampleNamesATAC))
# patientMetadata <- as.data.frame(patientMetadata)
# row.names(patientMetadata) <- patientMetadata$subject_ID
# brcaMetadata <- patientMetadata[sampleNames, c("subject_ID", "ER_status", "ER_percent",
#                                                "age_diagnosis", "days_to_death",
#                                                "days_to_last_follow_up")]
# brcaMetadata$ER_status = as.factor(brcaMetadata$ER_status)
# save(brcaMetadata, file = "brcaMetadata.RData", compress="xz")
