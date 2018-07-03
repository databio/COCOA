# script for generating package data
# restricting data to reduce how much memory the package takes up

library(simpleCache)
library(GenomicRanges)

setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))

# loading values
simpleCache("allMPCA_657", assignToVariable = "mPCA")
loadingMat = mPCA$rotation

# coordinates
coordinateDT = brcaMList[["coordinates"]]

# region sets
# some of the top ranked region sets for PC1 (in top 5)
er = GRList[rsName == "Human_MCF-7_ESR1_E2-45min_Brown.bed"][[1]]
gata = GRList[rsName == "wgEncodeAwgTfbsSydhMcf7Gata3UcdUniPk.narrowPeak"][[1]]
# some of the bottom ranked region sets for PC1 (in bottom 1% of >2200 region sets)
nrf1 = GRList[rsName == "wgEncodeAwgTfbsSydhHepg2Nrf1IggrabUniPk.narrowPeak"][[1]]
atf3 = GRList[rsName == "wgEncodeAwgTfbsSydhK562Atf3UniPk.narrowPeak"][[1]]


# restrict to chromosome 22 (reduce file sizes)
chr22Ind = coordinateDT$chr == "chr22"
coord22 = coordinateDT[chr22Ind,]
# subset loadings, also restrict PCs
load22 = loadingMat[chr22Ind, paste0("PC", 1:4)]

grList = GRangesList(er, gata, nrf1, atf3)
ol22 = lapply(X = grList, function(x) findOverlaps(query = x, subject = MIRA:::dtToGr(coord22)))


# restrict to chromosome 1 (reduce file sizes)
chr1Ind = coordinateDT$chr == "chr1"
coord1 = coordinateDT[chr1Ind,]
# subset loadings, also restrict PCs
load1 = loadingMat[chr1Ind, paste0("PC", 1:4)]


grList = GRangesList(er, gata, nrf1, atf3)
ol1 = lapply(X = grList, function(x) findOverlaps(query = x, subject = MIRA:::dtToGr(coord1)))

save("coord1", file = "coord1.RData") # reduce ~4x from in-memory size
save("coord1", file = "coord1.RData", compress = "xz") # reduce ~7x
