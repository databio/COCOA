# PCA Region Set Enrichment Analysis of Ewing Sarcoma

library(MIRA)
library(simpleCache)
library(data.table)
library(LOLA)
library(GenomicRanges) # GRangesList, resize
source(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/PRA.R"))


# reading the files in and creating annotation data.table
setwd(paste0(Sys.getenv("PROCESSED"), "ews_patients/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "ews_patients/analysis/plots/"))
setCacheDir(paste0(Sys.getenv("PROCESSED"), "ews_patients/RCache/"))
ewsFiles = list.files(pattern = "RRBS_cpgMethylation_EWS.+bed", recursive = TRUE)
mDTList = lapply(ewsFiles, BSreadBiSeq)
mDTList = addMethPropCol(mDTList)

########## preprocessing and finding out Cs that all samples share ############
# creating object with x and y chromosomes excluded because this is something 
# that I will want to have for several analyses
simpleCache("mDTList_noXY", {
    excludeXYM = function(x) {
        x[!(chr=="chrX" | chr=="chrY" | chr=="chrM"), ]
    }
    mDTList_noXY = lapply(mDTList, function(x) excludeXYM(x))
})

# add back on the chrBase column, which will be used to find shared cytosines between samples
mDTList_noXY = lapply(mDTList_noXY, function(x) x[, chrBase := paste0(chr, ".", start)])
# mDTList_noXY = lapply(mDTList_noXY, function(x) x[, chrBase := paste0(chr, ".", start)])

# restrict to C's that are in all samples
# might bias towards certain types of regions (cpg islands/promoters/genes)
# perhaps include less regulatory elements that were only covered in some samples
simpleCache("sharedCVec", {
    # I would first get vector of all cpgs (found in at least one sample)
    # however if they are not in the first sample, they are not in all samples 
    # so just take cpgs in first sample
    sharedCVec = mDTList_noXY[[1]]$chrBase
    
    # progressively elimate cpgs that are not in other samples
    for (i in 2:length(mDTList_noXY)) {
        
        # only keeps shared cytosines each time
        sharedCVec = sharedCVec[sharedCVec %in% mDTList_noXY[[i]]$chrBase]
        
    }
    
    # cytosine coordinate was not always in order of increasing coordinate
    sharedCVec = sort(sharedCVec, decreasing = FALSE)
    
})

# then subset each sample by these shared cytosines
consensusMethyl = function(x, y) {
    # selecting consensus cytosines then ordering in ascending order
    # according to chrBase column
    data.table::setorder(x[x$chrBase %in% y, ], chrBase)
}
simpleCache("sharedCDTList", {
    sharedCDTList = lapply(mDTList_noXY, consensusMethyl, y=sharedCVec)
    lapply(sharedCDTList, function(x) x[, chrBase := NULL])
    sharedCDTList
})

# separating methylProp and coverage into separate objects 
sharedCMethylProp = lapply(sharedCDTList, function(x) x$methylProp)
sharedCCoverage = lapply(sharedCDTList, function(x) x$coverage)

# merge all samples into one data.table
# make column names the sample IDs
simpleCache("bigSharedC", {
    # rows already in the same order
    bigSharedC = list()
    bigSharedC[[1]] = sharedCDTList[[1]][, c("chr", "start")]
    bigSharedC[[2]] = do.call(cbind, sharedCMethylProp)
    bigSharedC[[3]] = do.call(cbind, sharedCCoverage)
    setattr(bigSharedC, "names", c("coordinates", "methylProp", "coverage"))
    bigSharedC
})

########################################

# reading in the region sets
# load LOLA database
lolaPath = paste0(Sys.getenv("REGIONS"), "LOLACore/hg38/")
#lolaPath = system.file("extdata", "hg19", package="LOLA")
regionSetDB = loadRegionDB(lolaPath)
loRegionAnno = regionSetDB$regionAnno
# a549Ind = grep("a549", loRegionAnno$cellType, ignore.case = TRUE)
sheff_dnaseInd = grep("sheffield_dnase", loRegionAnno$collection, ignore.case = TRUE)
# mcf7Ind = grep("mcf-7", loRegionAnno$cellType, ignore.case = TRUE)
# k562Ind = grep("k562", loRegionAnno$cellType,  ignore.case = TRUE)
# GRList = GRangesList(regionSetDB$regionGRL[c(a549Ind, mcf7Ind)])
GRList = GRangesList(regionSetDB$regionGRL[-sheff_dnaseInd])
# adding ER Chipseq dataset
erSet = fread(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/inst/extdata/",
                     "GSM2305313_MCF7_E2_peaks_hg38.bed"))
setnames(erSet, c("V1", "V2", "V3"), c("chr", "start", "end"))
GRList = c(GRangesList(MIRA:::dtToGr(erSet)), GRList)


#################################################################
# doing PCA of the methylation data
mData = bigSharedC$methylProp

simpleCache("allMPCA", {
    prcomp(t(mData), center = TRUE)
})
allMPCAWeights = as.data.table(allMPCA$rotation)
mIQR = apply(mData, 1, IQR)
simpleCache("top10MPCA", {
    prcomp(t(mData[mIQR >= quantile(mIQR, 0.9), ]), center = TRUE)
})
coordinates = bigSharedC[["coordinates"]]
top10Coord = coordinates[mIQR >= quantile(mIQR, 0.9), ]
top10PCWeights = as.data.table(top10MPCA$rotation)
# top10TSNE = Rtsne(X = top10MPCA$x[, 1:50], pca = FALSE, max_iter=5000,
#                   perplexity = 30)
# plot(top10TSNE$Y)
# run PC region set enrichment analysis
simpleCache("rsEnrichment", {
    rsEnrichment = pcRegionSetEnrichment(loadingMat=allMPCAWeights, coordinateDT = coordinates, 
                                         GRList, 
                                         PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"), permute=FALSE)
    # rsNames = c("Estrogen_Receptor", loRegionAnno$filename[c(a549Ind, mcf7Ind)])
    rsNames = c("Estrogen_Receptor", loRegionAnno$filename[-sheff_dnaseInd])
    rsEnrichment[, rsNames:= rsNames]
    rsEnrichment
})

View(rsEnrichment[order(PC1,decreasing = TRUE)])

# check whether is enrichment is specific to this region set by
# seeing if loading values have a spike in the center of these region sets
# compared to surrounding genome 
GRList = lapply(GRList, resize, width = 10000, fix="center")
simpleCache("pcProfAllM", {
    pcProfAllM = pcEnrichmentProfile(loadingMat = allMPCAWeights, coordinateDT = coordinates,
                                     GRList=GRList, PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                                     binNum = 21)
    
    length(rsNames) == length(pcProfAllM)
    rsNames = c("Estrogen_Receptor", loRegionAnno$filename[-sheff_dnaseInd])
    names(pcProfAllM) <- rsNames 
    pcProfAllM
})


pcP = pcProfAllM
# plot(pcP$PC1, type="l")
# plot(pcP$PC2, type="l")
# plot(pcP$PC3, type="l")
# plot(pcP$PC4, type="l")
# plot(pcP$PC5, type="l")


grDevices::pdf(paste0(Sys.getenv("PLOTS"), "allMPCProfiles.pdf"))
for (i in 1:length(pcProfAllM)) {
    plot(pcP[[i]]$PC1, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC2, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC3, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC4, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC5, type="l") + title(rsNames[i])
}
dev.off()




# running the enrichment analysis



