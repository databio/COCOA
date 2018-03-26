# PCA Region Set Enrichment Analysis of Ewing Sarcoma

library(MIRA)
library(simpleCache)
library(data.table)

# reading the files in and creating annotation data.table
setwd(paste0(Sys.getenv("PROCESSED"), "ews_patients/"))
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




# doing PCA of the methylation data




# running the enrichment analysis



