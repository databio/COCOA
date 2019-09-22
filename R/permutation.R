# functions that have to do with the permutation test
# and getting p values from that test with the gamma distribution


# permutation test by shuffling sample labels
# expected inputs: 
#' permutation test by shuffling sample labels
#' 
#' For reproducibility, set seed with 'set.seed()' function before running.
#' @param nPerm numeric. The number of permutations to do.
#' @param sampleLabels data.frame/matrix. Sample labels/values that 
#' you are running COCOA to find region sets associated with. These 
#' values will be shuffled for the permutation test. Rows are samples.
#' Each column is a sample label.
#' @param genomicSignal
#' @param signalCoord
#' @param GRList
#' @param colsToAnnotate character. The column names of `sampleLabels` that
#' you want to test.
#' @param dataID character. A unique identifier for this dataset 
#' (for saving results)
#' @param variationMetric character. Either "cor" (correlation), "pcor" (partial
#' correlation), or "cov" (covariation)
#' @param useSimpleCache logical. 
#' @param cacheDir character.
#' @param correctionMethod character. P value correction method. For acceptable 
#' arguments see ?p.adjust() (method parameter) 
#' @param resultType character. "pval" or "zscore"
#' @param ... character. Optional additional arguments for simpleCache
#' 
# # for visualization
#' @param realRSScores
#' @return Returns a list where each item is a data.frame of COCOA results 
#' from a separate permutation
#' 
#' @export

runCOCOAPerm <- function(genomicSignal,
                                 signalCoord,
                                 GRList,
                                 signalCol=c("PC1", "PC2"),
                                 signalCoordType = "default",
                                 scoringMetric="default",
                                 variationMetric="cor",
                                 nPerm=300,
                                 useSimpleCache=TRUE,
                                 cacheDir=getwd(),
                                 dataID="tmp", ...) {


    indList <- list()
    # generate random indices for shuffling of samples
    for (i in 1:nPerm) {
        indList[[i]] <- sample(1:nrow(sampleLabels), nrow(sampleLabels))
    }
    
    # # plug random indices into parallelized function
    # rsPermScores <- COCOA:::lapplyAlias(X= indList, FUN=function(x) corPerm(randomInd=x,
    #                                                         genomicSignal=methylMat,
    #                                                         signalCoord=signalCoord,
    #                                                         GRList=GRList,
    #                                                         calcCols=colsToAnnotate,
    #                                                         sampleLabels=latentFactors))
    
    # create the main permutation cache
    simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), {
        
        rsPermScores <- list()
        for (i in seq_along(indList)) {
            # for (i in (length(rsPermScores) + 1):nPerm) {
            onePermCacheName <- paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID, "_Cache", i)
            # create sub caches, one for each permutation
            simpleCache(onePermCacheName, cacheSubDir = paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), {
                
                tmp <- corPerm(randomInd=indList[[i]],
                               genomicSignal=genomicSignal,
                               signalCoord=signalCoord,
                               GRList=GRList,
                               calcCols=colsToAnnotate,
                               sampleLabels=sampleLabels,
                               variationMetric = variationMetric)
                message(i) # must be ahead of object that is saved as cache, not after
                tmp
                
            })
        }
        
        # combining all individual permutations/caches into one object
        for (i in seq_along(indList)) {
            onePermCacheName <- paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID, "_Cache", i)
            rsPermScores[[i]] <- get(onePermCacheName)
        }
        rsPermScores
        
    }, assignToVariable="rsPermScores")
    
    .analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
    .plotSubdir = paste0(plotSubdir, "StatsPlots", .analysisID, "/")
    if (!dir.exists(ffPlot(.plotSubdir))) {
        dir.create(ffPlot(.plotSubdir))
    }
    
    # remove region sets that had no overlap
    keepInd = apply(rsPermScores[[1]], MARGIN = 1, FUN = function(x) !any(is.na(x)))
    
    nullDistList = lapply(X = 1:nrow(rsPermScores[[1]]),
                          FUN = function(x) extractNullDist(resultsList=rsPermScores, rsInd = x))
    
    # screen out region sets with no overlap
    nullDistList = nullDistList[keepInd]
    realRSScores = realRSScores[keepInd, ]
    
    simpleCache(paste0("permPValsUncorrected", .analysisID), {
        rsPVals = getPermStat(rsScores=realRSScores, nullDistList=nullDistList,
                              calcCols=colsToAnnotate, whichMetric = "pval")
        rsPVals
    }, recreate = TRUE, reload = TRUE)
    
    simpleCache(paste0("permZScores", .analysisID), {
        rsZScores = getPermStat(rsScores=realRSScores, nullDistList=nullDistList,
                                calcCols=colsToAnnotate, whichMetric = "zscore")
        rsZScores
        
    }, recreate = TRUE, reload=TRUE)
    
    
    #topRSInd = rsRankingIndex(rsScores = rsZScores, signalCol = colsToAnnotate)
    
    #################
    # simpleCache(paste0("rsScore_", dataID), assignToVariable = "realRSScores")
    # p-values based on fitted gamma distributions
    correctionMethod = "BH" # input to p.adjust
    gPValDF = getGammaPVal(scores = realRSScores[, colsToAnnotate, drop=FALSE], nullDistList = nullDistList, method = "mme", realScoreInDist = TRUE)
    gPValDF = apply(X = gPValDF, MARGIN = 2, FUN = function(x) p.adjust(p = x, method = correctionMethod))
    gPValDF = cbind(gPValDF, realRSScores[, colnames(realRSScores)[!(colnames(realRSScores) %in% colsToAnnotate)]])
    
    simpleCache(paste0("permPValsCorrected", .analysisID), {
        gPValDF
    }, recreate = TRUE, reload = TRUE)
    
}

# @param genomicSignal columns of dataMat should be samples/patients, rows should be genomic signal
# (each row corresponds to one genomic coordinate/range)
# @param sampleLabels Matrix or data.frame. Rows should be samples, 
# columns should be "features" 
# (whatever you want to get correlation with: eg PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns)
# @param calcCols character. the columns in `sampleLabels` for which to calculate
# correlation and then to run COCOA on
# @return data.frame. The output of runCOCOA for one permutation
corPerm <- function(randomInd, genomicSignal, 
                    signalCoord, GRList, calcCols,
                    sampleLabels, variationMetric = "cor") {
    
    # if vector is given, return error
    if (is.null(dim(sampleLabels))) {
        stop("`sampleLabels` should be a matrix or data.frame")
    }
    
    if (any(!(calcCols %in% colnames(sampleLabels)))) {
        stop("Not all specified columns are present in `sampleLabels`")
    }
    
    # subset to only calcCols
    sampleLabels = sampleLabels[, calcCols, drop=FALSE]
    
    # because names are dropped for a single column data.frame when indexing
    # single col data.frame is automatically converted to numeric
    featureNames = colnames(sampleLabels)
    # reorder the sample labels
    sampleLabels = data.frame(sampleLabels[randomInd, ])
    colnames(sampleLabels) = featureNames
    
    # calculate correlation
    featureLabelCor = createCorFeatureMat(dataMat = genomicSignal, 
                                          featureMat = sampleLabels, 
                                          centerDataMat = TRUE, 
                                          centerFeatureMat = TRUE,
                                          testType = variationMetric)
    
    # run COCOA
    thisPermRes = runCOCOA(signal=featureLabelCor, 
                           signalCoord=signalCoord, GRList=GRList, 
                           signalCol = calcCols, 
                           scoringMetric = "default", verbose = TRUE)
    
    # return
    return(thisPermRes)
    
}
# This function will take a list of results of permutation tests that included
# many region sets and return a data.frame/data.table with the null
# distribution for a single region set (row)
# @param resultsList each item in the list is a data.frame, one item for
# each permutation with the results of that permutation. Each row in the 
# data.frame is a region set. Rows in all the data.frames should be
# in the same order.
# @param rsInd numeric. The row number for the region set of interest.
extractNullDist <- function(resultsList, rsInd) {
    rowList = lapply(resultsList, FUN = function(x) x[rsInd, ])
    rsNullDist = as.data.frame(rbindlist(rowList))
    return(rsNullDist)
}

################################################################################
# p value functions 



#' Get p value after fitting a gamma distribution to the null distribution
#' @param nullDistList list of a data.frame. Each list item 
#' has null distributions for a single 
#' region set. Each column corresponds to a null distribution for that 
#' region set for a given variable/sample attribute.   
#' @param scores a data.frame. Has same columns as nullDistDF. One row per
#' region set (should be in same order as nullDistDF) The scores
#' that will be used to get p values.
#' @param realScoreInDist logical. Should the actual score (from 
#' test with no permutations) be included in the null distribution 
#' when fitting the gamma distribution
#' 
#' @return Returns a data.frame with p values, one column for each col in
#' scores and nullDistDF 

getGammaPVal <- function(scores, nullDistList, method="mme", realScoreInDist=FALSE, force=FALSE) {
    
    # make sure the same columns are present/in the same order
    
    
    colsToAnnotate = colnames(scores)
    
    if (realScoreInDist) {
        # to get a more accurate gamma distribution, include the score from unpermuted test.
        # add to each null distribution
        for (i in 1:nrow(scores)) {
            nullDistList[[i]] = rbind(nullDistList[[i]], as.numeric(scores[i, ]))
        }
    }
    
    
    # returns list, each item in list is also a list.
    # in sub list: each col in nullDistDF has one list item
    fittedDistList = lapply(X = nullDistList, function(x) fitGammaNullDistr(nullDistDF = x[, colsToAnnotate, drop=FALSE], 
                                                                            method=method, 
                                                                            force=force))
    
    pValList = list()
    # once for each region set
    for (i in seq_along(nullDistList)) {
        
        pValList[[i]] = as.data.frame(t(pGammaList(scoreVec = as.numeric(scores[i, ]), 
                                                   fitDistrList = fittedDistList[[i]])))
    }
    pValDF = as.data.frame(rbindlist(pValList))
    
    
    return(pValDF)
}


#' @param nullDistDF a data.frame. Has null distributions for a single 
#' region set. Each column corresponds to a null distribution for that 
#' region set for a given variable/sample attribute.   
#' 
#' @return Returns a list. Each list item is a "fitdist" object which is 
#' a fitted function 
#' for one of the columns in nullDistDF (output of fitdist() from
#' fitdistrplus.
#' These are gamma distributions and can be used to get p values for the 
#' null distribution so that a large number of permutations 
#' are not required.The list is in order of the columns and will
#' have the names of the data.frame columns. 

fitGammaNullDistr <- function(nullDistDF, method="mme", force=FALSE) {
    if (force) {
        modelList = apply(X = nullDistDF, 
                          MARGIN = 2, 
                          FUN = function(x) fitdistrplus::fitdist(data=x, 
                                                                  distr = "gamma", 
                                                                  method=method))
    } else {
        # try "method". if it fails, do method = "mme"
        modelList = apply(X = nullDistDF, 
                          MARGIN = 2, 
                          FUN = function(x) tryCatch({fitdistrplus::fitdist(data=x, 
                                                                            distr = "gamma", 
                                                                            method=method)}, 
                                                     error = function(e) {fitdistrplus::fitdist(data=x, 
                                                                                                distr = "gamma", 
                                                                                                method="mme")}))   
    }
    
    return(modelList)    
}

pGammaList <- function(scoreVec, fitDistrList) {
    pValVec = mapply(FUN = function(x, y) pgamma(q = y, 
                                                 shape = x$estimate["shape"], 
                                                 rate = x$estimate["rate"], 
                                                 lower.tail = FALSE), x = fitDistrList, y = scoreVec)
    return(pValVec)
}