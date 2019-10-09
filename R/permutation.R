# functions that have to do with the permutation test
# and getting p values from that test with the gamma distribution





#' Run COCOA permutations to get p-values
#' 
#' This is a convenience function that runs multiple steps of the 
#' permutation process together: it runs COCOA permutations, converts these
#' to null distributions, gets the empirical p value (which is limited by the
#' number of permutations), gets z scores, and fits a gamma distribution 
#' to each null distribution to estimate p values (not limited by the 
#' number of permutations),
#' Requires that the user has previously calculated the real COCOA scores. 
#' See these individual functions for more info on each step: corPerm, 
#' convertToFromNullDist, getPermStat, and getGammaPVal. 
#' 
#' 
#' For reproducibility, set seed with 'set.seed()' function before running.
#' @param nPerm numeric. The number of permutations to do.
#' @template genomicSignal
#' @template signalCoord
#' @template GRList
#' @template rsScores
#' @template sampleLabels
#' @param signalCol character. The column names of `sampleLabels` that
#' you want to test. These must also be columns in rsScores.
#' @template scoringMetric
#' @template absVal
#' @param dataID character. A unique identifier for this dataset 
#' (for saving results with simpleCache)
#' @param variationMetric character. Either "cor" (Pearson correlation), 
#' "pcor" (partial correlation), "spearmanCor (Spearman correlation) 
#' or "cov" (covariation). 
#' @param useSimpleCache logical. Whether to use save caches. Caches
#' will be created for each permutation so that if the function is disrupted
#' it can restart where it left off. The final results are also saved 
#' as a cache.  
#' @param cacheDir character.
#' @param correctionMethod character. P value correction method. Default
#' is "BH" for Benjamini and Hochberg false discovery rate. For acceptable 
#' arguments and more info see ?stats::p.adjust() (method parameter) 
#' @param gammaFitMethod character. method to use for fitting the gamma
#' distribution to null distribution. Options are 
#' "mme" (moment matching estimation), "mle" (maximum likelihood estimation), 
#' "qme" (quantile matching estimation), and "mge" (maximum goodness-of-fit 
#' estimation). See ?COCOA::getGammaPVal and 
#' ?fitdistrplus::fitdist() for more info.
#' @param realScoreInDist logical. Should the actual score (from 
#' test with no permutations) be included in the null distribution 
#' when fitting the gamma distribution. realScoreInDist=TRUE is 
#' recommended.
#' @template verbose
#' @param ... character. Optional additional arguments for simpleCache
#'
#' 
#' @return Returns a list where each item is a data.frame of COCOA results 
#' from a separate permutation
#' @examples data("brcaMCoord1")
#' data("brcaLoadings1")
#' data("esr1_chr1")
#' data("nrf1_chr1")
#' data("brcaMethylData1")
#' data("brcaPCScores657")
#' pcCor = corFeature
#' sampleLabels <- brcaPCScores657[colnames(brcaMethylData1), ]
#' sampleLabels$ER_Status <- scale(as.numeric(sampleLabels$ER_Status), 
#'                                center=TRUE, scale=FALSE)
#' # give the actual order of samples to randomInd to get the real scores
#' realRSScores <- corPerm(randomInd=1:4, genomicSignal=brcaMethylData1, 
#'         signalCoord=brcaMCoord1, GRList=GRangesList(esr1_chr1, nrf1_chr1),
#'         calcCols=c("PC1", "PC2"), sampleLabels=sampleLabels, 
#'         variationMetric="cor")
#' 
#' a=runCOCOAPerm(genomicSignal=brcaMethylData1, 
#'         signalCoord=brcaMCoord1, GRList=GRangesList(esr1_chr1, nrf1_chr1),
#'         rsScores=realRSScores, 
#'         sampleLabels=sampleLabels, signalCol=c("PC1", "PC2"),
#'         variationMetric="cor", nPerm = 10, useSimpleCache=FALSE)
#' 
#' @export

runCOCOAPerm <- function(genomicSignal,
                         signalCoord,
                         GRList,
                         rsScores,
                         sampleLabels,
                         signalCol=c("PC1", "PC2"),
                         signalCoordType = "default",
                         scoringMetric="default",
                         absVal=TRUE,
                         variationMetric="cor",
                         nPerm=300,
                         useSimpleCache=TRUE,
                         cacheDir=getwd(),
                         dataID="tmp",
                         correctionMethod="BH",
                         gammaFitMethod="mme"
                         verbose=TRUE, ...) {
    
    colsToAnnotate = signalCol
    allResultsList = list()
    
    
    indList <- list()
    # generate random indices for shuffling of samples
    for (i in 1:nPerm) {
        indList[[i]] <- sample(1:nrow(sampleLabels), nrow(sampleLabels))
    }

    if (useSimpleCache) {
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
                                   variationMetric = variationMetric,
                                   scoringMetric=scoringMetric,
                                   absVal=absVal)
                    message(i) # must be ahead of object that is saved as cache, not after
                    tmp
                    
                }, cacheDir=cacheDir, ...)
            }
            
            # combining all individual permutations/caches into one object
            for (i in seq_along(indList)) {
                onePermCacheName <- paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID, "_Cache", i)
                rsPermScores[[i]] <- get(onePermCacheName)
            }
            rsPermScores
            
        }, assignToVariable="rsPermScores", cacheDir=cacheDir, ...)    
    } else {
        
        rsPermScores <- list()
        for (i in seq_along(indList)) {
            # for (i in (length(rsPermScores) + 1):nPerm) {
                
            rsPermScores[[i]] <- corPerm(randomInd=indList[[i]],
                                         genomicSignal=genomicSignal,
                                         signalCoord=signalCoord,
                                         GRList=GRList,
                                         calcCols=colsToAnnotate,
                                         sampleLabels=sampleLabels,
                                         variationMetric = variationMetric,
                                         verbose=verbose)
            message(i) # must be ahead of object that is saved as cache, not after
            
        }
    }
    
    
    .analysisID = paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
    
    # # remove region sets that had no overlap
    # keepInd = apply(rsPermScores[[1]], MARGIN = 1, FUN = function(x) !any(is.na(x)))
    # 
    # # screen out region sets with no overlap
    # nullDistList = nullDistList[keepInd]
    # rsScores = rsScores[keepInd, ]
    
    nullDistList = convertToFromNullDist(resultsList=rsPermScores)
    if (useSimpleCache) {

        simpleCache(paste0("permPValsUncorrected", .analysisID), {
            rsPVals = getPermStat(rsScores=rsScores, nullDistList=nullDistList,
                                  calcCols=colsToAnnotate, whichMetric = "pval")
            rsPVals
        }, assignToVariable=rsPVals, cacheDir=cacheDir, ...)
        
        simpleCache(paste0("permZScores", .analysisID), {
            rsZScores = getPermStat(rsScores=rsScores, nullDistList=nullDistList,
                                    calcCols=colsToAnnotate, whichMetric = "zscore")
            rsZScores
            
        }, assignToVariable="rsZScores", cacheDir=cacheDir, ...)
    
        simpleCache(paste0("permPValsCorrected", .analysisID), {
            # p-values based on fitted gamma distributions
            gPValDF = getGammaPVal(scores = rsScores, 
                                   nullDistList = nullDistList, 
                                   calcCols = colsToAnnotate, 
                                   method = "mme", realScoreInDist = TRUE)
            gPValDF = apply(X = gPValDF, MARGIN = 2, 
                            FUN = function(x) p.adjust(p = x, method = correctionMethod))
            gPValDF = cbind(gPValDF, 
                            rsScores[, colnames(rsScores)[!(colnames(rsScores) 
                                                                    %in% colsToAnnotate)]])
            gPValDF
        }, assignToVariable="gPValDF", cacheDir=cacheDir, ...)    
        
    } else {
        rsPVals = getPermStat(rsScores=rsScores, nullDistList=nullDistList,
                              calcCols=colsToAnnotate, whichMetric = "pval")
        
        rsZScores = getPermStat(rsScores=rsScores, nullDistList=nullDistList,
                                calcCols=colsToAnnotate, whichMetric = "zscore")
        
        # p-values based on fitted gamma distributions
        gPValDF = getGammaPVal(scores = rsScores, 
                               nullDistList = nullDistList, 
                               calcCols = colsToAnnotate, 
                               method = "mme", realScoreInDist = TRUE)
        gPValDF = cbind(gPValDF, 
                        rsScores[, colnames(rsScores)[!(colnames(rsScores) 
                                                                %in% colsToAnnotate)]])
        
    }
    
    allResultsList$permRSScores = rsPermScores
    allResultsList$empiricalPVals = rsPVals
    allResultsList$zScores = rsZScores
    allResultsList$gammaPVal = gPValDF
    return(allResultsList)

    # return(named list)
    # definitely: rsPermScores or nullDistList, 
    # also uncorrected perm pvals, z scores, uncorrected gamma pvals
    # possible outputs:  
    # uncorrected pvals perm and gamma, corrected pvals perm and gamma
    # zscores
}

#' @param randomInd numeric. A vector of 1:(number of samples) but shuffled in a
#' random order. E.g. randomInd = sample(1:ncol(genomicSignal), ncol(genomicSignal))
#' where ncol(genomicSignal) is the number of samples. 
#' Set the seed with set.seed() before making randomInd to ensure reproducibility.
#' If the vector is unshuffled,
#' this will give the real COCOA results.
#' @template genomicSignal
#' @template signalCoord
#' @template GRList
#' @param calcCols character. the columns in `sampleLabels` for which to calculate
#' correlation and then to run COCOA on
#' @template sampleLabels
#' @param variationMetric character. 
#' @template scoringMetric
#' @template verbose
#' @template absVal
#' @return data.frame. The output of runCOCOA for one permutation
#' @examples data("brcaMCoord1")
#' data("brcaLoadings1")
#' data("esr1_chr1")
#' data("nrf1_chr1")
#' data("brcaMethylData1")
#' data("brcaPCScores657")
#' sampleLabels = brcaPCScores657[colnames(brcaMethylData1), ]
#' sampleLabels$ER_Status = scale(as.numeric(sampleLabels$ER_Status), 
#'                                center=TRUE, scale=FALSE)
#' # shuffling sample labels
#' randomInd = sample(1:nrow(sampleLabels), nrow(sampleLabels))
#' corPerm(randomInd=randomInd, genomicSignal=brcaMethylData1, 
#'         signalCoord=brcaMCoord1, GRList=GRangesList(esr1_chr1, nrf1_chr1),
#'         calcCols="ER_Status", sampleLabels=sampleLabels, 
#'         variationMetric="cor")
#' 
corPerm <- function(randomInd, genomicSignal, 
                    signalCoord, GRList, calcCols,
                    sampleLabels, variationMetric = "cor", 
                    scoringMetric="default", verbose=TRUE,
                    absVal=TRUE) {
    
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
                           scoringMetric = scoringMetric, verbose = verbose,
                           absVal = absVal)
    
    # return
    return(thisPermRes)
    
}


#' Converts COCOA permutation results to null distributions and vice versa
#' 
#' This function will take a list of results of permutation tests that included
#' many region sets and return a list of data.frames where each data.frame
#' contains the null distribution for a single region set.
#' The function can 
#' also convert in the reverse order from a list of null distributions to a 
#' list of COCOA results. 
#' @param resultsList each item in the list is a data.frame, one item for
#' each permutation with the results of that permutation. Each row in the 
#' data.frame is a region set. All data.frames should be the same size and
#' each data.frame's rows should be in the same order
#' @return a list of data.frames. If given a list where each item is 
#' a data.frame with results from one COCOA permutation, this function
#' will return a list of data.frames where each data.frame contains the
#' null distributions for a single region set. The output data.frames will
#' have the same columns as the input data.frames. If given a list where each
#' item is a data.frame with the null distribution/s for a single region
#' set, this function will return a list where each item is a data.frame
#' with one row for each region set (e.g. a data.frame with results for
#' a single COCOA permutation).
#' 
#' @examples
#' # six region sets (rows), 2 signals (columns)
#' fakePermScores = data.frame(abs(rnorm(6)), abs(rnorm(6)))
#' fakePermScores2 = data.frame(abs(rnorm(6)), abs(rnorm(6)))
#' # 2 fake COCOA results (i.e. nPerm=2)
#' permRSScores = list(fakePermScores, fakePermScores2)
#' convertToFromNullDist(permRSScores)
#' 
#' @export

convertToFromNullDist <- function(resultsList) {

    nullDistList = lapply(X = 1:nrow(resultsList[[1]]),
                          FUN = function(x) permListToOneNullDist(resultsList=resultsList, 
                                                                  rsInd = x))
    return(nullDistList)
}


# @param rsInd numeric. The row number for the region set of interest.
# do for only one region set
permListToOneNullDist <- function(resultsList, rsInd) {
    
    rowList = lapply(resultsList, FUN = function(x) x[rsInd, ])
    rsNullDist = as.data.frame(rbindlist(rowList))
    return(rsNullDist)
}





################################################################################
# p value functions 


#' Get a p-value for region set scores based on a gamma distribution. 
#' 
#' First fit a gamma distribution to each region set's null distribution/s
#' (nullDistList). Then use this gamma distribution to convert scores in
#' rsScores to p-values.
#' 
#' @template rsScores 
#' @param nullDistList list of data.frames. Each list item 
#' has null distributions for a single 
#' region set (list items should be in the same order as rows of rsScores). 
#' Has same score columns as rsScores. 
#' Each column corresponds to a null distribution for that 
#' region set for a given sample variable of interest (e.g. PC or sample phenotype).  
#' @param calcCols character.
#' @param method character. Has the method to use to fit the gamma 
#' distribution to the null distribution.
#' Options are 
#' "mme" (moment matching estimation), "mle" (maximum likelihood estimation), 
#' "qme" (quantile matching estimation), and "mge" (maximum goodness-of-fit 
#' estimation). See ?fitdistrplus::fitdist() for
#' available options and meaning.
#' @param realScoreInDist logical. Should the actual score (from 
#' test with no permutations) be included in the null distribution 
#' when fitting the gamma distribution. realScoreInDist=TRUE is 
#' recommended.
#' @param force logical.
# 
#' @return Returns a data.frame with p values, one column for each col in
#' rsScores 
#' 
#' @examples 
#' fakeOriginalScores = data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' fakePermScores = data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' fakePermScores2 = data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' fakePermScores3 = data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' permRSScores = list(fakePermScores, fakePermScores2, fakePermScores3)
#' nullDistList = convertToFromNullDist(permRSScores)
#' getGammaPVal(rsScores=fakeOriginalScores, nullDistList=nullDistList, calcCols=c("PC1", "PC2")) 
#' 
#' @export

getGammaPVal <- function(rsScores, nullDistList, calcCols, method="mme", realScoreInDist=TRUE, force=FALSE) {
    
    # make sure the same columns are present/in the same order
    
    
    colsToAnnotate = calcCols[calcCols %in% colnames(rsScores)]
    
    if (realScoreInDist) {
        # to get a more accurate gamma distribution, include the score from unpermuted test.
        # add to each null distribution
        for (i in 1:nrow(rsScores)) {
            nullDistList[[i]] = rbind(nullDistList[[i]], as.numeric(rsScores[i, ]))
        }
    }
    
    rsScores = rsScores[, colsToAnnotate, drop=FALSE]  
    
    # returns list, each item in list is also a list.
    # in sub list: each col in nullDistDF has one list item
    fittedDistList = lapply(X = nullDistList, function(x) fitGammaNullDist(nullDistDF = x[, colsToAnnotate, drop=FALSE], 
                                                                            method=method, 
                                                                            force=force))
    
    pValList = list()
    # once for each region set
    for (i in seq_along(nullDistList)) {
        
        pValList[[i]] = as.data.frame(t(pGammaList(scoreVec = as.numeric(rsScores[i, ]), 
                                                   fitDistrList = fittedDistList[[i]])))
    }
    pValDF = as.data.frame(rbindlist(pValList))
    
    
    return(pValDF)
}


# @param nullDistDF a data.frame. Has null distributions for a single 
# region set. Each column corresponds to a null distribution for that 
# region set for a given variable/sample attribute.   
# 
# @return Returns a list or NA. Each list item is a "fitdist" object which is 
# a fitted function. If there are any NA's in nullDistDF, then NA will 
# be returned instead.
# for one of the columns in nullDistDF (output of fitdist() from
# fitdistrplus.
# These are gamma distributions and can be used to get p values for the 
# null distribution so that a large number of permutations 
# are not required.The list is in order of the columns and will
# have the names of the data.frame columns. 

fitGammaNullDist <- function(nullDistDF, method="mme", force=FALSE) {
    
    if (any(is.na(nullDistDF))) {
        
        return(NA)
    }
    
    if (force) {
        # apply returns a list of "fitdist" objects, one list item for each column
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
    
    return(modelList)    #' rsScores <- runCOCOA(signal=brcaLoadings1, 
    #'                                  signalCoord=brcaMCoord1, 
    #'                                  GRList=GRangesList(esr1_chr1), 
    #'                                  signalCol=c("PC1", "PC2"), 
    #'                                  scoringMetric="regionMean")
}

pGammaList <- function(scoreVec, fitDistrList) {
    
    if (any(is.na(scoreVec))) {
        return(rep(NA, length(scoreVec)))
    }
    
    pValVec = mapply(FUN = function(x, y) pgamma(q = y, 
                                                 shape = x$estimate["shape"], 
                                                 rate = x$estimate["rate"], 
                                                 lower.tail = FALSE), x = fitDistrList, y = scoreVec)
    return(pValVec)
}


#' @template rsScores 
#' @param nullDistList list. one item per region set. Each item is a 
#' data.frame with the 
#' null distribution for a single region set. Each column in the data.frame
#' is for a single variable (e.g. PC or latent factor)
#' @param calcCols
#' @param testType character. "greater", "lesser", "two-sided" Whether to
#' create p values based on one sided test or not.
#' @param whichMetric character. Can be "pval" or "zscore"
#' @examples 
#' fakeOriginalScores = data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' fakePermScores = data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' fakePermScores2 = data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' fakePermScores3 = data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' permRSScores = list(fakePermScores, fakePermScores2, fakePermScores3)
#' nullDistList = convertToFromNullDist(permRSScores)
#' getPermStat(rsScores=fakeOriginalScores, nullDistList=nullDistList, 
#'             calcCols=c("PC1", "PC2"), whichMetric="pval") 
#' getPermStat(rsScores=fakeOriginalScores, nullDistList=nullDistList, 
#'             calcCols=c("PC1", "PC2"), whichMetric="zscore") 
#' 
#' @export

getPermStat <- function(rsScores, nullDistList, calcCols, testType="greater", whichMetric = "pval") {
    
    if (is(rsScores, "data.table")) {
        rsScores = as.data.frame(rsScores)
    }
    
    # do once for each region set
    thisStatList = list()
    for (i in 1:nrow(rsScores)) {
        thisStatList[[i]] = as.data.frame(t(getPermStatSingle(rsScore=rsScores[i, calcCols], 
                                                              nullDist = nullDistList[[i]],
                                                              calcCols = calcCols,
                                                              whichMetric=whichMetric)))
        colnames(thisStatList[[i]]) <- calcCols
    }
    thisStat = rbindlist(thisStatList)
    # add back on annotation info
    thisStat = cbind(thisStat, rsScores[, colnames(rsScores)[!(colnames(rsScores) %in% calcCols)]])
    
    # pVals = mapply(FUN = function(x, y) getPermPvalSingle(rsScore=x, 
    #                                            nullDist = y,
    #                                            calcCols = calcCols), 
    #        x = rsScores[, calcCols], y=nullDistList)
    
    
    return(thisStat)
} 


# get p values for a single region set (can get p val for multiple columns)
# @param rsScore a row of values for a single region set. One 
# value for each calcCols
getPermStatSingle <- function(rsScore, nullDist, 
                              calcCols, testType="greater", whichMetric = "pval") {
    
    if (is(nullDist, "data.table")) {
        nullDist = as.data.frame(nullDist)
    }
    
    # if score was NA, return NA
    if (any(is.na(rsScore))) {
        return(rep(NA, length(rsScore)))
    }
    
    
    if (whichMetric == "pval") {
        
        pVal = rep(-1, length(rsScore))
        if (testType == "greater") {
            
            for (i in seq_along(pVal)) {
                # only for one sided test (greater than)
                pVal[i] = 1 - ecdf(x = nullDist[, calcCols[i]])(rsScore[i])
            }
        }
        
        # (-abs(x-0.5) + 0.5) * 2
        
        thisStat = pVal
    }
    
    if (whichMetric == "zscore") {
        
        
        zScore = rep(NA, length(rsScore))
        for (i in seq_along(zScore)) {
            # only for one sided test (greater than)
            zScore[i] = (rsScore[i] - mean(nullDist[, calcCols[i]])) / sd(x = nullDist[, calcCols[i]])
        }
        
        thisStat = as.numeric(zScore)
    }
    
    
    return(thisStat)
}


##############################################################################
############### old ###################

# Run COCOA permutations to get null distributions for each region set.
# For one permutation, there are several steps: 
# First, we shuffle the sample labels. Second, we calculate the association
# between the epigenetic data and the shuffled sample labels using the 
# chosen metric (e.g. correlation). Third, the resulting feature coefficients
# are used as input to the runCOCOA function to score each region set.
# This process is repeated `nPerm` times.
# 
# @param groupByRS logical

# @return A list where each item is a data.frame with the null 
# distribution for a single region set. The length of the list
# is equal to the number of region sets. The number of rows of 
# each data.frame is equal to the number of permutations.
# getNullDist <- function(groupByRS=TRUE) {
#     
# }