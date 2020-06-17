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
#' See these individual functions for more info on each step: runCOCOA, 
#' convertToFromNullDist, getPermStat, and getGammaPVal. 
#' 
#' 
#' For reproducibility, set seed with 'set.seed()' function before running.
#' @param nPerm Numeric. The number of permutations to do.
#' @template genomicSignal
#' @template signalCoord
#' @template GRList
#' @template rsScores
#' @template targetVar
#' @templateVar usesTargetVar TRUE
#' @template signalCol
#' @template scoringMetric
#' @template absVal
#' @template olList
#' @param centerGenomicSignal Logical. Should rows in genomicSignal
#' be centered based on
#' their means? (subtracting row mean from each row)
#' @param centerTargetVar Logical. Should columns in targetVar be 
#' centered based
#' on their means? (subtract column mean from each column)
#' @param dataID Character. A unique identifier for this dataset 
#' (for saving results with simpleCache)
#' @template variationMetric
#' @param useSimpleCache Logical. Whether to use save caches. Caches
#' will be created for each permutation so that if the function is disrupted
#' it can restart where it left off. The final results are also saved 
#' as a cache. See simpleCache package for more details.
#' @param cacheDir Character. The path for the directory in which the
#' caches should be saved. 
# @param correctionMethod Character. P value correction method. Default
# is "BH" for Benjamini and Hochberg false discovery rate. For acceptable 
# arguments and more info see ?stats::p.adjust() (method parameter) 
#' @param testType Character. Parameter for `getPermStat`. Whether to
#' create p values based on one a two sided test or a lesser/greater one
#' sided test. Options are: "greater", "lesser", "two-sided" 
#' @param gammaFitMethod Character. method to use for fitting the gamma
#' distribution to null distribution. Options are 
#' "mme" (moment matching estimation), "mle" (maximum likelihood estimation), 
#' "qme" (quantile matching estimation), and "mge" (maximum goodness-of-fit 
#' estimation). See ?COCOA::getGammaPVal and 
#' ?fitdistrplus::fitdist() for more info.
#' @param realScoreInDist Logical. Should the actual score (from 
#' test with no permutations) be included in the null distribution 
#' when fitting the gamma distribution. realScoreInDist=TRUE is 
#' recommended.
#' @param force Logical. If force=TRUE, when fitting the gamma distribution
#' returns an error (as may happen when a method other than "mme"
#' is used) then allow the error. If force=FALSE, when fitting the 
#' gamma distribution returns an error then don't return an error but 
#' instead use the "mme" method
#' for fitting that specific gamma distribution.
#' @template verbose
#' @template returnCovInfo
#' @param ... Character. Optional additional arguments for simpleCache.
#'
#' 
#' @return Returns a list with the following 4 items: 1. a list of length nPerm
#' where each item is a data.frame of the COCOA scores from a single 
#' permutation. Each data.frame is the output of `runCOCOA()` 
#' 2. a data.table/data.frame of empirical p-values (the
#' output of `getPermStat`) 3. a 
#' data.table/data.frame of z-scores (the output of `getPermStat`. 
#' 4. a data.frame of p-values based on
#' the gamma approximation (the output of getGammaPVal(). 
#' @examples 
#' data("esr1_chr1")
#' data("nrf1_chr1")
#' data("brcaMethylData1")
#' data("brcaMCoord1")
#' pcScores <- prcomp(t(brcaMethylData1))$x
#' targetVarCols <- c("PC1", "PC2")
#' targetVar <- pcScores[, targetVarCols]
#' 
#' # give the actual order of samples to `runCOCOA` to get the real scores
#' correctSampleOrder=1:nrow(targetVar)
#' realRSScores <- runCOCOA(genomicSignal=brcaMethylData1,
#'                         signalCoord=brcaMCoord1,
#'                         GRList=GRangesList(esr1_chr1, nrf1_chr1),
#'                         signalCol=c("PC1", "PC2"),
#'                         targetVar=targetVar,
#'                         sampleOrder=correctSampleOrder,
#'                         variationMetric="cor")
#'         
#' # give random order of samples to get random COCOA scores 
#' # so you start building a null distribution for each region set 
#' # (see vignette for example of building a null distribution with `runCOCOA`)
#' randomOrder <- sample(1:nrow(targetVar), 
#'                       size=nrow(targetVar),
#'                       replace=FALSE)
#' randomRSScores <- runCOCOA(genomicSignal=brcaMethylData1,
#'                           signalCoord=brcaMCoord1,
#'                           GRList=GRangesList(esr1_chr1, nrf1_chr1),
#'                           signalCol=c("PC1", "PC2"),
#'                           targetVar=targetVar,
#'                           sampleOrder=randomOrder,
#'                           variationMetric="cor")
#' 
#' # runCOCOAPerm
#' permResults <- runCOCOAPerm(genomicSignal=brcaMethylData1,
#'                            signalCoord=brcaMCoord1,
#'                            GRList=GRangesList(esr1_chr1, nrf1_chr1),
#'                            rsScores=realRSScores,
#'                            targetVar=targetVar,
#'                            signalCol=c("PC1", "PC2"),
#'                            variationMetric="cor",
#'                            nPerm = 10,
#'                            useSimpleCache=FALSE)
#' permResults
#'   
#' 
#' @export

runCOCOAPerm <- function(genomicSignal,
                         signalCoord,
                         GRList,
                         rsScores,
                         targetVar,
                         signalCol=c("PC1", "PC2"),
                         scoringMetric="default",
                         absVal=TRUE,
                         olList=NULL,
                         centerGenomicSignal=TRUE,
                         centerTargetVar=TRUE,
                         variationMetric="cor",
                         nPerm=300,
                         useSimpleCache=TRUE,
                         cacheDir=getwd(),
                         dataID="",
                         testType="greater",
                         gammaFitMethod="mme",
                         realScoreInDist=TRUE,
                         force=FALSE,
                         verbose=TRUE,
                         returnCovInfo=FALSE, ...) {
    
    
    colsToAnnotate <- signalCol
    allResultsList <- list()
    
    if (is(rsScores, "data.table")) {
        rsScores = as.data.frame(rsScores)
    }
    rsScores = rsScores[, colsToAnnotate] # prevents error that occurs if extra column is factor
    
    # more efficient to only do once (not that high impact though)
    if (centerGenomicSignal) {
        cpgMeans <- rowMeans(genomicSignal, na.rm = TRUE)
        # centering before calculating correlation
        genomicSignal <- apply(X = genomicSignal, MARGIN = 2, function(x) x - cpgMeans)
        # don't do later
        centerGenomicSignal <- FALSE
    }
    if (centerTargetVar) {
        featureMeans <- colMeans(targetVar, na.rm = TRUE)
        # centering before calculating correlation (also, t() converts to matrix)
        targetVar <- t(apply(X = t(targetVar), MARGIN = 2, function(x) x - featureMeans))
        if (dim(targetVar)[1] == 1) {
            targetVar <- t(targetVar)
        }
        # don't do later
        centerTargetVar <- FALSE
    }
    

    checkConvertInputClasses(signal=genomicSignal,
                             signalCoord=signalCoord,
                             regionSet=NULL,
                             signalCol=signalCol,
                             rsOL=NULL)
    
        # detect signalCoordType
        # when signalCoord is a GRanges object
        if (any(start(signalCoord) != end(signalCoord))) {
            signalCoordType <- "multiBase"
        } else {
            signalCoordType <- "singleBase"
        }

    
    # if "default" scoring method is given, choose based on signalCoordType
    if (scoringMetric == "default") {
        if (signalCoordType == "singleBase") {
            scoringMetric <- "regionMean"   
        } else if (signalCoordType == "multiBase") {
            scoringMetric <- "proportionWeightedMean"
        } else {
            stop(cleanws("signalCoordType not recognized. 
                         Check spelling/capitalization."))
        }
    }
    
    #################
    if (is.null(olList)) {
        
        #######
        # must take out NA rows before getting OL list. Otherwise later calculations 
        # will use wrong indices. 
        
        # what happens if there are NAs or Inf in `signal`?
        # any NAs that overlap the regionSet will cause the score to be NA
        if (is(genomicSignal, "data.table")) {
            naRows = apply(X = genomicSignal[, , with=FALSE, drop=FALSE], 
                           MARGIN = 1, FUN = function(x) any(is.na(x)))
        } else {
            naRows = apply(X = genomicSignal[, , drop=FALSE], 
                           MARGIN = 1, FUN = function(x) any(is.na(x)))    
        }
        
        if (any(naRows)) {
            genomicSignal <- genomicSignal[!naRows, ]
            signalCoord <- signalCoord[!naRows]
            warning("Removing rows with NA from `genomicSignal`")
        }
        
        #################################################################
        
        # calculate overlaps only once
        # region set must be subject to fit with scoring functions
        olList <- lapply(X = GRList, FUN = function(x) findOverlaps(query = signalCoord, 
                                                                 subject = x))
        totalRegionNumber = sapply(X = GRList, length)
        meanRegionSize = sapply(X = GRList, function(x) round(mean(width(x))))
    }
    
    # also calculate coverage info
    # @param rsOL 
    calculateCovInfo <- function(rsOL, 
                                 scoringMetric=scoringMetric, 
                                 pOlap=NULL) {
        
        covInfo <- data.frame(signalCoverage=length(unique(queryHits(rsOL))), 
                              regionSetCoverage=length(unique(subjectHits(rsOL))))
        
        if (scoringMetric == "proportionWeightedMean") {
            
            covInfo$sumProportionOverlap <- sum(pOlap)
        }
        
    }
    
    covInfo <- lapply(X = olList, FUN = function(x) calculateCovInfo(rsOL=x, 
                                                          scoringMetric = scoringMetric))
    
    
    if (scoringMetric == "proportionWeightedMean") {
        getPOlap <- function(rsOL, signalGR, regionSet) {
            olap  <- pintersect(signalGR[queryHits(rsOL)],
                                regionSet[subjectHits(rsOL)])
            pOlap <- width(olap) / width(regionSet[subjectHits(rsOL)])
            return(pOlap)
        }
        
        # list
        pOlapList <- mapply(FUN = function(x, y) getPOlap(rsOL = x, 
                                                     signalGR=signalCoord, 
                                                     regionSet = y), 
                       x=olList, y=GRList, SIMPLIFY = FALSE)
        covInfo$sumProportionOverlap <- sapply(X = pOlapList, FUN = sum)
    } else {
        pOlapList <- NULL
    }
    
    #################
    
        
    indList <- list()
    # generate random indices for shuffling of samples
    for (i in 1:nPerm) {
        indList[[i]] <- sample(1:nrow(targetVar), nrow(targetVar))
    }
    # replicate(10, sample(seq_len(nrow(targetVar)))

    if (useSimpleCache) {
        
        # create the main permutation cache
        simpleCache(paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), {
            
        helperFun <- function(x, y, ...) {
            # for (i in (length(rsPermScores) + 1):nPerm) {
            onePermCacheName <- paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID, "_Cache", y)
            # create sub caches, one for each permutation
            simpleCache(onePermCacheName, cacheSubDir = paste0("rsPermScores_", nPerm, "Perm_", variationMetric, "_", dataID), {
                
                tmp <- runCOCOA(sampleOrder=x,
                                genomicSignal=genomicSignal,
                                signalCoord=signalCoord,
                                GRList=GRList,
                                signalCol=colsToAnnotate,
                                targetVar=targetVar,
                                variationMetric = variationMetric,
                                scoringMetric=scoringMetric,
                                absVal=absVal,
                                centerGenomicSignal = centerGenomicSignal,
                                centerTargetVar = centerTargetVar,
                                verbose=verbose,
                                olList=olList,
                                pOlapList=pOlapList,
                                returnCovInfo = returnCovInfo)
                message(y) # must be ahead of object that is saved as cache, not after
                tmp
                
            }, cacheDir=cacheDir, assignToVariable="tmp", ...)
            return(tmp)
        }
        rsPermScores = mapply(FUN = helperFun, x=indList, y=seq_along(indList), ..., SIMPLIFY=FALSE)
        rsPermScores
            
        }, assignToVariable="rsPermScores", cacheDir=cacheDir, ...)    
    } else {
        
        helperFun <- function(x) {
            tmp <- runCOCOA(sampleOrder=x,
                          genomicSignal=genomicSignal,
                          signalCoord=signalCoord,
                          GRList=GRList,
                          signalCol=colsToAnnotate,
                          targetVar=targetVar,
                          variationMetric = variationMetric,
                          scoringMetric=scoringMetric,
                          absVal=absVal,
                          verbose=verbose,
                          olList=olList,
                          pOlapList=pOlapList,
                          returnCovInfo = returnCovInfo)
            message(".")
            return(tmp)
        }
        
        rsPermScores <- lapply(X = indList, helperFun)
 
    }
    
    
    .analysisID <- paste0("_", nPerm, "Perm_", variationMetric, "_", dataID)
    
    # # remove region sets that had no overlap
    # keepInd <- apply(rsPermScores[[1]], MARGIN = 1, FUN = function(x) !any(is.na(x)))
    # 
    # # screen out region sets with no overlap
    # nullDistList <- nullDistList[keepInd]
    # rsScores <- rsScores[keepInd, ]
    
    nullDistList <- convertToFromNullDist(rsPermScores)
    if (useSimpleCache) {

        simpleCache(paste0("empiricalPValsUncorrected", .analysisID), {
            rsPVals <- getPermStat(rsScores=rsScores, nullDistList=nullDistList,
                                  signalCol=colsToAnnotate, whichMetric = "pval",
                                  testType=testType)
            rsPVals
        }, assignToVariable="rsPVals", cacheDir=cacheDir, ...)
        
        simpleCache(paste0("permZScores", .analysisID), {
            rsZScores <- getPermStat(rsScores=rsScores, nullDistList=nullDistList,
                                    signalCol=colsToAnnotate, whichMetric = "zscore")
            rsZScores
            
        }, assignToVariable="rsZScores", cacheDir=cacheDir, ...)
    
        simpleCache(paste0("gammaPValsUncorrected", .analysisID), {
            # p-values based on fitted gamma distributions
            gPValDF <- getGammaPVal(rsScores = rsScores, 
                                   nullDistList = nullDistList, 
                                   signalCol = colsToAnnotate, 
                                   method = gammaFitMethod, 
                                   realScoreInDist = realScoreInDist,
                                   force=force)
            # gPValDF <- apply(X = gPValDF, MARGIN = 2, 
            #                 FUN = function(x) p.adjust(p = x, method = correctionMethod))
            gPValDF <- cbind(gPValDF, 
                            rsScores[, colnames(rsScores)[!(colnames(rsScores) 
                                                                    %in% colsToAnnotate)]])
            gPValDF
        }, assignToVariable="gPValDF", cacheDir=cacheDir, ...)    
        
    } else {
        rsPVals <- getPermStat(rsScores=rsScores, nullDistList=nullDistList,
                              signalCol=colsToAnnotate, whichMetric = "pval",
                              testType = testType)
        
        rsZScores <- getPermStat(rsScores=rsScores, nullDistList=nullDistList,
                                signalCol=colsToAnnotate, whichMetric = "zscore")
        
        # p-values based on fitted gamma distributions
        gPValDF <- getGammaPVal(rsScores = rsScores, 
                               nullDistList = nullDistList, 
                               signalCol = colsToAnnotate, 
                               method = gammaFitMethod, 
                               realScoreInDist = realScoreInDist,
                               force=force)
        gPValDF <- cbind(gPValDF, 
                        rsScores[, colnames(rsScores)[!(colnames(rsScores) 
                                                                %in% colsToAnnotate)]])
        
    }
    
    allResultsList$permRSScores <- rsPermScores
    allResultsList$empiricalPVals <- rsPVals
    allResultsList$zScores <- rsZScores
    allResultsList$gammaPVal <- gPValDF
    return(allResultsList)

    # return(named list)
    # definitely: rsPermScores or nullDistList, 
    # also uncorrected perm pvals, z scores, uncorrected gamma pvals
    # possible outputs:  
    # uncorrected pvals perm and gamma, corrected pvals perm and gamma
    # zscores
}

#' Run COCOA: quantify inter-sample variation, score region sets
#' 
#' This is a convenience function that does the two steps of COCOA: 
#' quantifying the epigenetic variation and scoring the region sets. 
#' This function will return the real COCOA scores if using the default
#' `sampleOrder` parameter values. This
#' function also makes it easy to generate null distributions in order to
#' evaluate the statistical significance of the real COCOA results.
#' You can use the sampleOrder parameter to shuffle the samples,
#' then run COCOA to get fake scores for each region set. By doing 
#' this many times, you can build a null distribution for each 
#' region set composed of the region set's random scores from each
#' permutation. There are multiple options for quantifying the
#' epigenetic variation, specified by the `variationMetric` parameter.
#' Quantifying the variation for the real/non-permuted COCOA 
#' scores should be done with the same 
#' variation metric as is used for the random permutations. For an
#' unsupervised analysis using dimensionality reduction, first, the
#' dimensionality reduction is done outside `runCOCOA`, then the
#' latent factors/principal components are input to `runCOCOA` as the
#' sample labels (targetVar parameter) when calculating both the real and 
#' also the permutated region set scores. For a supervised analysis, 
#' the target variables/phenotypes are the targetVar.
#' See the vignettes for examples.  
#' 
#' @param sampleOrder numeric. A vector of length (number of samples). If
#' sampleOrder is 1:(number of samples) then this function will return the
#' real COCOA scores.
#' To generate random COCOA scores in order to make 
#' null distributions, shuffle the samples in a random order.
#' E.g. sampleOrder = sample(1:ncol(genomicSignal), ncol(genomicSignal))
#' where ncol(genomicSignal) is the number of samples. 
#' Set the seed with set.seed() before making sampleOrder to ensure reproducibility.
#' @template genomicSignal
#' @template signalCoord
#' @template GRList
#' @templateVar usesTargetVar TRUE
#' @template signalCol
#' @template targetVar
#' @template variationMetric
#' @template scoringMetric
#' @template verbose
#' @template absVal
#' @template olList
#' @template pOlapList
#' @param centerGenomicSignal Logical. Should rows in genomicSignal
#' be centered based on
#' their means? (subtracting row mean from each row)
#' @param centerTargetVar Logical. Should columns in targetVar be 
#' centered based
#' on their means? (subtract column mean from each column)
#' @template returnCovInfo
#' @return data.frame. The output of aggregateSignalGRList for one permutation.
#' @examples
#' data("esr1_chr1")
#' data("nrf1_chr1")
#' data("brcaMethylData1")
#' data("brcaMCoord1")
#' pcScores <- prcomp(t(brcaMethylData1))$x
#' targetVarCols <- c("PC1", "PC2")
#' targetVar <- pcScores[, targetVarCols]
#' 
#' # give the actual order of samples to `runCOCOA` to get the real scores
#' correctSampleOrder=1:nrow(targetVar)
#' realRSScores <- runCOCOA(genomicSignal=brcaMethylData1,
#'                         signalCoord=brcaMCoord1,
#'                         GRList=GRangesList(esr1_chr1, nrf1_chr1),
#'                         signalCol=c("PC1", "PC2"),
#'                         targetVar=targetVar,
#'                         sampleOrder=correctSampleOrder,
#'                         variationMetric="cor")
#' realRSScores
#'         
#' # give random order of samples to get random COCOA scores 
#' # so you start building a null distribution for each region set 
#' # (see vignette for example of building a null distribution with `runCOCOA`)
#' randomOrder <- sample(1:nrow(targetVar), 
#'                       size=nrow(targetVar),
#'                       replace=FALSE)
#' randomRSScores <- runCOCOA(genomicSignal=brcaMethylData1,
#'                           signalCoord=brcaMCoord1,
#'                           GRList=GRangesList(esr1_chr1, nrf1_chr1),
#'                           signalCol=c("PC1", "PC2"),
#'                           targetVar=targetVar,
#'                           sampleOrder=randomOrder,
#'                           variationMetric="cor")
#' randomRSScores
#' @export
runCOCOA <- function(genomicSignal, 
                    signalCoord, GRList, signalCol,
                    targetVar, 
                    sampleOrder=1:nrow(targetVar),
                    variationMetric = "cor", 
                    scoringMetric="default", verbose=TRUE,
                    absVal=TRUE, olList=NULL, pOlapList=NULL,
                    centerGenomicSignal=TRUE,
                    centerTargetVar=TRUE, 
                    returnCovInfo=TRUE) {
    
    # if vector is given, return error
    if (is.null(dim(targetVar))) {
        stop("`targetVar` should be a matrix or data.frame")
    }
    
    if (any(!(signalCol %in% colnames(targetVar)))) {
        stop("Not all specified columns are present in `targetVar`")
    }
    
    
    # subset to only signalCol
    targetVar <- targetVar[, signalCol, drop=FALSE]
    
    # because names are dropped for a single column data.frame when indexing
    # single col data.frame is automatically converted to numeric
    featureNames <- colnames(targetVar)
    # reorder the sample labels
    targetVar <- data.frame(targetVar[sampleOrder, ])
    colnames(targetVar) <- featureNames
    
    # calculate correlation
    featureLabelCor <- createCorFeatureMat(dataMat = genomicSignal, 
                                          featureMat = targetVar, 
                                          centerDataMat = centerGenomicSignal, 
                                          centerFeatureMat = centerTargetVar,
                                          testType = variationMetric)
    
    # more efficient to do only once instead of for each region set later on
    if (absVal) {
        featureLabelCor <- abs(featureLabelCor)
        absVal <- FALSE    
    }
    
    
    # run COCOA
    thisPermRes <- aggregateSignalGRList(signal=featureLabelCor, 
                           signalCoord=signalCoord, GRList=GRList, 
                           signalCol = signalCol, 
                           scoringMetric = scoringMetric, verbose = verbose,
                           absVal = absVal, olList = olList, pOlapList=pOlapList,
                           returnCovInfo=returnCovInfo)
    
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
#' @param rsScoresList each item in the list is a data.frame, one item for
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
#' fakePermScores <- data.frame(abs(rnorm(6)), abs(rnorm(6)))
#' fakePermScores2 <- data.frame(abs(rnorm(6)), abs(rnorm(6)))
#' # 2 fake COCOA results (i.e. nPerm=2)
#' permRSScores <- list(fakePermScores, fakePermScores2)
#' convertToFromNullDist(permRSScores)
#' 
#' @export

convertToFromNullDist <- function(rsScoresList) {

    nullDistList <- lapply(X = 1:nrow(rsScoresList[[1]]),
                          FUN = function(x) permListToOneNullDist(resultsList=rsScoresList, 
                                                                  rsInd = x))
    return(nullDistList)
}


# @param rsInd numeric. The row number for the region set of interest.
# do for only one region set
permListToOneNullDist <- function(resultsList, rsInd) {
    
    rowList <- lapply(resultsList, FUN = function(x) x[rsInd, ])
    rsNullDist <- as.data.frame(rbindlist(rowList))
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
#' region set for a given sample variable of interest/target variable
#' (e.g. PC or sample phenotype).  
#' @templateVar usesRSScores TRUE
#' @template signalCol
#' @param method Character. Has the method to use to fit the gamma 
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
#' @param force logical. If force=TRUE, when fitting the gamma distribution
#' returns an error (as may happen when a method other than "mme"
#' is used) then allow the error. If force=FALSE, when fitting the 
#' gamma distribution returns an error then don't return an error but 
#' instead use the "mme" method
#' for fitting that specific gamma distribution.
# 
#' @return Returns a data.frame with p values, one column for each signalCol in
#' rsScores 
#' 
#' @examples 
# fakeOriginalScores <- data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
# fakePermScores <- data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
# fakePermScores2 <- data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
# fakePermScores3 <- data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
# permRSScores <- list(fakePermScores, fakePermScores2, fakePermScores3)
# nullDistList <- convertToFromNullDist(permRSScores)
# getGammaPVal(rsScores=fakeOriginalScores,
#              nullDistList=nullDistList,
#              signalCol=c("PC1", "PC2"))
#' 
#' @export

getGammaPVal <- function(rsScores, nullDistList, signalCol, method="mme", 
                         realScoreInDist=TRUE, force=FALSE) {
    
    
    # take absolute value since gamma distribution cannot have negative values
    if (any(rsScores < 0, na.rm = TRUE)) {
        rsScores <- abs(rsScores)
    }
    conditionalAbs <- function(x) {
        if (any(x < 0, na.rm = TRUE)) {
            return(abs(x))
        } else {
            return(x)
        }
    }
    nullDistList <- lapply(X = nullDistList, FUN = conditionalAbs)
    
    
    # make sure the same columns are present/in the same order
    
    
    colsToAnnotate <- signalCol[signalCol %in% colnames(rsScores)]
    
    if (realScoreInDist) {
        # to get a more accurate gamma distribution, include the score from unpermuted test.
        # add to each null distribution
        for (i in 1:nrow(rsScores)) {
            nullDistList[[i]] <- rbind(nullDistList[[i]], as.numeric(rsScores[i, ]))
        }
    }
    
    rsScores <- rsScores[, colsToAnnotate, drop=FALSE]  
    
    # returns list, each item in list is also a list.
    # in sub list: each col in nullDistDF has one list item
    fittedDistList <- lapply(X = nullDistList, 
                             function(x) fitGammaNullDist(nullDistDF = 
                                                              x[, colsToAnnotate, 
                                                                drop=FALSE], 
                                                          method=method, 
                                                          force=force))
    
    # region sets with no coverage have NA for fittedDistList value
    # pGammaList cannot assign these column names so assign them ahead of time
    # and don't run pGammaList on those region sets
    naInd = vapply(X = fittedDistList, 
                   FUN = function(x) all(is.na(x)), 
                   FUN.VALUE = TRUE)
    
    pValList <- list()
    
    # assign na items
    naEntry = as.data.frame(t(rep(NA, length(colsToAnnotate))))
    colnames(naEntry) = colsToAnnotate
    pValList[naInd] = list(naEntry)
    
    # once for each region set
    for (i in seq_along(nullDistList)[!naInd]) {
        
        pValList[[i]] <- as.data.frame(t(pGammaList(scoreVec = as.numeric(rsScores[i, ]), 
                                                   fitDistrList = fittedDistList[[i]])))
    }
    pValDF <- as.data.frame(rbindlist(pValList))
    
    
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
        modelList <- apply(X = nullDistDF, 
                          MARGIN = 2, 
                          FUN = function(x) fitdistrplus::fitdist(data=x, 
                                                                  distr = "gamma", 
                                                                  method=method))
    } else {
        # try "method". if it fails, do method = "mme"
        modelList <- apply(X = nullDistDF, 
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
    
    if (any(is.na(scoreVec))) {
        naVec <- rep(NA, length(scoreVec))
        names(naVec) <- names(scoreVec)
        return(naVec)
    }
    
    pValVec <- mapply(FUN = function(x, y) pgamma(q = y, 
                                                 shape = x$estimate["shape"], 
                                                 rate = x$estimate["rate"], 
                                                 lower.tail = FALSE), x = fitDistrList, y = scoreVec)
    return(pValVec)
}

#' Get p-value or z-score based on permutation results
#' 
#' This function starts with real COCOA scores for each
#' region set and null distributions for each
#' region set that come
#' from running COCOA on permuted data. Then this function uses the
#' null distributions to get an empirical p-value or z-score for
#' each region set. See vignettes for the workflow that leads to
#' this function. The calculation of the p-value/z-score does not 
#' include the real region set score in the null distribution.
#'
#' @template rsScores 
#' @param nullDistList List. one item per region set. Each item is a 
#' data.frame with the 
#' null distribution/s for a single region set. Each column in the data.frame
#' is for a target variable (e.g. PC or phenotype), which is given
#' by the `signalCol` parameter (each target variable has a different
#' null distribution for a given region set).
#' @templateVar usesRSScores
#' @template signalCol 
#' @param testType Character. "greater", "lesser", "two-sided" Whether to
#' create p values based on one sided test or not. Only applies when
#' whichMetric="pval".
#' @param whichMetric Character. Can be "pval" or "zscore"
#' @return A data.table/data.frame. 
#' If whichMetric="pval", returns the empirical p-value for
#' each region set in `rsScores`. If the region set score is more extreme
#' than all scores in the null distribution, a p-value of 0 is returned but
#' this simply means the p-value is the minimum detectable p-value with
#' the given number of permutations used to make the null distributions. If
#' whichMetric="zscore", the function returns a z-score for each region set
#' score: ((region set score) - mean(null distribution)) / sd(null distribution)
#' @examples 
#' fakeOriginalScores <- data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' fakePermScores <- data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' fakePermScores2 <- data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' fakePermScores3 <- data.frame(PC1=abs(rnorm(6)), PC2=abs(rnorm(6)))
#' permRSScores <- list(fakePermScores, fakePermScores2, fakePermScores3)
#' nullDistList <- convertToFromNullDist(permRSScores)
#' getPermStat(rsScores=fakeOriginalScores, nullDistList=nullDistList, 
#'             signalCol=c("PC1", "PC2"), whichMetric="pval") 
#' getPermStat(rsScores=fakeOriginalScores, nullDistList=nullDistList, 
#'             signalCol=c("PC1", "PC2"), whichMetric="zscore") 
#' 
#' @export

getPermStat <- function(rsScores, nullDistList, signalCol, 
                        testType="greater", whichMetric = "pval") {
    
    if (is(rsScores, "data.table")) {
        rsScores <- as.data.frame(rsScores)
    }
    
    # do once for each region set
    thisStatList <- list()
    for (i in 1:nrow(rsScores)) {
        thisStatList[[i]] <- as.data.frame(t(getPermStatSingle(rsScore=rsScores[i, signalCol], 
                                                              nullDist = nullDistList[[i]],
                                                              signalCol = signalCol,
                                                              whichMetric=whichMetric)))
        colnames(thisStatList[[i]]) <- signalCol
    }
    thisStat <- rbindlist(thisStatList)
    # add back on annotation info
    thisStat <- cbind(thisStat, rsScores[, colnames(rsScores)[!(colnames(rsScores) %in% signalCol)]])
    
    # pVals <- mapply(FUN = function(x, y) getPermPvalSingle(rsScore=x, 
    #                                            nullDist = y,
    #                                            signalCol = signalCol), 
    #        x = rsScores[, signalCol], y=nullDistList)
    
    
    return(thisStat)
} 


# get p values for a single region set (can get p val for multiple columns)
# @param rsScore a row of values for a single region set. One 
# value for each signalCol
getPermStatSingle <- function(rsScore, nullDist, 
                              signalCol, testType="greater", whichMetric = "pval") {
    
    if (is(nullDist, "data.table")) {
        nullDist <- as.data.frame(nullDist)
    }
    
    # if score was NA, return NA
    if (any(is.na(rsScore))) {
        return(rep(NA, length(rsScore)))
    }
    
    
    if (whichMetric == "pval") {
        
        pVal <- rep(-1, length(rsScore))
        if (testType == "greater") {
            for (i in seq_along(pVal)) {
                # only for one sided test (greater than)
                pVal[i] <- 1 - ecdf(x = nullDist[, signalCol[i]])(rsScore[i])
            }
        } else if (testType == "lesser") {
            for (i in seq_along(pVal)) {
                # only for one sided test (less than)
                pVal[i] <- ecdf(x = nullDist[, signalCol[i]])(rsScore[i])
            }
        } else if (testType == "two-sided") {
            
            for (i in seq_along(pVal)) {
                cdfVal <- ecdf(x = nullDist[, signalCol[i]])(rsScore[i])
                # convert to two sided pval
                pVal[i] <- -2 * abs(cdfVal - 0.5) + 1
            }

        }
        
        thisStat <- pVal
    }
    
    if (whichMetric == "zscore") {
        
        zScore <- rep(NA, length(rsScore))
        for (i in seq_along(zScore)) {
            # greater than/less than hypothesis does not matter for z-score
            zScore[i] <- (rsScore[i] - mean(nullDist[, signalCol[i]])) / sd(x = nullDist[, signalCol[i]])
        }
        
        thisStat <- as.numeric(zScore)
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
# are used as input to the aggregateSignalGRList function to score each region set.
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