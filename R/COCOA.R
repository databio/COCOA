# PACKAGE DOCUMENTATION
#' Coordinate Covariation Analysis (COCOA)
#'
#' 
#' COCOA is a method for understanding variation among samples.
#' COCOA can be used with data that includes 
#' genomic coordinates such as DNA methylation. 
#' To describe the method on a high level, COCOA uses a database of 
#' "region sets" and principal component analysis (PCA) of your data 
#' to identify sources of variation among samples. A region set is a set of 
#' genomic regions that share a biological annotation, 
#' for instance transcription factor (TF) binding regions, 
#' histone modification regions, or open chromatin regions. 
#' In contrast to some other common techniques, COCOA is unsupervised, 
#' meaning that samples do not have to be divided into groups 
#' such as case/control or healthy/disease, although COCOA works in 
#' those situations as well. Also, COCOA focuses on continuous variation 
#' between samples instead of having cutoffs. Because of this, COCOA can 
#' be used as a complementary method alongside "differential" methods 
#' that find discrete differences between groups of samples and 
#' it can also be used in situations where there are no groups.  
#' COCOA can identify biologically meaningful 
#' sources of variation between samples and increase understanding of 
#' variation in your data. 
#'
#' @docType package
#' @name COCOA
#' @author John Lawson
#' @author Nathan Sheffield
#'
#' @references \url{http://github.com/databio}
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_jitter geom_line
#'             theme_classic xlab ylab geom_hline ylim scale_color_discrete
#'             scale_x_discrete scale_fill_brewer scale_color_manual
#'             scale_color_brewer
#' @importFrom ComplexHeatmap Heatmap draw
#' @import BiocGenerics S4Vectors IRanges GenomicRanges simpleCache
#' @importFrom data.table ":=" setDT data.table setkey fread setnames 
#'             setcolorder rbindlist setattr setorder copy is.data.table 
#'             setorderv as.data.table
#' @importFrom Biobase sampleNames
#' @importFrom stats lm coefficients poly wilcox.test ecdf pgamma p.adjust
#' @importFrom methods is
#' @importFrom MIRA binRegion
#' @importFrom tidyr gather
#' @importFrom grid grid.newpage grid.grabExpr grid.draw popViewport 
#'             pushViewport unit viewport
#' @importFrom grDevices dev.off
#' @importFrom methods hasArg
#' @importFrom fitdistrplus fitdist
#' @importFrom simpleCache simpleCache
#' @importFrom ppcor pcor
NULL

# now package lists GenomicRanges in "Depends" instead of "Imports" in 
# DESCRIPTION, still import package with @import though
# @importFrom GenomicRanges GRanges GRangesList elementMetadata strand
#             seqnames granges

# Because of some issues, 
# (see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation, 
# in order to pass some R CMD check NOTES.
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        ".", "..calcCols", "bin", "binID", "chr", "id", "colsToAnnotate", 
        "coordinateDT",
        "coverage", "pOlap", "regionGroupID", "regionID", "theme", 
        "meanRegionSize", "regionSetCoverage", "rowIndex", "rsIndex",
        "rsRegionID", "totalRegionNumber", "signalCoverage", ".SD",
        "sumProportionOverlap")) 
}

#########################################################################


#' Use PCA loadings to score a region set
#' 
#' First, this function identifies which loadings are within the region set. 
#' Then the loadings are used to score the region set 
#' according to the `scoringMetric` parameter.  
#' 
#' @param signal matrix of loadings (the coefficients of 
#' the linear combination that defines each PC). One named column for each PC.
#' One row for each original dimension/variable (should be same order 
#' as original data/signalCoord). The x$rotation output of prcomp().
#' @template signalCoord
#' @template signalCoordType
#' @templateVar refGenomeWarning TRUE
#' @template regionSet 
#' @template signalCol
#' @template scoringMetric
#' @template verbose
# Useful when using 
# aggregateSignal with 'apply' to do many region sets at a time.
# @param wilcox.conf.int logical. Only applies when using "rankSum" scoring
# method. returns a 95% confidence interval from the Wilcoxon rank sum test
# instead of p value.
#' @template absVal

#' @return a data.frame with one row and the following 
#' columns: one column for each item of signalCol with names given
#' by signalCol. These columns have scores for the region set for each PC.
#' Other columns: signalCoverage (formerly cytosine_coverage) which
#' has number of cytosines that overlapped with regionSet (or for "multiBase"
#' data, 
#' the number of regions from signalCoord that overlapped regionSet) 
#' regionSetCoverage which has number of regions from regionSet
#' that overlapped any coordinates from signalCoord, 
#' totalRegionNumber that has
#' number of regions in regionSet, meanRegionSize that has average
#' size in base pairs of regions in regionSet, the average is based on
#' all regions in regionSet and not just ones that overlap.
#' For "multiBase" data, if the "proportionWeightedMean" scoring metric 
#' is used, then the output will also have a "sumProportionOverlap" column.
#' During this scoring method, the proportion overlap between each signalCoord
#' region and overlapping regionSet region is calculated. This column is
#' the sum of all those proportion overlaps and is another method to quantify
#' coverage of regionSet in addition to regionSetCoverage.
#' 
#' @examples
#' data("brcaMCoord1")
#' data("brcaLoadings1")
#' data("esr1_chr1")
#' rsScores <- aggregateSignal(signal=brcaLoadings1, 
#'                                  signalCoord=brcaMCoord1, 
#'                                  regionSet=esr1_chr1, 
#'                                  signalCol=c("PC1", "PC2"), 
#'                                  scoringMetric="regionMean")
#' @export

aggregateSignal <- function(signal,
                              signalCoord,
                              regionSet,
                              signalCol = c("PC1", "PC2"),
                              signalCoordType = "default",
                              scoringMetric = "default",
                              verbose = FALSE,
                              absVal=TRUE) {
    
    ################### checking inputs  #################################
    
    ########## check that inputs are the correct class
    # exports coordinateDT to this environment (converts signalCoord)
    checkConvertInputClasses(signal=signal,
                             signalCoord=signalCoord,
                             regionSet=regionSet,
                             signalCol = signalCol)
    
    ########## check that dimensions of inputs are consistent
    # length of signal coord = nrow of signal
    if (nrow(coordinateDT) != nrow(signal)) {
        stop(cleanws("The number of coordinates in 
            signalCoord (length(signalCoord)) does not equal the number of 
                     rows in `signal`"))
    } 
    
    ######### check that appropriate columns are present
    # signalCol are column names of signal
    if (!all(signalCol %in% colnames(signal))) {
        missingCols = signalCol[!(signalCol %in% colnames(signal))]
        stop(cleanws(paste0("Some signalCol are not 
                            columns of signal: ", missingCols)))
    }
    
    ######## check that scoringMetric is appropriate
    
    if (!(scoringMetric %in% c("default", "regionMean", "simpleMean", 
                               # "meanDiff", "rankSum", 
                               "proportionWeightedMean"))) {
        stop(cleanws("scoringMetric was not recognized. 
                      Check spelling and available options."))
    }
    
    ###### check that signalCoordType is appropriate
    if (!(signalCoordType %in% c("default", "singleBase", "multiBase"))) {
        stop(cleanws("signalCoordType not recognized. 
                         Check spelling/capitalization."))
    }
    
    #######
    # what happens if there are NAs or Inf in `signal`?
    
    #################################################################
    
    # detect signalCoordType
    if (signalCoordType == "default") {
        if (is(signalCoord, "data.frame")) {
            if ("end" %in% colnames(signalCoord)) {
                if (all(signalCoord$start == signalCoord$end)) {
                    signalCoordType <- "singleBase"
                } else {
                    signalCoordType <- "multiBase"
                }
            } else {
                signalCoordType <- "singleBase"
            }
        } else {
            # when signalCoord is a GRanges object
            if (all(start(signalCoord) == end(signalCoord))) {
                signalCoordType <- "singleBase" 
            } else {
                signalCoordType <- "multiBase"
            }
        }
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
    
    # make sure that scoringMetric is consistent with signalCoordType
    if (signalCoordType == "singleBase") {
        if (!(scoringMetric %in% c("regionMean", "simpleMean" #, 
                                   # "rankSum", "meanDiff"
                                   ))) {
            stop("The scoringMetric you selected is not available for
                 this data's signalCoordType")
        }
    } else if (signalCoordType == "multiBase") {
        if (!(scoringMetric %in% c("proportionWeightedMean", "simpleMean"))) {
            stop("The scoringMetric you selected is not available for
                 this data's signalCoordType")
        }
    }
    
    numOfRegions <- length(regionSet)
    totalCpGs    <- nrow(signal)
    
    # extreme positive or negative values both give important information
    # take absolute value or not
    if (absVal) {
        signal <- abs(signal) # required for later code
    }
    
    loadingDT  <- as.data.table(signal)
    
    # reformat into data.table with chromosome location and weight
    # also restricting to signalCol so unnecessary computations 
    # are not done
    loadingDT <- data.table(coordinateDT, 
                            loadingDT[, signalCol, with=FALSE])
    # naming does not work if only using one PC so add this line for that case
    setnames(loadingDT, c(colnames(coordinateDT), signalCol))
    
    # would rounding speed up aggregation?, potentially make a sparse matrix
    # if a lot of entries became 0
    
    # specify aggregation operation
    # will be done separately for each PC specified
    aggrCommand <- buildJ(signalCol, 
                                rep("mean", length(signalCol)))
    

    
    if (signalCoordType == "singleBase") {
        # do the actual aggregation
        if (scoringMetric == "regionMean") {
            
            
                # previously used BSAggregate from RGenomeUtils but now using local, 
                # modified copy
                loadAgMain <- BSAggregate(BSDT = loadingDT, 
                                          regionsGRL = GRangesList(regionSet),
                                          jExpr = aggrCommand,
                                          byRegionGroup = TRUE,
                                          splitFactor = NULL,
                                          returnOLInfo = TRUE)
            
            # if no cytosines from loadings were included in regionSet, result is NA
            if (is.null(loadAgMain)) {
                results <- as.data.table(t(rep(NA, length(signalCol))))
                setnames(results, signalCol)
                results[, signalCoverage := 0]
                results[, regionSetCoverage := 0]
                results[, totalRegionNumber := numOfRegions]
                results[, meanRegionSize := round(mean(width(regionSet)), 1)]
            } else {
                results <- loadAgMain[, .SD, .SDcols = signalCol]
                results[, signalCoverage := loadAgMain[, .SD, .SDcols = "signalCoverage"]]
                results[, regionSetCoverage := loadAgMain[, .SD, .SDcols = "regionSetCoverage"]]
                results[, totalRegionNumber := numOfRegions]
                results[, meanRegionSize := round(mean(width(regionSet)), 1)]
            }
        } else if (scoringMetric == "simpleMean") {
            # average of loadings for all CpGs within region set
            loadMetrics <- signalOLMetrics(dataDT=loadingDT, regionSet=regionSet, 
                                           signalCol = signalCol,
                                           metrics="mean", 
                                           alsoNonOLMet=FALSE)
            if (is.null(loadMetrics)) {
                results <- as.data.table(t(rep(NA, length(signalCol))))
                setnames(results, signalCol)
                results[, signalCoverage := 0]
                results[, regionSetCoverage := 0]
                results[, totalRegionNumber := numOfRegions]
                results[, meanRegionSize := round(mean(width(regionSet)), 1)]
            } else {
                # simple mean 
                results <- as.data.table(t(loadMetrics$mean_OL))
                colnames(results) <- loadMetrics$testCol
                
                # add information about degree of overlap
                results <- cbind(results, 
                                 loadMetrics[1, .SD, 
                                             .SDcols = c("signalCoverage", 
                                                         "regionSetCoverage", 
                                                         "totalRegionNumber", 
                                                         "meanRegionSize")]) 
            }
        # } else if (scoringMetric == "meanDiff") {
        #     # if (is.null(pcLoadAv)) {
        #     #     # calculate (should already be absolute)
        #     #     pcLoadAv <- apply(X = loadingDT[, signalCol, with=FALSE], 
        #     #                       MARGIN = 2, FUN = mean)
        #     # }
        #     loadMetrics <- signalOLMetrics(dataDT=loadingDT, regionSet=regionSet,
        #                                    signalCol = signalCol,
        #                                    metrics=c("mean", "sd"), 
        #                                    alsoNonOLMet=TRUE)
        #     if (is.null(loadMetrics)) {
        #         results <- as.data.table(t(rep(NA, length(signalCol))))
        #         setnames(results, signalCol)
        #         results[, signalCoverage := 0]
        #         results[, regionSetCoverage := 0]
        #         results[, totalRegionNumber := numOfRegions]
        #         results[, meanRegionSize := round(mean(width(regionSet)), 1)]
        #     } else {
        #         # calculate mean difference
        #         # pooled standard deviation
        #         sdPool <- sqrt((loadMetrics$sd_OL^2 + loadMetrics$sd_nonOL^2) / 2)
        #         
        #         # mean difference
        #         # error if signalCoverage > (1/2) * totalCpGs
        #         meanDiff <- (loadMetrics$mean_OL - loadMetrics$mean_nonOL) / 
        #             (sdPool * sqrt((1 / loadMetrics$signalCoverage) - (1 / (totalCpGs - loadMetrics$signalCoverage))))
        #         results <- as.data.table(t(meanDiff))
        #         colnames(results) <- loadMetrics$testCol
        #         
        #         # add information about degree of overlap
        #         results <- cbind(results, loadMetrics[1, .SD, .SDcols = c("signalCoverage", "regionSetCoverage", "totalRegionNumber", "meanRegionSize")]) 
        #     }
        #     
        # } else if (scoringMetric == "rankSum") {
        #     
        #     # if (wilcox.conf.int) {
        #     #     # returns confidence interval
        #     #     wRes <- rsWilcox(dataDT = loadingDT, regionSet=regionSet, 
        #     #                      signalCol = signalCol, 
        #     #                      conf.int = wilcox.conf.int)
        #     #     
        #     #     if (is.null(wRes)) {
        #     #         results <- as.data.table(t(rep(NA, length(signalCol) * 2)))
        #     #         setnames(results, paste0(rep(signalCol, each=2), 
        #     #                                  c("_low", "_high")))
        #     #         results[, signalCoverage := 0]
        #     #         results[, regionSetCoverage := 0]
        #     #         results[, totalRegionNumber := numOfRegions]
        #     #         results[, meanRegionSize := round(mean(width(regionSet)), 1)]
        #     #     } else {
        #     #         results <- as.data.table(wRes)
        #     #     }    
        #     #     
        #     # } else {
        #         # returns p value
        #         # one sided test since I took the absolute value of the loadings 
        #         wRes <- rsWilcox(dataDT = loadingDT, regionSet=regionSet, 
        #                          signalCol = signalCol, alternative="greater")
        #         
        #         if (is.null(wRes)) {
        #             results <- as.data.table(t(rep(NA, length(signalCol))))
        #             setnames(results, signalCol)
        #             results[, signalCoverage := 0]
        #             results[, regionSetCoverage := 0]
        #             results[, totalRegionNumber := numOfRegions]
        #             results[, meanRegionSize := round(mean(width(regionSet)), 1)]
        #         } else {
        #             results <- as.data.table(wRes)
        #         }    
        #     #}
        }  
    } else {

        # signalCoordType == "multiBase"
        # for ATAC-seq
        if (scoringMetric == "proportionWeightedMean") {
            loadAgMain <- regionOLWeightedMean(signalDT = loadingDT, 
                                               signalGR = dtToGr(coordinateDT),
                                               regionSet = regionSet,
                                               calcCols= signalCol)
            
            # setnames(loadAgMain, c("signalCoverage", "regionSetCoverage"),
            #          c("numCpGsOverlapping", "numRegionsOverlapping"))
        } else if (scoringMetric == "simpleMean") {
            loadAgMain <- regionOLMean(signalDT = loadingDT, 
                                       signalGR = dtToGr(coordinateDT),
                                       regionSet = regionSet,
                                       calcCols= signalCol)
            # setnames(loadAgMain, c("signalCoverage", "regionSetCoverage"),
            #          c("numCpGsOverlapping", "numRegionsOverlapping"))
        }

        
        # in the case of no overlaps
        if (is.null(loadAgMain)) {
            results <- as.data.table(t(rep(NA, length(signalCol))))
            setnames(results, signalCol)
            results[, signalCoverage := 0]
            results[, regionSetCoverage := 0]
            results[, totalRegionNumber := numOfRegions]
            results[, meanRegionSize := round(mean(width(regionSet)), 1)]
            # this column is only added by this scoring method
            if (scoringMetric == "proportionWeightedMean") {
                results[, sumProportionOverlap := 0]
            }
        } else {
            loadAgMain <- as.data.table(loadAgMain)
            results <- loadAgMain[, .SD, .SDcols = signalCol]
            results[, signalCoverage := loadAgMain[, .SD, .SDcols = "signalCoverage"]]
            results[, regionSetCoverage := loadAgMain[, .SD, .SDcols = "regionSetCoverage"]]
            results[, totalRegionNumber := numOfRegions]
            results[, meanRegionSize := round(mean(width(regionSet)), 1)]
            # this column is only added by this scoring method
            if (scoringMetric == "proportionWeightedMean") {
                results[, sumProportionOverlap := loadAgMain[, .SD, .SDcols = "sumProportionOverlap"]]
            }
        }
    }
    
    # signalOLMetrics() # make sure it works with no overlap
    if (verbose) {
        message("|")
    }
        
    return(as.data.frame(results))
}



#' Do COCOA with many region sets
#' 
#' This function will give each region set a score for each PC
#' in `signalCol` based on
#' the `scoringMetric` parameter. Based on these scores, you can determine
#' which region sets out of a region set database (given by GRList) 
#' are most associated with the top PCs. See the vignette "Introduction
#' to Coordinate Covariation Analysis" for help interpreting your 
#' results. 
#'
#' @param signal matrix of loadings (the coefficients of 
#' the linear combination that defines each PC). One named column for each PC.
#' One row for each original dimension/variable (should be same order 
#' as original data/signalCoord). The x$rotation output of prcomp().
#' @template signalCoord
#' @template signalCoordType
#' @template GRList
#' @template signalCol
#' @param scoringMetric A character object with the scoring metric.
#' There are different scoring metrics available for 
#' signalCoordType="singleBase" vs  signalCoordType="multiBase".
#' For "singleBase", the available scoring methods are "regionMean", 
#' "simpleMean", and "rankSum". The default method is "regionMean".
#' For "multiBase", the scoring methods are "proportionWeightedMean" and 
#' "simpleMean". The default is "proportionWeightedMean".
#' "regionMean" is a weighted
#' average of the signal, weighted by region (absolute value of signal 
#' if absVal=TRUE). First the signal is
#' averaged within each regionSet region, 
#' then all the regions are averaged. With
#' "regionMean" score, be cautious in interpretation for
#' region sets with low number of regions that overlap signalCoord. 
#' The "simpleMean"
#' method is just the unweighted average of all (absolute) signal values that
#' overlap the given region set. For multiBase data, this includes
#' signal regions that overlap a regionSet region at all (1 base
#' overlap or more) and the signal for each overlapping region is
#' given the same weight for the average regardless of how much it overlaps. 
#' "proportionWeightedMean" is a weighted average of all signalCoord 
#' regions that overlap with regionSet regions. For each signalCoord region
#' that overlaps with a regionSet region, we calculate what proportion
#' of the regionSet region is covered. Then this proportion is used to
#' weight the signal value when calculating the mean. 
#' The denominator of the mean
#' is the sum of all the proportion overlaps. 
#' Wilcoxon rank sum test ("rankSum") is also supported but is
#' skewed toward ranking large region sets highly and is
#' significantly slower than the "regionMean" method. 
#' For the ranksum method, the absolute loadings for loadings that
#' overlap the given region set are taken as a group and all the
#' loadings that do not overlap the region set are taken as
#' the other group. Then p value is then given as the score.
#' It is a one sided test, with the alternative hypothesis
#' that the loadings in the region set will be greater than
#' the loadings not in the region set.
#' @template verbose
# @param wilcox.conf.int logical. Only applies when using "rankSum" scoring
# method. returns a 95% confidence interval from the Wilcoxon rank sum test
# instead of p value.
#' @template absVal
#' @return data.frame of results, one row for each region set. 
#' It has the following columns:
#' one column for each item of signalCol with names given
#' by signalCol. These columns have scores for the region set for each PC.
#' Other columns: signalCoverage (formerly cytosine_coverage) which
#' has number of cytosines that overlapped with regionSet (or for "multiBase"
#' data, 
#' the number of regions from signalCoord that overlapped regionSet) 
#' regionSetCoverage which has number of regions from regionSet
#' that overlapped any coordinates from signalCoord, 
#' totalRegionNumber that has
#' number of regions in regionSet, meanRegionSize that has average
#' size in base pairs of regions in regionSet, the average is based on
#' all regions in regionSet and not just ones that overlap.
#' For "multiBase" data, if the "proportionWeightedMean" scoring metric 
#' is used, then the output will also have a "sumProportionOverlap" column.
#' During this scoring method, the proportion overlap between each signalCoord
#' region and overlapping regionSet region is calculated. This column is
#' the sum of all those proportion overlaps and is another method to quantify
#' coverage of regionSet in addition to regionSetCoverage.
#' 
#' 
#' @examples 
#' data("brcaMCoord1")
#' data("brcaLoadings1")
#' data("esr1_chr1")
#' rsScores <- runCOCOA(signal=brcaLoadings1, 
#'                                  signalCoord=brcaMCoord1, 
#'                                  GRList=GRangesList(esr1_chr1), 
#'                                  signalCol=c("PC1", "PC2"), 
#'                                  scoringMetric="regionMean")
#' 
#' @export

runCOCOA <- function(signal,
                     signalCoord,
                     GRList,
                     signalCol = c("PC1", "PC2"),
                     signalCoordType = "default",
                     scoringMetric = "default",
                     verbose = TRUE,
                     absVal=TRUE) {

    ################### checking inputs  #################################
    
    ########## check that inputs are the correct class
    # exports coordinateDT to this environment (converts signalCoord)
    checkConvertInputClasses(signal=signal,
                             signalCoord=signalCoord,
                             regionSet=NULL,
                             signalCol = signalCol,
                             GRList=GRList)
    
    ########## check that dimensions of inputs are consistent
    # length of signal coord = nrow of signal
    if (nrow(coordinateDT) != nrow(signal)) {
        stop(cleanws("The number of coordinates in 
            signalCoord (length(signalCoord)) does not equal the number of 
                     rows in `signal`"))
    } 
    
    ######### check that appropriate columns are present
    # signalCol are column names of signal
    if (!all(signalCol %in% colnames(signal))) {
        missingCols = signalCol[!(signalCol %in% colnames(signal))]
        stop(cleanws(paste0("Some signalCol are not 
                            columns of signal: ", missingCols)))
    }
    
    ######## check that scoringMetric is appropriate
    
    if (!(scoringMetric %in% c("default", "regionMean", "simpleMean", 
                               # "meanDiff", "rankSum", 
                               "proportionWeightedMean"))) {
        stop(cleanws("scoringMetric was not recognized. 
                      Check spelling and available options."))
    }
    
    ###### check that signalCoordType is appropriate
    if (!(signalCoordType %in% c("default", "singleBase", "multiBase"))) {
        stop(cleanws("signalCoordType not recognized. 
                     Check spelling/capitalization."))
    }
    
    #######
    # what happens if there are NAs or Inf in `signal`?
    
    #################################################################
    
    # detect signalCoordType
    if (signalCoordType == "default") {
        if (is(signalCoord, "data.frame")) {
            if ("end" %in% colnames(signalCoord)) {
                if (all(signalCoord$start == signalCoord$end)) {
                    signalCoordType <- "singleBase"
                } else {
                    signalCoordType <- "multiBase"
                }
            } else {
                signalCoordType <- "singleBase"
            }
        } else {
            # when signalCoord is a GRanges object
            if (all(start(signalCoord) == end(signalCoord))) {
                signalCoordType <- "singleBase" 
            } else {
                signalCoordType <- "multiBase"
            }
        }
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
    
    
    # apply over the list of region sets
    resultsList <- lapplyAlias(GRList,
                               function(x) aggregateSignal(
                                signal = signal,
                                signalCoord = coordinateDT,
                                signalCoordType = signalCoordType,
                                regionSet = x,
                                signalCol = signalCol,
                                scoringMetric = scoringMetric,
                                verbose = verbose,
                                absVal = absVal))
                                #wilcox.conf.int = wilcox.conf.int,
    resultsDT <- do.call(rbind, resultsList) 

    # # add names if they are present
    # if (!is.null(names(GRList))) {
    #     # resultsDT[, rsName := names(GRList)]
    #     row.names(resultsDT) <- names(GRList)
    # }

    # rsName <- row.names(resultsDT)
    resultsDF <- as.data.frame(resultsDT) #, row.names = rsName)

    return(resultsDF)
}    

# @param dataMat columns of dataMat should be samples/patients, rows should be genomic signal
# (each row corresponds to one genomic coordinate/range)
# @param featureMat Rows should be samples, columns should be "features" 
# (whatever you want to get correlation with: eg PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns)
# @param centerDataMat logical object. Should rows in dataMat be centered based on
# their means? (subtracting row means from each row)
# @param centerFeatureMat. logical.
# @param testType character object. Can be "cor" (Pearson correlation),
# "spearmanCor (Spearman correlation)
# "pcor" (partial correlation), "cov" (covariance (Pearson)), "spearmanCov" 
# (Spearman covariance)
# @param covariate
#
# If a row in dataMat has 0 stand. deviation, correlation will be set to 0
# instead of NA as would be done by cor()
#
# @return returns a matrix where rows are the genomic signal (eg a CpG or region) and
# columns are the columns of featureMat
# @examples dataMat = matrix(rnorm(50), 5, 10)
# featureMat = matrix(rnorm(20), 10, 2)
createCorFeatureMat <- function(dataMat, featureMat, 
                               centerDataMat=TRUE, centerFeatureMat = TRUE, 
                               testType="cor", covariate=NULL) {
    
    featureMat <- as.matrix(featureMat)
    featureNames <- colnames(featureMat)
    nFeatures <- ncol(featureMat)
    nDataDims <- nrow(dataMat)
    
    if (centerDataMat) {
        cpgMeans <- rowMeans(dataMat, na.rm = TRUE)
        # centering before calculating correlation
        dataMat <- apply(X = dataMat, MARGIN = 2, function(x) x - cpgMeans)
        
    }
    
    if (centerFeatureMat) {
        featureMeans <- colMeans(featureMat, na.rm = TRUE)
        # centering before calculating correlation
        featureMat <- t(apply(X = t(featureMat), MARGIN = 2, function(x) x - featureMeans))
        if (dim(featureMat)[1] == 1) {
            featureMat <- t(featureMat)
        }
        featureMat <- as.matrix(featureMat)
    }
    
    # avoid this copy and/or delay transpose until after calculating correlation?
    dataMat <- as.data.frame(t(dataMat))
    
    
    if (testType == "cor") {
        # create feature correlation matrix with PCs (rows: features/CpGs, columns:PCs)
        # how much do features correlate with each PC?
        
        # put epigenetic data first in cor()
        featurePCCor <- cor(dataMat, featureMat, use="pairwise.complete.obs", method="pearson")

    } else if (testType == "spearmanCor") {
        # xtfrm(x) ranking
        featurePCCor <- cor(dataMat, featureMat, use="pairwise.complete.obs", method="spearman")

    } else if (testType == "pcor") {
        # partial correlation (account for covariates), ppcor package
        
        featurePCCor <- apply(X = featureMat, MARGIN = 2, function(y) apply(X = dataMat, 2, 
                                                                           FUN = function(x) pcor.test(x = x, y=y,
                                                                                                       z=covariate,
                                                                                                       method="pearson")$estimate))
        
    } else if (testType == "cov") {
        featurePCCor <- cov(dataMat, featureMat, use="pairwise.complete.obs")

    } else {
        stop("invalid testType")
    }
    
    
    # if standard deviation of the data was zero, NA will be produced
    # set to 0 because no standard deviation means no correlation with attribute of interest
    featurePCCor[is.na(featurePCCor)] <- 0
    colnames(featurePCCor) <- featureNames
    
    return(featurePCCor)
    # corLoadRatio <- signal[, signalCol] / featurePCCor 
    # hist(corLoadRatio[, "PC10"])
}


#' Create a "meta-region" loading profile 
#' 
#' This loading profile can show enrichment 
#' of genomic signals with high loading values in region set but not in
#' surrounding genome, suggesting that variation is linked specifically
#' to that region set. 
#' 
#' All regions in a given region set 
#' are combined into a single aggregate profile. Regions should be
#' expanded on each side to include a wider area of the genome around
#' the regions of interest (see example and vignettes). 
#' To make the profile, first we take 
#' the absolute value of the loadings. Then each region is
#' split into `binNum` bins. All loadings in each bin are 
#' averaged to get one value per bin. Finally, corresponding bins from
#' the different regions are averaged (eg all bin1's averaged with each other, 
#' all bin2's averaged with each other, etc.) to get a single "meta-region"
#' loading profile. Since DNA strand information is not considered, 
#' the profile is averaged symmetrically around the center.
#' A peak in the middle of this profile suggests
#' that variability is specific to the region set of interest and is 
#' not a product of the surrounding genome. A region set can still be
#' significant even if it does not have a peak. For example, some
#' histone modification region sets may be in large genomic blocks
#' and not show a peak, despite having variation across samples.
#'
#'
#' @param signal matrix of loadings (the coefficients of 
#' the linear combination that defines each PC). One named column for each PC.
#' One row for each original dimension/variable (should be same order 
#' as original data/signalCoord). Given by prcomp(x)$rotation.
#' @template signalCoord
#' @template regionSet
#' @template signalCol
#' @template signalCoordType
#' @param binNum Number of bins to split each region into when
#' making the aggregate loading profile. More bins will
#' give a higher resolution but perhaps more noisy profile.
#' @param verbose A "logical" object. Whether progress 
#' of the function should be shown, one
#' bar indicates the region set is completed. Useful when
#' using `lapply` to get the loading profiles of many region sets.
#' @param aggrMethod character. A character object with the aggregation method.
#' Similar to `scoringMetric`.
#' There are different aggregation methods available for 
#' signalCoordType="singleBase" vs  signalCoordType="multiBase".
#' For "singleBase", the available scoring methods are "regionMean" and
#' "simpleMean". The default method is "regionMean".
#' For "multiBase", the scoring methods are "proportionWeightedMean" and 
#' "simpleMean". The default is "proportionWeightedMean".
#' "regionMean" is a weighted
#' average of the signal, weighted by region (absolute value of signal 
#' if absVal=TRUE). First the signal is
#' averaged within each regionSet region, 
#' then all the regions are averaged. With
#' "regionMean" score, be cautious in interpretation for
#' region sets with low number of regions that overlap signalCoord. 
#' The "simpleMean"
#' method is just the unweighted average of all (absolute) signal values that
#' overlap the given region set. For multiBase data, this includes
#' signal regions that overlap a regionSet region at all (1 base
#' overlap or more) and the signal for each overlapping region is
#' given the same weight for the average regardless of how much it overlaps. 
#' "proportionWeightedMean" is a weighted average of all signalCoord 
#' regions that overlap with regionSet regions. For each signalCoord region
#' that overlaps with a regionSet region, we calculate what proportion
#' of the regionSet region is covered. Then this proportion is used to
#' weight the signal value when calculating the mean. 
#' The denominator of the mean
#' is the sum of all the proportion overlaps. 
#' @template absVal
#' @return A data.frame with the binned loading profile,
#' one row per bin. columns: binID and one column for each PC
#' in signalCol. The function will return NULL if there
#' is no overlap between regionSet and signalCoord.
#' 
#' @examples 
#' data("brcaMCoord1")
#' data("brcaLoadings1")
#' data("esr1_chr1")
#' esr1_chr1_expanded <- resize(esr1_chr1, 14000, fix="center")
#' getMetaRegionProfile(signal=brcaLoadings1, 
#'                     signalCoord=brcaMCoord1, 
#'                     regionSet=esr1_chr1_expanded, 
#'                     signalCol=c("PC1", "PC2"), 
#'                     binNum=25)
#' @export

getMetaRegionProfile <- function(signal, signalCoord, regionSet,
                    signalCol = c("PC1", "PC2"),
                    signalCoordType = "default", 
                    binNum = 25,
                    verbose=TRUE,  
                    aggrMethod = "default", absVal=TRUE) {
    
    ################### checking inputs  #################################
    
    ########## check that inputs are the correct class
    # exports coordinateDT to this environment (converts signalCoord)
    checkConvertInputClasses(signal=signal,
                             signalCoord=signalCoord,
                             regionSet=regionSet,
                             signalCol = signalCol)
    
    ########## check that dimensions of inputs are consistent
    # length of signal coord = nrow of signal
    if (nrow(coordinateDT) != nrow(signal)) {
        stop(cleanws("The number of coordinates in 
            signalCoord (length(signalCoord)) does not equal the number of 
                     rows in `signal`"))
    } 
    
    ######### check that appropriate columns are present
    # signalCol are column names of `signal`
    if (!all(signalCol %in% colnames(signal))) {
        missingCols = signalCol[!(signalCol %in% colnames(signal))]
        stop(cleanws(paste0("Some signalCol are not 
                            columns of signal: ", missingCols)))
    }
    
    
    ######## check that aggregation method is appropriate
    
    if (!(aggrMethod %in% c("default", "regionMean", "simpleMean", 
                               "proportionWeightedMean"))) {
        stop(cleanws("scoringMetric was not recognized. 
                     Check spelling and available options."))
    }
    
    ###### check that signalCoordType is appropriate
    if (!(signalCoordType %in% c("default", "singleBase", "multiBase"))) {
        stop(cleanws("signalCoordType not recognized. 
                     Check spelling/capitalization."))
    }
    
    #######
    # what happens if there are NAs or Inf in `signal`?
    
    #################################################################
    
    # detect signalCoordType
    if (signalCoordType == "default") {
        if (is(signalCoord, "data.frame")) {
            if ("end" %in% colnames(signalCoord)) {
                if (all(signalCoord$start == signalCoord$end)) {
                    signalCoordType <- "singleBase"
                } else {
                    signalCoordType <- "multiBase"
                }
            } else {
                signalCoordType <- "singleBase"
            }
        } else {
            # when signalCoord is a GRanges object
            if (all(start(signalCoord) == end(signalCoord))) {
                signalCoordType <- "singleBase" 
            } else {
                signalCoordType <- "multiBase"
            }
        }
    }
    
    # if "default" aggregation method is given, choose based on signalCoordType
    if (aggrMethod == "default") {
        if (signalCoordType == "singleBase") {
            aggrMethod <- "regionMean"   
        } else if (signalCoordType == "multiBase") {
            aggrMethod <- "proportionWeightedMean"
        } else {
            stop(cleanws("signalCoordType not recognized. 
                         Check spelling/capitalization."))
        }
    }
    
    
    
    
    ##################################################################
    # take absolute value or not
    if (absVal) {
        loadingDT <- as.data.table(abs(signal))
    } else {
        loadingDT <- as.data.table(signal)
    }
    loadingDT <- cbind(coordinateDT, loadingDT)
    
    GRDT <- grToDt(regionSet)
    
    loadProf <- BSBinAggregate(BSDT = loadingDT, 
                               rangeDT = GRDT, 
                               binCount = binNum,
                               byRegionGroup = TRUE, 
                               splitFactor = NULL,
                               signalCol = signalCol,
                               aggrMethod = aggrMethod)
    
    # if loadProf is NULL, return NULL from function, otherwise make symmetrical
    # it will be NULL when there was no overlap between data and any of the bins
    if (!is.null(loadProf)) {
        loadProf <- makeSymmetric(loadProf)
        loadProf[, regionGroupID := seq_len(binNum)][]
        setnames(loadProf, old = "regionGroupID", new="binID")
        loadProf <- as.data.frame(loadProf)
    } else {
        warning("Insufficient overlap between regionSet and signalCoord")
    }
    
    return(loadProf)
}

makeSymmetric <- function(prof) {
    symProf <- apply(prof, 2, .makeSymmetric)
    symProf <- as.data.table(symProf)
    return(symProf)
}

.makeSymmetric <- function(vec) {
    symVec <- (vec + rev(vec)) / 2
    return(symVec)
}

# Produced originally for binning Ewing RRBS data across various region sets
#
# @param BSDT A data.table. For COCOA, a data.table of loading values
# with the PCs to be annotated. One column for the loadings of each PC
# and also has columns with the coordinates for CpGs that the loadings
# are for: chr (chromosome) and start column
# @param rangeDT A data table with the sets of regions to be binned, 
# with columns named start, end
# @param binCount Number of bins across the region
# @param byRegionGroup Pass along to binCount (see ?binCount)
# @template signalCol
# @param verbose A "logical" object. Whether progress 
# of the function should be shown, one
# bar indicates the region set is completed.
# useful when using BSBinAggregate with 'apply' to do many 
# region sets at a time.
# @param aggrMethod 
# "unweightedMean" determines whether single, or multiple, genomic signal(s) 
# that span(s) more than a single position overlap(s) with a region set 
# location. This method determines the unweighted mean signal of any genomic 
# signal that overlaps in any amount with a region set region.
# "proportionWeightedMean" determines whether single, or multiple, 
# genomic signal(s) that span(s) more than a single position overlap(s) with a 
# region set location. This method weights the signals by the proportion 
# overlap with the corresponding region set regions.
BSBinAggregate <- function(BSDT, rangeDT, binCount,
                           byRegionGroup = TRUE,
                           splitFactor = NULL,
                           signalCol,
                           verbose = FALSE,
                           aggrMethod) {
    if (!is(rangeDT, "data.table")) {
    stop("rangeDT must be a data.table")
}
    seqnamesColName <- "seqnames"  # default column name
    if (! "seqnames" %in% colnames(rangeDT)) {
        if ("chr" %in% colnames(rangeDT)) {
            # message("seqnames column name set to: chr")
            seqnamesColName <- "chr"
        } else {
            # Got neither.
            stop("rangeDT must have a seqnames column")
        }
    }
    
    # message("Binning...")
    binnedDT <- rangeDT[, MIRA::binRegion(start, end, 
                                          binCount, get(seqnamesColName))]
    # output is a list of GRanges objects, does not play well with vapply
    # one GRanges object for each bin, containing a segment of each original rangeDT region 
    binnedGR <- sapply(split(binnedDT, binnedDT$binID), dtToGr)
    # message("Aggregating...")
    
    if (aggrMethod == "proportionWeightedMean") {

        binMeansList <- lapply(X = binnedGR, 
                               FUN = function(x) regionOLWeightedMean(signalDT = BSDT, 
                                                                      signalGR = dtToGr(BSDT[, .(chr, start, end)]), 
                                                                      regionSet = x, 
                                                                      calcCols = signalCol))
        
        # any bins that had no overlap with data will be NULL
        # if any bins had no data, return NULL
        if (any(vapply(X = binMeansList, FUN = is.null, FUN.VALUE = TRUE))) {
            return(NULL)
        }
        
        # rbindlist of all NULL's will still return an (empty) data.table
        # otherwise any NULL items will just be skipped and other rows will
        # be concatenated
        binnedBSDT <- rbindlist(binMeansList)
        regionGroupID = 1:length(binMeansList)

        # regionGroupID = regionGroupID[!vapply(X = binMeansList, FUN = is.null, FUN.VALUE = TRUE)]
        binnedBSDT[, regionGroupID := regionGroupID]

    } else if (aggrMethod == "simpleMean") {

        binMeansList <- lapply(X = binnedGR, 
                               FUN = function(x) regionOLMean(signalDT = BSDT, 
                                                              signalGR = dtToGr(BSDT[, .(chr, start, end)]), 
                                                              regionSet = x, 
                                                              calcCols = signalCol))
        # any bins that had no overlap with data will be NULL
        # if any bins had no data, return NULL
        if (any(vapply(X = binMeansList, FUN = is.null, FUN.VALUE = TRUE))) {
            return(NULL)
        }
        
        binnedBSDT <- rbindlist(binMeansList)
        regionGroupID = 1:length(binMeansList)

         
        # regionGroupID = regionGroupID[!vapply(X = binMeansList, FUN = is.null, FUN.VALUE = TRUE)]
        binnedBSDT[, regionGroupID := regionGroupID]

    } else if (aggrMethod == "regionMean") { # aggrMethod == "regionMean"
        
        # what is output if a region set has no overlap?
        binnedBSDT <- BSAggregate(BSDT,
                                  regionsGRL=GRangesList(binnedGR),
                                  jExpr=buildJ(signalCol,
                                               rep("mean", length(signalCol))),
                                  byRegionGroup = byRegionGroup,
                                  splitFactor = splitFactor)
    }
    # RGenomeUtils::BSAggregate

    # # If we aren't aggregating by bin, then don't restrict to min reads!
    # if (byRegionGroup) {
    #     binnedBSDT <- binnedBSDT[readCount > minReads,]
    # }
    if (verbose) {
        message("|")
    }
    
    return(binnedBSDT)
}

# modification of BSAggregate to just return mean per region
# 
# @param signal matrix of loadings (the coefficients of 
# the linear combination that defines each PC). One named column for each PC.
# One row for each original dimension/variable (should be same order 
# as original data/signalCoord). The x$rotation output of prcomp().
# @template signalCoord
# @template regionSet
# Must be from the same reference genome
# as the coordinates for the actual data/samples (signalCoord).
# @template signalCol
# @param returnQuantile "logical" object. If FALSE, return region averages. If TRUE,
# for each region, return the quantile of that region's average value
# based on the distribution of individual genomic signal/feature values
# @template absVal
# @return a data.table with region coordinates and average loading 
# values for each region. Has columns chr, start, end, and a column for each
# PC in signalCol. Regions are not in order along the rows of the data.table.
#
# @example averagePerRegion(BSDT = BSDT, regionsGRL, 
#          jCommand = MIRA:::buildJ(cols = "methylProp", "mean")) 
# Devel note: I could add a column for how many cytosines are in each region 

averagePerRegion <- function(signal,
                             signalCoord,
                             regionSet,
                             signalCol = c("PC1", "PC2"),
                             returnQuantile = FALSE,
                             absVal=TRUE) {

    ################### checking inputs  #################################
    
    ########## check that inputs are the correct class
    # exports coordinateDT to this environment (converts signalCoord)
    checkConvertInputClasses(signal=signal,
                             signalCoord=signalCoord,
                             regionSet=regionSet,
                             signalCol = signalCol)
    
    ########## check that dimensions of inputs are consistent
    # length of signal coord = nrow of signal
    if (nrow(coordinateDT) != nrow(signal)) {
        stop(cleanws("The number of coordinates in 
            signalCoord (length(signalCoord)) does not equal the number of 
                     rows in `signal`"))
    } 
    
    ######### check that appropriate columns are present
    # signalCol are column names of signal
    if (!all(signalCol %in% colnames(signal))) {
        missingCols = signalCol[!(signalCol %in% colnames(signal))]
        stop(cleanws(paste0("Some signalCol are not 
                            columns of signal: ", missingCols)))
    }
    
    #######
    # what happens if there are NAs or Inf in `signal`?
    
    #################################################################

    # determine whether coordinates are single base or a range
    if (!("end" %in% colnames(coordinateDT))) {
        dataCoordType <- "singleBase"
    } else {
        if (all(coordinateDT$start == coordinateDT$end)) {
            dataCoordType <- "singleBase"
        } else {
            dataCoordType <- "multiBase"
        }
    }

    # take absolute value or not
    if (absVal) {
        signalDT <- as.data.table(abs(signal))
    } else {
        signalDT <- as.data.table(signal)
    }

    # use different function for single base data and for region data
    if (dataCoordType == "singleBase") {
        # linking coordinates to loading values, has columns chr start, signalCol
        BSDT  <- cbind(coordinateDT, signalDT[, .SD, .SDcols = signalCol])
        jExpr <- buildJ(signalCol, rep("mean", length(signalCol)))

        avPerRegion <- BSAggregate(BSDT = BSDT,
                                   regionsGRL = regionSet,
                                   excludeGR = NULL,
                                   regionsGRL.length = NULL,
                                   splitFactor = NULL, 
                                   keepCols = NULL,
                                   sumCols = NULL,
                                   jExpr = jExpr,
                                   byRegionGroup = FALSE,
                                   keep.na = FALSE,
                                   returnSD = FALSE,
                                   returnOLInfo = FALSE,
                                   meanPerRegion = TRUE,
                                   returnQuantile = returnQuantile) 

    } else if (dataCoordType == "multiBase") {

        avPerRegion <- weightedAvePerRegion(signalDT = signalDT,
                             signalCoord=signalCoord,
                             regionSet=regionSet,
                             calcCols = signalCol,
                             returnQuantile = returnQuantile) 

    } else {
        stop("dataCoordType not recognized.")
    }

    return(avPerRegion)

}

# signalDT 
# average by region for region-based data (eg ATAC-seq)
# @return data.table. chr, start, end columns. One column for each calcCols
# that has the average value in each region for that col
weightedAvePerRegion <- function(signalDT,
                                 signalCoord,
                                 regionSet,
                                 calcCols = c("PC1", "PC2"),
                                 returnQuantile = FALSE) {
    
    if (is(signalCoord, "data.frame")) {
        signalCoord <- dtToGr(signalCoord)
    }

    hits  <- findOverlaps(query = signalCoord, subject = regionSet)
    # if no overlap, return NULL
    if (length(hits) == 0) {
        return(NULL)
    }

    olap  <- pintersect(signalCoord[queryHits(hits)],
                        regionSet[subjectHits(hits)])
    polap <- width(olap) / width(regionSet[subjectHits(hits)])

    # get total proportion overlap per region
    # aggregate polap by region
    pOlapDT <- data.table(signalDT[queryHits(hits), calcCols, with=FALSE],
                          rsRegionID = subjectHits(hits),
                          pOlap = polap)
    pOlapByRegionDT <- pOlapDT[, .(regionPOlap = sum(pOlap)), by=rsRegionID]

    # specify aggregation operation
    # will be done separately for each PC specified
    aggrCommand <- paste("list(", paste(paste0(calcCols, "=", "sum", "(",
                                calcCols, " * pOlap)"), collapse = ", "), ")")
    weightedSumByRegionDT <- pOlapDT[, eval(parse(text=aggrCommand)), by=rsRegionID]

    # weightedSumByRegionDT and pOlapByRegionDT should be in the same order 
    regionInd <- weightedSumByRegionDT$rsRegionID
    olCoord <- grToDt(regionSet)[regionInd, .(chr, start, end)]
    
    # divide by total proportion overlap to get mean value 
    jCommand <- paste("list(", paste(paste0(calcCols, "=", calcCols, " / regionPOlap"), collapse = ", "), ")")
    meanPerRegion <- cbind(pOlapByRegionDT, weightedSumByRegionDT)
    meanPerRegion <- meanPerRegion[, eval(parse(text=jCommand))]
    meanPerRegion <- cbind(olCoord, meanPerRegion)
    
    if (returnQuantile) {
        for (i in seq_along(calcCols)) {
            # perhaps this could be more efficient with mapply
            # ecdf example: ecdf(distributionData)(getPercentileOfThis)
            meanPerRegion[, c(calcCols[i]) := ecdf(signalDT[[calcCols[i]]])(meanPerRegion[[calcCols[i]]])]
            # I tried set() to improve performance but it took about the same time
        }
    }

    # meanPerRegion <-  pOlapDT[, .(regionMean = sum(score * (pOlap/sum(pOlap))), by=rsRegionID]

    return(meanPerRegion)
}




#' Get regions that are most associated with PCs of interest
#'
#' Get a GRanges with top regions from the region set based on average
#' loadings for the regions or the quantile of the region's loading.
#' Returns average loading or quantile as GRanges metadata.
#' 
#' @param signal matrix of loadings (the coefficients of 
#' the linear combination that defines each PC). One named column for each PC.
#' One row for each original dimension/variable (should be same order 
#' as original data/signalCoord). The x$rotation output of prcomp().
#' @template signalCoord
#' @template regionSet
#' @template signalCol
#' @param returnQuantile "logical" object. If FALSE, return region averages. If TRUE,
#' for each region, return the quantile of that region's average value
#' based on the distribution of individual genomic signal/feature values
#' @return a GRanges object with region coordinates for regions with
#' scores/quantiles above "cutoff" for any PC in signalCol. The scores/quantiles
#' for signalCol are given as metadata in the GRanges.

# Are regions in order along the rows of the data.table?
#
#' @examples 
#' data("brcaLoadings1")
#' data("brcaMCoord1")
#' data("esr1_chr1")
#' getTopRegions(signal=brcaLoadings1,
#' signalCoord=brcaMCoord1, regionSet=esr1_chr1, returnQuantile = TRUE)

getTopRegions <- function(signal, 
                          signalCoord, 
                          regionSet, 
                          signalCol = c("PC1", "PC2"), cutoff = 0.8, 
                          returnQuantile=TRUE) {
    
    regionLoadDT <- averagePerRegion(signal=signal,
                            signalCoord=signalCoord, regionSet=regionSet, 
                            signalCol = signalCol,
                            returnQuantile = returnQuantile)[]
    
    keepInd <- regionLoadDT[, signalCol, with=FALSE] >= cutoff
    
    # keep region if it is above cutoff in any of the PCs in signalCol
    keepInd <- apply(X = keepInd, MARGIN = 1, FUN = any)
    
    highGR <- dtToGr(regionLoadDT[keepInd, ])
    
    values(highGR) <- as.data.frame(regionLoadDT[keepInd, signalCol, with=FALSE])
    
    return(highGR)
}


# Get regions that are most associated with PCs of interest
#
# Get a GRanges with top regions from the region set based on average
# loadings for the regions or the quantile of the region's loading.
# Returns average loading or quantile as GRanges metadata.
# 
# @param loadingMat matrix of loadings (the coefficients of 
# the linear combination that defines each PC). One named column for each PC.
# One row for each original dimension/variable (should be same order 
# as original data/signalCoord). The x$rotation output of prcomp().
# @param signalCoord a GRanges object or data frame with coordinates 
# for the genomic signal/original data (eg DNA methylation) 
# included in the PCA. Coordinates should be in the 
# same order as the original data and the loadings 
# (each item/row in signalCoord
# corresponds to a row in loadingMat). If a data.frame, 
# must have chr and start columns. If end is included, start 
# and end should be the same. Start coordinate will be used for calculations.
# @param regionSet A GRanges object with regions corresponding
# to the same biological annotation.
# @param PCsToAnnotate A character vector with principal components to  
# include. eg c("PC1", "PC2") These should be column names of loadingMat.
# @param returnQuantile "logical" object. If FALSE, return region averages. If TRUE,
# for each region, return the quantile of that region's average value
# based on the distribution of individual genomic signal/feature values
# @return a GRanges object with region coordinates for regions with
# scores/quantiles above "cutoff" for any PC in PCsToAnnotate. The scores/quantiles
# for PCsToAnnotate are given as metadata in the GRanges.

# Are regions in order along the rows of the data.table?
#
# @examples data("brcaLoadings1")
# data("brcaMCoord1")
# data("esr1_chr1")
# COCOA:::getTopRegions(loadingMat=brcaLoadings1,
# signalCoord=brcaMCoord1, regionSet=esr1_chr1, returnQuantile = TRUE)

getTopRegions <- function(loadingMat, 
                          signalCoord, 
                          regionSet, 
                          PCsToAnnotate = c("PC1", "PC2"), cutoff = 0.8, 
                          returnQuantile=TRUE) {
    
    
    regionLoadDT = COCOA:::averageByRegion(loadingMat=loadingMat,
                            signalCoord=signalCoord, regionSet=regionSet, 
                            PCsToAnnotate = PCsToAnnotate,
                            returnQuantile = returnQuantile)[]
    
    keepInd = regionLoadDT[, PCsToAnnotate, with=FALSE] >= cutoff
    
    # keep region if it is above cutoff in any of the PCs in PCsToAnnotate
    keepInd = apply(X = keepInd, MARGIN = 1, FUN = any)
    
    highGR = COCOA:::dtToGr(regionLoadDT[keepInd, ])
    
    values(highGR) <- as.data.frame(regionLoadDT[keepInd, PCsToAnnotate, with=FALSE])
    
    return(highGR)

}


# different scoring metrics
# remniscent of LOLA: 
# support (number of regions, for us regions that have at least 1 cytosine and can be scored)
# mean loading value (or could do ratio of peak to surrounding regions)
# signal to noise ratio, how big is peak compared to noise of surrounding area
# with SNR, even a small peak could have a high SNR

# @param loadingProf
# loadingProfileSNR = function() {
#     # magnitude of the peak divided by standard deviation of 
#     # noise (signal in surrounding areas)
# }


# support is just number of regions from each input region set that overlap at all with
# the cytosines that have loading values 
# findOverlaps(regionSet, dtToGr(coordinateDT)) 

###########################################################################


# BSaggregate -- Aggregate a BSDT across regions or region groups,
# for multiple samples at a time.
# This function is as BScombineByRegion, but can handle not only multiple
# samples in BSDT, but also simultaneously multiple region sets by passing
# a regionsGRL (GRangesList object).
# you can use jExpr to do other functions.

# Given a bisulfite data table as input, with an identifier column for
# different samples; plus a GRanges objects with regions to aggregate.
#
# @param BSDT The bisulfite data.table (output from one of the parsing
# functions for methylation calls) that you wish to aggregate. It can
# be a combined table, with individual samples identified by column passed
# to splitFactor.
# @param regionsGRL Regions across which you want to aggregate. Should be 
# from a single region set. eg GRangesList(regionSet)
# @param excludeGR A GRanges object with regions you want to 
# exclude from the aggregation function. These regions will be eliminated
# from the input table and not counted.
# @param jExpr You can pass a custom command in the j slot to data.table
# specifying which columns to aggregate, and which functions to use. You
# can use buildJ() to build a jExpr argument easily.
# @param byRegionGroup You can aggregate by regionID or by regionGroupID; 
# this reflects the regionsGRL that you pass; by default, BSAggregate will
# aggregate each region individually -- scores will then be contiguous, and
# the output is 1 row per region.
# Turn on this flag to aggregate across all region groups, making the result
# uncontiguous, and resulting in 1 row per *region group*.
# @param returnSD Whether the standard deviation of the columns of interest 
# should be returned. Standard deviation for rows that overlap with region set
# and also for rows that do not overlap with region set. Not currently functional.
# @param returnOLInfo if true, include number of overlapping cpgs and number 
# of overlapping regions in the result
# @param meanPerRegion Will return the mean value in each region in a data.table
# instead of averaging all regions into a single value, returnOLInfo does
# nothing if meanPerRegion=TRUE (function exits before that code, which
# expects a different data structure)
# @param returnQuantile Only used if meanPerRegion=TRUE, instead of mean
# return the quantile/percentile of the mean of each region
# in relation to the distribution of original values in BSDT
# 
BSAggregate <- function(BSDT, regionsGRL, excludeGR=NULL, 
                        regionsGRL.length = NULL, splitFactor=NULL, 
                        keepCols=NULL, sumCols=NULL, 
                        jExpr=NULL, byRegionGroup=FALSE, 
                        keep.na=FALSE, returnSD=FALSE, 
                        returnOLInfo=FALSE, meanPerRegion=FALSE,
                        returnQuantile=FALSE) {

    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if(is(regionsGRL, "GRanges")) {
        regionsGRL <- GRangesList(regionsGRL)
    } else if (!is(regionsGRL, "GRangesList")) {
        stop("regionsGRL is not a GRanges or GRangesList object")
    }

    # will cause error if BSDT is only a data.frame
    if (is(BSDT, "data.frame") & !is(BSDT, "data.table")) {
        BSDT <- as.data.table(BSDT)
    } 
    if (!is(BSDT, "data.table"))  {
        stop("BSDT must be a data.table")
    }

    if(! is.null(excludeGR)) {
        BSDT <- BSFilter(BSDT, minReads=0, excludeGR)
    }

    if (returnQuantile) {
        # keep all data so quantiles can be calculated
        # later code will change BSDT
        origBSDT <- data.table::copy(BSDT)
    }
    # TODO: BSdtToGRanges needs to include end coordinate!!!
    bsgr <- BSdtToGRanges(list(BSDT))

    additionalColNames <- setdiff(colnames(BSDT), 
                                  c("chr","start", "end",
                                    "hitCount","readCount", splitFactor))

    # specify that "character" outcome is expected from mode by 
    # supplying "a" as last vapply argument (any character object would work)
    colModes <- vapply(BSDT, mode, "a")
    if (is.null(sumCols)) {
        sumCols <- setdiff(colnames(BSDT),c("chr", "start", "end", 
                                            "strand", splitFactor, keepCols))
        # Restrict to numeric columns.		
        sumCols <- intersect(sumCols, 
                             names(colModes[which(colModes == "numeric")]))
    }
    # It's required to do a findoverlaps on each region individually,
    # Not on a GRL, because of the way overlaps with GRLs work. So,
    # we must convert the GRL to a GR, but we must keep track of which
    # regions came from which group.
    regionsGR <- unlist(regionsGRL)

    if(is.null(regionsGRL.length)) {
        # if (length(regionsGRL) > 100) {
        #     message("BSAggregate: Calculating sizes. You can speed this up by supplying a regionsGRL.length vector...", appendLF=FALSE)
        # }
        # specify that output should be numeric with vapply 
        # (any number would work instead of 1)
        regionsGRL.length <- vapply(regionsGRL, length, 1)
        # message("Done counting regionsGRL lengths.")
    }

    # Build a table to keep track of which regions belong to which group
    # BIOC note: sapply returns a list where each item is of different length
    # therefore, I'm not using vapply
    region2group <- data.table(
        regionID=seq_along(regionsGR), 
        chr=as.vector(seqnames(regionsGR)), 
        start=as.vector(start(regionsGR)), 
        end=as.vector(end(regionsGR)),
        withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
        regionGroupID=rep(seq_along(regionsGRL), regionsGRL.length))
    setkey(region2group, regionID)

    #TODO: if overlap method is single, but "end" is present, find center
    # and set that to start!
    if (all(start(bsgr[[1]]) != end(bsgr[[1]]))) {
        stop("BSDT start and end coordinates are not the same. Choose a different aggrMethod.")
    } else {
        fo <- findOverlaps(query = bsgr[[1]], subject = regionsGR)
    }

    ### use info from findOverlaps to see how many individual
    # cytosines (or input regions) and how many regions (in region sets) 
    # overlap with one another
    if (returnOLInfo) {
        signalCoverage    <- length(unique(queryHits(fo)))
        regionSetCoverage <- length(unique(subjectHits(fo)))
    }

    if("end" %in% colnames(BSDT)){
        setkey(BSDT, chr, start, end)
    } else{
        setkey(BSDT, chr, start) 
    }

    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))

    # message("Setting regionIDs...")
    BSDT <- BSDT[queryHits(fo),] #restrict the table to CpGs (or input region) in any region set.

    if (NROW(BSDT) < 1) {
        warning("No BSDT sites in the given region list; please expand your regionsGRL")
        return(NULL)
    }

    BSDT[,regionID:=subjectHits(fo)] #record which region they overlapped.
    #BSDT[queryHits(fo),regionID:=subjectHits(fo)]
    #if (!keep.na) {
    #    BSDT = BSDT[queryHits(fo),]
    #}

    if (is.null(jExpr)) {
        cols=c(sumCols, keepCols)
        funcs <- c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
        jExpr <- buildJ(cols, funcs)
    }
    # message("jExpr: ", jExpr)

    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        agCol <- "regionGroupID"
    } else {
        agCol <- "regionID" # Default
    }

    # Build the by string
    if (is.null(splitFactor)) {
        byString <- paste0("list(regionID)")
    } else {
        byString <- paste0("list(", paste("regionID", 
                                          paste0(splitFactor, ""), 
                                          collapse=", ", sep=", "), ")")
    }

    # Now actually do the aggregate:
    # message("Combining...")
    bsCombined <- BSDT[,eval(parse(text=jExpr)), by=eval(parse(text=byString))]
    setkey(bsCombined, regionID)

    if (meanPerRegion) {
        setkey(region2group, regionID)

        avPerRegion <- merge(bsCombined, region2group)
        avPerRegion[, c("regionID", "withinGroupID", "regionGroupID") := NULL]

        if (returnQuantile) {

            for (i in seq_along(sumCols)) {
                # perhaps this could be more efficient with mapply
                # ecdf example: ecdf(distributionData)(getPercentileOfThis)
                avPerRegion[, c(sumCols[i]) := ecdf(origBSDT[[sumCols[i]]])(avPerRegion[[sumCols[i]]])]
                # I tried set() to improve performance but it took about the same time
            }

        }
        return(avPerRegion)
    }
    # Now aggregate across groups.
    # I do this in 2 steps to avoid assigning regions to groups,
    # which takes awhile. I think this preserves memory and is faster.

    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        # must set allow=TRUE here in case there are multiple IDs (splitCol)
        bsCombined[region2group, regionGroupID:=regionGroupID, allow=TRUE]
        if (! is.null(splitFactor) ) { 
            byStringGroup <- paste0("list(", paste("regionGroupID", 
                                                   paste0(splitFactor, 
                                                          collapse=", "), 
                                                   sep=", "), ")")
        } else {
            byStringGroup <- "list(regionGroupID)"
        }
        bsCombined=bsCombined[,eval(parse(text=jExpr)), 
                              by=eval(parse(text=byStringGroup))]
        if (returnOLInfo) {
            bsCombined[, signalCoverage := signalCoverage]
            bsCombined[, regionSetCoverage := regionSetCoverage]
        }
        return(bsCombined)
    } else {
        e <- region2group[bsCombined,]
        setkey(e, regionID)
        return(e)
    }
    # WARNING: There are now 2^2 ways to aggregate, sum vs mean
    # at each level: across regions, then across region sets. This
    # doesn't give you a choice at this point. 
}



# Given a BSDT (bisulfite data.table), remove any entries that overlap
# regions given in the excludeGR argument and/or filter out sites
# that have lower than a minimum number of reads.
# @param BSDT Bisulfite data.table to filter
# @param minReads Require at least this level of coverage at a cpg.
# @param excludeGR GRanges object with regions to filter.
# 
# @return The BSDT with appropriate regions removed.
BSFilter <- function(BSDT, minReads = 10, excludeGR = NULL) {
    # First, filter for minimum reads.
    if (minReads > 0) {
        BSDT <- BSDT[coverage >= minReads, ]
    }
    if (NROW(BSDT) == 0) { return(data.table(NULL)) }
    # Now, filter entries overlapping a region in excludeGR.
    if (!is.null(excludeGR)) {
        gr <- dtToGr(BSDT)
        fo <- findOverlaps(gr, excludeGR)
        qh <- unique(queryHits(fo))
        length(qh)
        nrow(BSDT)
        BSDT <- BSDT[-qh, ]
    }
    return(BSDT)
}


#' Get indices for top scored region sets 
#' 
#' For each PC, get index of original region sets but ordered by rsScores
#' ranking for each PC. The original index refers to that region set's position
#' in the `GRList` param given to `runCOCOA` which is also that region set's
#' row index in the COCOA output. The first number in a given column 
#' of this function's output will be the
#' original index of the region set ranked first for that PC. Second row for a
#' column will be the original index of the region set that ranked second
#' for that PC, etc. You can use this function to make it easier 
#' when you want to select the top region sets for further analysis or
#' just for sorting the results. Region set scores are sorted in decreasing
#' or increasing order according to the `decreasing` parameter.
#' 
#' @template rsScores
#' @templateVar isRSRankingIndex TRUE
#' @template signalCol
#' @param decreasing logical. Whether to sort rsScores in decreasing 
#' or increasing order. 
#' @param newColName character. The names of the columns of the output data.frame.
#' The order should correspond to the order of the
#'  input columns given by signalCol.
#' @return A data.frame with columns signalCol. Each column has been 
#' sorted by score for region sets for that PC (order given by `decreasing`
#' param).
#' Original indices for region sets that were used to create rsScores
#' are given. Region sets with a score of NA are counted as having the 
#' lowest scores and indices for these region sets will be at the bottom of the
#' returned data.frame (na.last=TRUE in sorting) 
#' @examples data("rsScores")
#' rsRankInd = rsRankingIndex(rsScores=rsScores, 
#'                            signalCol=c("PC1", "PC2"))
#' # region sets sorted by score for PC1
#' rsScores[rsRankInd$PC1, ]
#' # region sets sorted by score for PC2
#' rsScores[rsRankInd$PC2, ]
#' 
#' @export
#' 
rsRankingIndex <- function(rsScores, signalCol, decreasing=TRUE, newColName = signalCol) {
    
    if (!(is(rsScores, "data.frame") || is(rsScores, "matrix"))) {
        stop("rsScores should be a data.frame. Check object class.")
    }
    rsScores <- as.data.table(rsScores)
    
    if (!(is(signalCol, "character") | is(signalCol, "list"))) {
        stop("signalCol should be a character vector or list of character vectors.")
    }
    if (is(signalCol, "list")) {
        if (!all(vapply(X = signalCol, FUN = class, FUN.VALUE = "a") %in% "character")) {
            stop("Items of signalCol should be character vectors.")
        }
        if (!length(unique(vapply(signalCol, FUN = length, FUN.VALUE = 2)))) {
            stop("Items of signalCol should be the same length as each other.")
        }
        if (!all(unlist(signalCol) %in% colnames(rsScores))) {
            stop("Some column names in signalCol are not present in rsScores.")
        }
        
        # the first item in the list is taken for newColName
        if (is(newColName, "list")) {
            newColName <- signalCol[[1]]    
        }
        
    } else { # signalCol is character
        if (!all(signalCol %in% colnames(rsScores))) {
            stop("Some column names in signalCol are not present in rsScores.")
        }
    }
    
    
    dtOrder <- rep(-99L, length(decreasing))
    # how to sort scores
    # -1 for decreasing order of scores
    dtOrder[decreasing] <- -1L
    # +1 for increasing order of scores
    dtOrder[!decreasing] <- 1L
    
    # so by references changes will not be a problem
    rsScores <- copy(rsScores)
    rsScores[, rsIndex := seq_len(nrow(rsScores))]
    
    if (is(signalCol, "list")) {
        if (length(newColName) != length(signalCol[[1]])) {
            stop("newColName is not the same length as columns given in signalCol.")
        }
        
        rsEnSortedInd <- subset(rsScores, select= signalCol[[1]])
        setnames(rsEnSortedInd, newColName)
        
        colNameMat <- do.call(rbind, signalCol) 
        
        # then scores by each PC and make a column with the original index for sorted region sets
        # this object will be used to pull out region sets that were top hits for each PC
        for (i in seq_along(signalCol[[1]])) {
            theseOrderCols <- colNameMat[, i]
            
            setorderv(rsScores, cols = theseOrderCols, order=dtOrder, na.last=TRUE)
            
            rsEnSortedInd[, newColName[i] := rsScores[, rsIndex]]
        }
    } else if (is(signalCol, "character")) {
        
        if (length(newColName) != length(signalCol)) {
            stop("newColName is not the same length as columns given in signalCol.")
        }
        
        signalCol <- signalCol[signalCol %in% colnames(rsScores)]
        
        rsEnSortedInd <- subset(rsScores, select= signalCol)
        setnames(rsEnSortedInd, newColName)
        
        # then scores by each PC and make a column with the original index for sorted region sets
        # this object will be used to pull out region sets that were top hits for each PC
        for (i in seq_along(signalCol)) {
            
            
            setorderv(rsScores, cols = signalCol[i], order=dtOrder, na.last=TRUE)
            
            rsEnSortedInd[, newColName[i] := rsScores[, rsIndex]]
        }
        
    } else {
        stop("signalCol should be a character vector or list of character vectors.")
    }
    
    # reset order
    # setorderv(rsScores, cols = "rsIndex", order=1L)
    return(as.data.frame(rsEnSortedInd))
}

#################### Metric functions ########################################
# scores, metrics, or statistical tests

# Instead of averaging within regions first as BSAggregate does,
# this function does a simple average and standard deviation
# for all CpGs that overlap
# with regions of a region set, also does average and 
# standard deviation for non overlapping CpGs. Created to 
# get metrics of loading values for each PC.
# 
# Faster if given total average for each column of interest
# 
# @param dataDT a data.table with chr, start, end columns as well
# as columns to get metrics of eg (PC1, PC2). All columns
# except chr, start, and end will be considered 
# columns to get the metrics from so no unnecessary columns should be
# included.
# @template regionSet 
# Metrics will be calculated on
# only coordinates within this region set (and optionally separately
# on those outside this region set with alsoNonOLMet parameter)
# @param signalCol the columns to calculate the metrics on. The
# metrics will be calculated on each one of these columns separately.
# @param metrics character vector with the name of a function or functions
# to calculate on selected cytosines. Function should only require one
# input which should be the values of the cytosines.
# @param alsoNonOLMet also include same metrics
# for non overlapping CpGs (still also returns metrics for overlapping CpGs)
# param columnMeans Not functional/is deprecated. The idea was to use
# the mean of the column to speed up calculations (eg when calculating
# mean of overlapping CpGs, use that info and column mean to get
# mean for non overlapping CpGs without manually calculating it)
# 
signalOLMetrics <- function(dataDT,
                            regionSet,
                            signalCol = colnames(dataDT)[!(colnames(dataDT) %in% c("chr", "start", "end"))],
                            metrics=c("mean", "sd"),
                            alsoNonOLMet=TRUE) {
    
    # convert DT to GR for finding overlaps
    dataGR <- BSdtToGRanges(list(dataDT))[[1]]
    
    OL <- findOverlaps(query = regionSet, subject = dataGR)
    
    # if no overlap, exit
    if (length(OL) == 0) {
        return(NULL)
    }
    
    # get indices for overlapping and non overlapping CpGs
    olCpG <- subjectHits(OL)
    
    # region set info
    totalRegionNumber <- length(regionSet)
    meanRegionSize    <- round(mean(width(regionSet)), 1)
    
    
    # get info on degree of overlap
    # number of CpGs that overlap
    signalCoverage <- length(unique(olCpG))
    # number of regions that overlap
    regionSetCoverage <- length(unique(queryHits(OL)))
    
    
    
    
    nonOLCpG <- (seq_len(nrow(dataDT)))[-olCpG]
    
    # gets metrics for all columns except chr, start, end
    jExpr <- buildJ(cols=rep(signalCol, each=length(metrics)),
                    funcs=rep(metrics, length(signalCol)),
                    newColNames = paste0(rep(signalCol,
                                             each=length(metrics)),
                                         "_", metrics))
    
    # getting the metrics
    olMetrics <- as.data.frame(dataDT[olCpG, eval(parse(text=jExpr))])
    
    # if no OL for region set, don't calculate for non region set 
    # TODO make conditional on region set having any overlap
    nonOLMetrics <- as.data.frame(dataDT[nonOLCpG, eval(parse(text=jExpr))])
    
    # calculate average of nonOLCpGs based on columnMean if given
    # if (!is.null())
    #
    # formatting so there is one row per PC/testCol
    # output is a matrix with ncol = length(metrics)
    # for vapply, FUN.VALUE should have length equal to a single output of FUN
    olResults <- vapply(X = metrics,
                        FUN = function(x) as.numeric(olMetrics[, grepl(pattern = x, colnames(olMetrics))]),
                        as.numeric(seq_along(signalCol)))
    olResults <- as.data.table(olResults)
    setnames(olResults, old = colnames(olResults), new = paste0(colnames(olResults), "_OL"))
    
    nonOLResults <- vapply(X = metrics,
                           FUN = function(x) as.numeric(nonOLMetrics[, grepl(pattern = x, colnames(nonOLMetrics))]),
                           as.numeric(seq_along(signalCol)))
    nonOLResults <- as.data.table(nonOLResults)
    setnames(nonOLResults, old = colnames(nonOLResults), new = paste0(colnames(nonOLResults), "_nonOL"))
    
    metricDT <- cbind(data.table(testCol=signalCol),
                      olResults,
                      nonOLResults,
                      data.table(signalCoverage,
                                 regionSetCoverage,
                                 totalRegionNumber,
                                 meanRegionSize))
    
    return(metricDT)
}


# Wilcoxon rank sum test for a region set
# @param dataDT a data.table with chr, start, end columns as well
# as columns to get metrics of eg (PC1, PC2). All columns
# except chr, start, and end will be considered 
# columns to get the metrics from so no unnecessary columns should be
# included.
# @template regionSet
# @param signalCol the columns of interest. You will do ranksum test separately 
# on each of these columns (given test only uses info in one column)
# @param ... Additional parameters of wilcox.test function. See ?wilcox.test.
# For instance specify alternative hypothesis: alternative = "greater".
# @return A vector with a p value for each column other than chr, start or end. 

# @examples data("brcaLoadings1")
# data("brcaMCoord1")
# data("nrf1_chr1")
# dataDT = as.data.table(cbind(brcaMCoord1, brcaLoadings1))
# rsWilcox(dataDT = dataDT, regionSet = nrf1_chr1, conf.int=TRUE)

rsWilcox <- function(dataDT,
                     regionSet,
                     signalCol = colnames(dataDT)[!(colnames(dataDT) %in% c("chr", "start", "end"))], 
                     conf.int=FALSE,
                     ...) {
    
    # region set info
    totalRegionNumber <- length(regionSet)
    meanRegionSize    <- round(mean(width(regionSet)), 1)
    
    
    # convert DT to GR for finding overlaps
    dataGR <- BSdtToGRanges(list(dataDT))[[1]]
    
    OL <- findOverlaps(query = regionSet, subject = dataGR)
    
    # if no overlap, exit
    if (length(OL) == 0) {
        return(NULL)
    }
    
    # get indices for overlapping and non overlapping CpGs
    olCpG <- unique(subjectHits(OL))
    nonOLCpG <- (seq_len(nrow(dataDT)))[-olCpG]
    
    # get info on degree of overlap
    # number of CpGs that overlap
    signalCoverage   <- length(unique(olCpG))
    # number of regions that overlap
    regionSetCoverage   <- length(unique(queryHits(OL)))
    

    
    # each confidence interval has length of 2: [low, high]
    
    if (conf.int) {
        # calculate Wilcoxon rank sum test for each column
        # additional parameters given with ...
        # confIntervals will be [low1, high1, low2, high2, etc.]
        confIntervals <- as.numeric(vapply(X = signalCol, FUN = function(x) wilcox.test(x = as.numeric(as.matrix(dataDT[olCpG, x, with=FALSE])),
                                                                      y = as.numeric(as.matrix(dataDT[nonOLCpG, x, with=FALSE])), 
                                                                      conf.int = conf.int, ...)$conf.int, c(1, 1)))
        
        names(confIntervals) <- paste0(rep(signalCol, each=2), c("_low", "_high"))
        wRes <- data.frame(t(confIntervals),
                           signalCoverage,
                           regionSetCoverage,
                           totalRegionNumber,
                           meanRegionSize)        
            
    } else {
        # calculate Wilcoxon rank sum test for each column
        # additional parameters given with ...
        pVals <- vapply(X = signalCol, FUN = function(x) wilcox.test(x = as.numeric(as.matrix(dataDT[olCpG, x, with=FALSE])),
                                                                      y = as.numeric(as.matrix(dataDT[nonOLCpG, x, with=FALSE])), ...)$p.value, 1)
        wRes <- data.frame(t(pVals),
                           signalCoverage,
                           regionSetCoverage,
                           totalRegionNumber,
                           meanRegionSize)
    }
    
    return(wRes)
}


# @param signalDT Data to be aggregated (e.g. raw data: ATAC-seq,
# region based DNA methylation or loading values)
# @param signalGR GRanges object with coordinates for signalDT
# @template regionSet 
# The region set to score.
# @param calcCols character object. Column names. A weighted sum will be done 
# for each of these columns (columns should be numeric).
# @value Returns data.frame with columns 'calcCols', signalCoverage col has
# number of signalGR regions that overlapped with any regionSet regions, 
# regionSetCoverage has the sum of all proportion overlaps of regions from 
# signalGR with regionSet (regionSet region is denominator)
# containing weighted mean for each col.
# Returns NULL if there is no overlap between signalGR and regionSet

regionOLWeightedMean <- function(signalDT, signalGR, 
                                 regionSet, calcCols) {
    
    
    hits  <- findOverlaps(query = signalGR, subject = regionSet)
    # if no overlap, return NULL
    if (length(hits) == 0) {
        return(NULL)
    }
    
    regionSetCoverage = length(unique(subjectHits(hits))) 
    
    olap  <- pintersect(signalGR[queryHits(hits)],
                        regionSet[subjectHits(hits)])
    polap <- width(olap) / width(regionSet[subjectHits(hits)])
    
    # some rows may be duplicated if a signalDT region overlapped multiple
    # regions from signalGR but that is ok
    signalDT <- signalDT[queryHits(hits), ]
    
    # weight the signalDT values by the proportion overlap (weighted average)
    weightedSum <- t(as.matrix(polap)) %*% as.matrix(signalDT[, calcCols, with=FALSE])
    
    # weighted average
    denom <- sum(polap)
    weightedAve <- as.data.frame(weightedSum / denom)
    # names(weightedAve) <- colsOfInterest
    
    # add columns for coverage info
    weightedAve$signalCoverage = length(unique(queryHits(hits)))
    weightedAve$regionSetCoverage = regionSetCoverage
    weightedAve$sumProportionOverlap = denom

    return(weightedAve)
}


# @param signalDT Data to be aggregated (e.g. raw data: ATAC-seq,
# region based DNA methylation or loading values)
# @param signalGR GRanges object with coordinates for signalDT
# @template regionSet 
# The region set to score.
# @param calcCols character object. Column names. A mean will be calculated for 
# each of these columns (columns should be numeric).
# @value Returns data.frame with columns 'calcCols', signalCoverage col has
# number of signalGR regions that overlapped with any regionSet regions, 
# regionSetCoverage has the number of regions from 
# signalGR that overlapped with regionSet
# Returns NULL if there is no overlap between signalGR and regionSet

regionOLMean <- function(signalDT, signalGR, regionSet, calcCols) {

    hits  <- findOverlaps(query = signalGR, subject = regionSet)
    # if no overlap, return NULL
    if (length(hits) == 0) {
        return(NULL)
    }

    # some rows may be duplicated if a signalDT region overlapped multiple
    # regions from signalGR but that is ok
    signalDT  <- signalDT[queryHits(hits), ]

    # mean of the overlapping signalDT values
    signalAve <- as.data.frame(t(colMeans(signalDT[,..calcCols])))

    # add columns for coverage info
    signalAve$signalCoverage    <- length(unique(queryHits(hits)))
    signalAve$regionSetCoverage <- length(unique(subjectHits(hits)))

    return(signalAve)
}



