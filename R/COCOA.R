# PACKAGE DOCUMENTATION
#' Coordinate Covariation Analysis (COCOA)
#'
#' 
#' COCOA is a method for understanding epigenetic variation among samples.
#' COCOA can be used with epigenetic data that includes 
#' genomic coordinates and an epigenetic signal,
#' such as DNA methylation and chromatin accessibility 
#' data. 
#' To describe the method on a high level, COCOA quantifies 
#' inter-sample variation with either a supervised or unsupervised 
#' technique then uses a database of "region sets" to 
#' annotate the variation among samples. A region set is a set of 
#' genomic regions that share a biological annotation, 
#' for instance transcription factor (TF) binding regions, 
#' histone modification regions, or open chromatin regions. 
#' COCOA can identify region sets that are associated with
#' epigenetic variation between samples and
#' increase understanding of variation in your data.
# 
# In contrast to some other common techniques, COCOA offers both
# supervised (known groups/phenotype) and unsupervised (no known groups/
# phenotype) analyses. Also, COCOA focuses on continuous variation 
# between samples instead of having cutoffs. Because of this, COCOA can 
# be used as a complementary method alongside "differential" methods 
# that find discrete differences between groups of samples and 
# it can also be used in situations where there are no groups.  
# COCOA can identify biologically meaningful 
# sources of variation between samples and 
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
#'             scale_color_brewer theme element_line element_text geom_point
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
NULL

# @importFrom ppcor pcor.test


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
        "coverage", "Group", "pOlap", "regionGroupID", "regionID", "theme", 
        "meanRegionSize", "regionSetCoverage", "rowIndex", "rsIndex",
        "rsRegionID", "totalRegionNumber", "signalCoverage", ".SD",
        "sumProportionOverlap")) 
}

#########################################################################


#' Score a region set using feature contribution scores
#' 
#' First, this function identifies which epigenetic features 
#' overlap the region set. 
#' Then the region set is scored using the feature contribution scores 
#' (`signal` input) 
#' according to the `scoringMetric` parameter.  
#' 
#' @template signal
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
#' @template rsOL
#' @param pOlap Numeric vector. Only used if rsOL is given and scoringMetric
#' is "proportionWeightedMean". This vector should contain the proportion of 
#' each regionSet region that is overlapped by a signalCoord region. The 
#' order of pOlap should be the same as the overlaps in rsOL.
#' @template returnCovInfo
#' @template checkInput

#' @return A data.frame with one row and the following 
#' columns: one column for each item of signalCol with names given
#' by signalCol. These columns have scores for the region set for each signalCol.
#' Other columns: signalCoverage (formerly cytosine_coverage) which
#' has number of epigenetic features that overlapped at all with regionSet,
#' regionSetCoverage which has number of regions from regionSet
#' that overlapped any of the epigenetic features, 
#' totalRegionNumber that has
#' number of regions in regionSet, meanRegionSize that has average
#' size in base pairs of regions in regionSet, the average is based on
#' all regions in regionSet and not just ones that overlap.
#' For "multiBase" data, if the "proportionWeightedMean" scoring metric 
#' is used, then the output will also have a "sumProportionOverlap" column.
#' During this scoring method, the proportion overlap between each signalCoord
#' region and overlapping regionSet region is calculated. This column is
#' the sum of all those proportion overlaps and is another way to quantify
#' coverage of regionSet in addition to regionSetCoverage.
#' 
#' @examples
#' data("brcaATACCoord1")
#' data("brcaATACData1")
#' data("esr1_chr1")
#' featureContributionScores <- prcomp(t(brcaATACData1))$rotation
#' rsScores <- aggregateSignal(signal=featureContributionScores, 
#'                                  signalCoord=brcaATACCoord1, 
#'                                  regionSet=esr1_chr1, 
#'                                  signalCol=c("PC1", "PC2"), 
#'                                  scoringMetric="default")
#' @export

aggregateSignal <- function(signal,
                            signalCoord,
                            regionSet,
                            signalCol = c("PC1", "PC2"),
                            signalCoordType = "default",
                            scoringMetric = "default",
                            verbose = FALSE,
                            absVal=TRUE, 
                            rsOL=NULL, pOlap=NULL,
                            returnCovInfo=TRUE, 
                            .checkInput=FALSE) {
    
    ################### checking inputs  #################################
    
    # if it was already checked outside this function, don't need to re-check
    if (.checkInput) {
        ########## check that inputs are the correct class
        checkConvertInputClasses(signal=signal,
                                 signalCoord=signalCoord,
                                 regionSet=regionSet,
                                 signalCol=signalCol,
                                 rsOL=rsOL)
        
        ########## check that dimensions of inputs are consistent
        # length of signal coord = nrow of signal
        if (length(signalCoord) != nrow(signal)) {
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
        
        if (!(scoringMetric %in% getScoringMethods("both"))) {
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
        # any NAs that overlap the regionSet will cause the score to be NA
        if (is(signal, "data.table")) {
            naRows = apply(X = signal[, signalCol, with=FALSE, drop=FALSE], 
                           MARGIN = 1, FUN = function(x) any(is.na(x)))
        } else {
            naRows = apply(X = signal[, signalCol, drop=FALSE], 
                           MARGIN = 1, FUN = function(x) any(is.na(x)))    
        }
        
        if (any(naRows)) {
            signal <- signal[!naRows, ]
            signalCoord <- signalCoord[!naRows]
            warning("Removing rows with NA from `signal`")
        }
        
        #################################################################
        
        # detect signalCoordType
        if (signalCoordType == "default") {
            
            # when signalCoord is a GRanges object
            if (any(start(signalCoord) != end(signalCoord))) {
                signalCoordType <- "multiBase"
            } else {
                signalCoordType <- "singleBase"
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
            if (!(scoringMetric %in% getScoringMethods("singleBase"))) {
                stop("The scoringMetric you selected is not available for
                 this data's signalCoordType")
            }
        } else if (signalCoordType == "multiBase") {
            if (!(scoringMetric %in% getScoringMethods("multiBase"))) {
                stop("The scoringMetric you selected is not available for
                 this data's signalCoordType")
            }
        }
        
    }
    ################### finished checking inputs #########################
        
    
    
    #### UPDATE: only do this once, in outermost function possible #####
    # extreme positive or negative values both give important information
    # take absolute value or not
    if (absVal) {
        signal <- abs(signal) # required for later code
    }
    
    # XX copies unnecessarily:reformat into data.table with chromosome location and weight
    
    # make sure `signal` is the correct type for further steps
    # (proportionWeightedMean method requires a matrix)
    if (!is(signal, "data.table") && (scoringMetric != "proportionWeightedMean")) {
        signal <- as.data.table(signal)
    } else if (!is(signal, "matrix") && (scoringMetric == "proportionWeightedMean")) {
        signal <- as.matrix(signal)
    }

    # restricting to signalCol so unnecessary computations 
    # are not done
    if (is(signal, "data.table")) {
        loadingDT <- signal[, signalCol, with=FALSE]
        # # naming does not work if only using one PC so add this line for that case
        # setnames(loadiangDT, signalCol)
    } else {
        loadingDT  <- signal[, signalCol, drop=FALSE]
    }
    ######## UPDATE: can above code be done only once, like input checking?
    ########

    
    ###########################################################################
    # scoring

    # works for both singleBase and multiBase
    if (scoringMetric == "simpleMean") {

        loadAgMain <- as.data.table(regionOLMean(signalDT = loadingDT, 
                                                 signalGR = signalCoord,
                                                 regionSet = regionSet,
                                                 calcCols= signalCol,
                                                 metric = "mean",
                                                 rsOL = rsOL,
                                                 returnCovInfo=returnCovInfo))
        results <- .formatResults(loadAgMain, scoringMetric = scoringMetric, 
                                  regionSet=regionSet, signalCol = signalCol,
                                  returnCovInfo=returnCovInfo)

    } else if (scoringMetric == "simpleMedian") {
        # scoring singleBase and multiBase both with this function for
        # simpleMedian
        loadAgMain <- as.data.table(regionOLMean(signalDT = loadingDT, 
                                                 signalGR = signalCoord,
                                                 regionSet = regionSet,
                                                 calcCols= signalCol,
                                                 metric = "median", 
                                                 rsOL=rsOL,
                                                 returnCovInfo = returnCovInfo))
        
        results <- .formatResults(loadAgMain, scoringMetric = scoringMetric, 
                                  regionSet=regionSet, signalCol = signalCol,
                                  returnCovInfo=returnCovInfo)
        
        
    } else if (signalCoordType == "singleBase") {
        # do the actual aggregation
        if (scoringMetric == "regionMean") {

            # specify aggregation operation
            # will be done separately for each PC specified
            aggrCommand <- buildJ(signalCol, 
                                  rep("mean", length(signalCol)))
            # previously used BSAggregate from RGenomeUtils but now using local, 
            # modified copy
            loadAgMain <- BSAggregate(BSDT = loadingDT, 
                                      regionsGRL = regionSet,
                                      BSCoord = signalCoord,
                                      jExpr = aggrCommand,
                                      byRegionGroup = TRUE,
                                      splitFactor = NULL,
                                      returnOLInfo = returnCovInfo, 
                                      rsOL=rsOL)
            
            results <- .formatResults(loadAgMain, 
                                      scoringMetric = scoringMetric, 
                                      regionSet=regionSet, signalCol = signalCol,
                                      returnCovInfo=returnCovInfo)

        } else if (scoringMetric == "regionMedian") {
            
            aggrCommand <- buildJ(signalCol, 
                                  rep("median", length(signalCol)))
            loadAgMain <- BSAggregate(BSDT = loadingDT, 
                                      regionsGRL = regionSet,
                                      BSCoord = signalCoord,
                                      jExpr = aggrCommand,
                                      byRegionGroup = TRUE,
                                      splitFactor = NULL,
                                      returnOLInfo = returnCovInfo)
            
            results <- .formatResults(loadAgMain, 
                                      scoringMetric = scoringMetric,
                                      regionSet=regionSet, signalCol = signalCol,
                                      returnCovInfo=returnCovInfo)

        } 

    } else {
        # signalCoordType == "multiBase"
        # for ATAC-seq
        if (scoringMetric == "proportionWeightedMean") {
            loadAgMain <- regionOLWeightedMean(signalMat = loadingDT, 
                                               signalGR = signalCoord,
                                               regionSet = regionSet,
                                               calcCols= signalCol, 
                                               rsOL = rsOL,
                                               pOlap = pOlap,
                                               returnCovInfo=returnCovInfo)
            results <- .formatResults(as.data.table(loadAgMain), 
                                      scoringMetric = scoringMetric,
                                      regionSet=regionSet, signalCol = signalCol,
                                      returnCovInfo=returnCovInfo)
            
        } 
    }
    
    # signalOLMetrics() # make sure it works with no overlap
    if (verbose) {
        message(":", appendLF=FALSE)
    }
        
    return(as.data.frame(results))
}

# format results of scoring functions in aggregateSignal()
.formatResults <- function(loadAgMain, scoringMetric, 
                           regionSet, signalCol, returnCovInfo=TRUE) {
    
    numOfRegions <- length(regionSet)
        
    # if no cytosines from loadings were included in regionSet, result is NA
    if (is.null(loadAgMain)) {
        results <- as.data.table(t(rep(NA, length(signalCol))))
        setnames(results, signalCol)
        if (returnCovInfo) {
            results[, signalCoverage := 0]
            results[, regionSetCoverage := 0]
            results[, totalRegionNumber := numOfRegions]
            results[, meanRegionSize := round(mean(width(regionSet)), 1)]
            
            # this column is only added by this scoring method
            if (scoringMetric == "proportionWeightedMean") {
                results[, sumProportionOverlap := 0]
            }
        }
    } else {
        # regionMean, regionMedian, simpleMean for region data
        results <- loadAgMain[, .SD, .SDcols = signalCol]
        
        if (returnCovInfo) {
            results[, signalCoverage := loadAgMain[, .SD, .SDcols = "signalCoverage"]]
            results[, regionSetCoverage := loadAgMain[, .SD, .SDcols = "regionSetCoverage"]]
            results[, totalRegionNumber := numOfRegions]
            results[, meanRegionSize := round(mean(width(regionSet)), 1)]
            ################################
            # proportionWeightedMean
            # this column is only added by this scoring method
            if (scoringMetric == "proportionWeightedMean") {
                results[, sumProportionOverlap := loadAgMain[, .SD, .SDcols = "sumProportionOverlap"]]
            }
        }

        ##############################
        # simple mean for base pair data
    }
    return(results)
}


#' Score many region sets
#' 
#' This function will give each region set a score for each target variable
#' given by `signalCol` based on
#' the `scoringMetric` parameter. Based on these scores, you can determine
#' which region sets out of a region set database (given by `GRList`) 
#' are most associated with the target variables. See the vignette "Introduction
#' to Coordinate Covariation Analysis" for help interpreting your 
#' results. 
#'
#' @template signal
#' @template signalCoord
#' @template signalCoordType
#' @template GRList
#' @template signalCol
#' @template scoringMetric
#' @template verbose
# @param wilcox.conf.int logical. Only applies when using "rankSum" scoring
# method. returns a 95% confidence interval from the Wilcoxon rank sum test
# instead of p value.
#' @template absVal
#' @template rsMatList
#' @template returnCovInfo
#' @return Data.frame of results, one row for each region set. 
#' It has the following columns:
#' one column for each item of signalCol with names given
#' by signalCol. These columns have scores for the region set for each signalCol.
#' Other columns: signalCoverage (formerly cytosine_coverage) which
#' has number of epigenetic features that overlapped at all with regionSet,
#' regionSetCoverage which has number of regions from regionSet
#' that overlapped any of the epigenetic features, 
#' totalRegionNumber that has
#' number of regions in regionSet, meanRegionSize that has average
#' size in base pairs of regions in regionSet, the average is based on
#' all regions in regionSet and not just ones that overlap.
#' For "multiBase" data, if the "proportionWeightedMean" scoring metric 
#' is used, then the output will also have a "sumProportionOverlap" column.
#' During this scoring method, the proportion overlap between each signalCoord
#' region and overlapping regionSet region is calculated. This column is
#' the sum of all those proportion overlaps and is another way to quantify
#' coverage of regionSet in addition to regionSetCoverage.
#' 
#' 
#' @examples 
#' data("brcaATACCoord1")
#' data("brcaATACData1")
#' data("esr1_chr1")
#' data("nrf1_chr1")
#' featureContributionScores <- prcomp(t(brcaATACData1))$rotation
#' GRList <- GRangesList(esr1_chr1, nrf1_chr1)
#' rsScores <- aggregateSignalGRList(signal=featureContributionScores, 
#'                                  signalCoord=brcaATACCoord1, 
#'                                  GRList= GRList,
#'                                  signalCol=c("PC1", "PC2"), 
#'                                  scoringMetric="default")
#' 
#' @export

aggregateSignalGRList <- function(signal,
                     signalCoord,
                     GRList,
                     signalCol = c("PC1", "PC2"),
                     signalCoordType = "default",
                     scoringMetric = "default",
                     verbose = TRUE,
                     absVal=TRUE, 
                     rsMatList=NULL, 
                     signalList=NULL, rsInfo=NULL,
                     returnCovInfo=TRUE) {

    ################### checking inputs  #################################
    
    ########## check that inputs are the correct class
    checkConvertInputClasses(signal=signal,
                             signalCoord=signalCoord,
                             regionSet=NULL,
                             signalCol = signalCol,
                             GRList=GRList)#,
                             #olList=olList)
    
    ########## check that dimensions of inputs are consistent
    # length of signal coord = nrow of signal
    if (length(signalCoord) != nrow(signal)) {
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
    
    if (!(scoringMetric %in% getScoringMethods("both"))) {
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
    # any NAs that overlap the regionSet will cause the score to be NA
    if (is(signal, "data.table")) {
        naRows = apply(X = signal[, signalCol, with=FALSE, drop=FALSE], 
                       MARGIN = 1, FUN = function(x) any(is.na(x)))
    } else {
        naRows = apply(X = signal[, signalCol, drop=FALSE], 
                       MARGIN = 1, FUN = function(x) any(is.na(x)))    
    }
    
    if (any(naRows)) {
        signal <- signal[!naRows, ]
        signalCoord <- signalCoord[!naRows]
        warning("Removing rows with NA from `signal`")
    }
    
    #################################################################
    
    # detect signalCoordType
    if (signalCoordType == "default") {
        # when signalCoord is a GRanges object
        if (any(start(signalCoord) != end(signalCoord))) {
            signalCoordType <- "multiBase"
        } else {
            signalCoordType <- "singleBase"
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
    
    if (is.null(signalList)) {
        # data.table is only used for non-matrix COCOA (legacy)
        
        # convert object class outside aggregateSignal to extra prevent copying
        # (one scoring method needs `signal` as a matrix though)
        if (!is(signal, "data.table") && (scoringMetric != "proportionWeightedMean")) {
            signal <- as.data.table(signal)
        } else if (!is(signal, "matrix") && (scoringMetric == "proportionWeightedMean")) {
            signal <- as.matrix(signal)
        }
        
    } else {
        # needed for gamma normalization
        if (!is(signal, "matrix")) {
            signal <- as.matrix(signal)
        }
    }


    # take absolute value outside aggregateSignal to prevent extra copying
    if (absVal) {
        if (is(signal, "data.table")) {
            signal[, signalCol] <- abs(signal[, signalCol, with=FALSE])
        } else {
            signal[, signalCol] <- abs(signal[, signalCol])    
        }
        absVal <- FALSE
    }
    
    # ## if this function is only being run once, use data.table calculations
    # # create region set overlap matrix
    # # this code should be after code modifying "signal"
    # if (scoringMetric %in% c("simpleMean", 
    #                                               "regionMean", 
    #                                               "proportionWeightedMean")) {
    #     if (is.null(rsMatList)) {
    #         olMatRes <- olToMat(signalCoord = signalCoord,
    #                             GRList = GRList, 
    #                             scoringMetric = scoringMetric)
    #         rsMatList <- olMatRes[[1]]
    #         rsInfo <- olMatRes[[2]]
    #     }
    # 
    #     if (is.null(signalList)) {
    #         signalMatList <- splitSignal(signal = signal, 
    #                                      maxRow = nrow(rsMatList[[1]]))
    #     }
    # }
    
    
    if (!is.null(rsMatList) && (scoringMetric %in% c("simpleMean", 
                                                   "regionMean", 
                                                   "proportionWeightedMean"))) {
        resultsDF <- matScore(rsMatList = rsMatList, 
                                signalMatList=signalList, 
                              rsInfo=rsInfo)
    } else { # not matrix COCOA. Old COCOA with data.table
        # apply over the list of region sets
        resultsList <- lapplyAlias(GRList,
                                   function(x) aggregateSignal(
                                       signal = signal,
                                       signalCoord = signalCoord,
                                       signalCoordType = signalCoordType,
                                       regionSet = x,
                                       signalCol = signalCol,
                                       scoringMetric = scoringMetric,
                                       verbose = verbose,
                                       absVal = absVal, pOlap=NULL,
                                       returnCovInfo=returnCovInfo))
    }

    if (!is.null(rsMatList) && (scoringMetric %in% c("simpleMean", 
                                                   "regionMean", 
                                                   "proportionWeightedMean"))) {
        resultsDF[, signalCol] <- gammaNormalize(rsScores=resultsDF[, signalCol],
                                                 featureVals=signal[, signalCol])
        return(as.data.frame(resultsDF))
    } else {
        resultsDT <- do.call(rbind, resultsList) 
        resultsDF <- as.data.frame(resultsDT) 
        if (is(signal, "data.frame")) {
            signal <- as.matrix(signal)
        }
        resultsDF[, signalCol] <- gammaNormalize(rsScores=resultsDF[, signalCol],
                                                 featureVals=signal[, signalCol])
        
        return(resultsDF)   
    }
}    

# @param dataMat If noNA=FALSE, rows of dataMat should be samples/patients, 
# columns should be genomic signal
# (each column corresponds to one genomic coordinate/range). If noNA=TRUE, rows
# and columns should be flipped.
# @param featureMat Matrix. Rows should be samples, columns should be "features" 
# (whatever you want to get correlation with: eg PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns)
# @param testType character object. Can be "cor" (Pearson correlation),
# "spearmanCor (Spearman correlation)
# "pcor" (partial correlation), "cov" (covariance (Pearson)),
# @param covariate
# @param alreadyCenteredDM logical. Whether dataMat has already been centered. If
# known, can do matrix multiplication to calculate covariance (and correlation
# if scaling has been done)
# @param noNA logical. Assume there might be NA, Inf or NaN. If TRUE, the 
# function will use matrix multiplication which is faster than cor/cov
#
# If a row in dataMat has 0 stand. deviation, correlation will be set to 0
# instead of NA as would be done by cor()
#
# @return returns a matrix where rows are the genomic signal (eg a CpG or region) and
# columns are the columns of featureMat
# @examples dataMat = matrix(rnorm(50), 5, 10)
# featureMat = matrix(rnorm(20), 10, 2)
createCorFeatureMat <- function(dataMat, featureMat,
                               testType="cor", covariate=NULL, 
                               alreadyCenteredDM=FALSE, alreadyScaledDM=FALSE,
                               alreadyCenteredFM=FALSE, alreadyScaledFM=FALSE, 
                               noNA=FALSE) {
    
    # featureMat <- as.matrix(featureMat) # copies?
    featureNames <- colnames(featureMat)

    
    # avoid this copy and/or delay transpose until after calculating correlation?
     # dataMat <- as.data.frame(t(dataMat)) # copies, expect transposed form as input instead
    
    # spearman not implemented yet
    if (noNA & (testType != "spearmanCor")) {
        # for matrix correlation, dataMat will be given with rows as features 
        # and cols as samples
        # featureMat has the same orientation for noNA=T or F
        
        # must normalize covariation and correlation by sample number (n-1)
        nSamples <- nrow(featureMat)
        
        if (!alreadyCenteredFM) {
            featureMeans <- colMeans(featureMat, na.rm = TRUE)
            # centering before calculating correlation(also, t() converts to matrix)
            # apply also converts to matrix here I believe
            featureMat <- t(apply(X = featureMat, MARGIN = 1, function(x) x - featureMeans))
            if (dim(featureMat)[1] == 1) {
                featureMat <- t(featureMat)
            }
        }
        if (!alreadyCenteredDM) {
            # more efficient to only do once
            cpgMeans <- rowMeans(dataMat, na.rm = TRUE)
            # centering before calculating correlation
            dataMat <- apply(X = dataMat, MARGIN = 2, function(x) x - cpgMeans)
        }
        if (!is(featureMat, "matrix")) {
            featureMat <- as.matrix(featureMat)
        }
        if (!is(dataMat, "matrix")) {
            dataMat <- as.matrix(dataMat)
        }
        
        ##### cor or cov calculations
        if (testType == "cor") {
            if (!alreadyScaledFM) {
                featureMat <- scale(x = featureMat, center = FALSE, scale = TRUE)
            }
            if (!alreadyScaledDM) {
                featSD <- apply(X = dataMat, MARGIN = 1, FUN = sd)
                featSD[featSD == 0] = 1 # 0 sd if no variation, 1 will leave unchanged and not cause error
                dataMat <- dataMat / featSD
            }
            
            # create feature correlation matrix with PCs (rows: features/CpGs, columns:PCs)
            # how much do features correlate with each PC?
            
            # put epigenetic data first, inner dimensions are samples
            featurePCCor <- (dataMat %*% featureMat) / (nSamples - 1)
            
        } else if (testType == "cov") {
            featurePCCor <- (dataMat %*% featureMat) / (nSamples - 1)
        } else {
            stop("invalid testType")
        }
    } else {
        if (testType == "cor") {
            # create feature correlation matrix with PCs (rows: features/CpGs, columns:PCs)
            # how much do features correlate with each PC?
            
            # put epigenetic data first in cor()
            featurePCCor <- cor(dataMat, featureMat, use="pairwise.complete.obs", method="pearson")
            
        } else if (testType == "spearmanCor") {
            # xtfrm(x) ranking
            featurePCCor <- cor(dataMat, featureMat, use="pairwise.complete.obs", method="spearman")
            
            # } else if (testType == "pcor") {
            #     # partial correlation (account for covariates), ppcor package
            #     
            #     featurePCCor <- apply(X = featureMat, MARGIN = 2, function(y) apply(X = dataMat, 2, 
            #                                                                        FUN = function(x) pcor.test(x = x, y=y,
            #                                                                                                    z=covariate,
            #                                                                                                    method="pearson")$estimate))
            #     
        } else if (testType == "cov") {
            featurePCCor <- cov(dataMat, featureMat, use="pairwise.complete.obs")
            
        } else {
            stop("invalid testType")
        }
        # if standard deviation of the data was zero, NA will be produced
        # set to 0 because no standard deviation means no correlation with attribute of interest
        featurePCCor[is.na(featurePCCor)] <- 0
    }
    
    colnames(featurePCCor) <- featureNames
    
    return(featurePCCor)
    # corLoadRatio <- signal[, signalCol] / featurePCCor 
    # hist(corLoadRatio[, "PC10"])
}


#' Create a "meta-region" profile 
#' 
#' This profile can show enrichment 
#' of genomic signals with high feature contribution scores 
#' in the region set but not in the
#' surrounding genome, suggiiiesting that variation is linked specifically
#' to that region set. 
#' 
#' All regions in a given region set 
#' are combined into a single aggregate profile. Regions in `regionSet` 
#' should be
#' expanded on each side to include a wider area of the genome around
#' the regions of interest (see example and vignettes). 
#' To make the profile, first we optionally take 
#' the absolute value of `signal` (`absVal` parameter). 
#' Then each expanded regionSet region is
#' split into `binNum` bins. The corresponding 
#' bins from each region
#' (e.g. all bin1's, all bin2's, etc.) are grouped.  
#' All overlapping values from `signal` are 
#' aggregated in each bin group according to the `aggrMethod` parameter to 
#' get a meta-region profile. Since DNA strand information is not considered, 
#' the profile is averaged symmetrically around the center.
#' A peak in the middle of this profile suggests
#' that variability is specific to the region set of interest and is 
#' not a product of the surrounding genome. A region set can still be
#' significant even if it does not have a peak. For example, some
#' histone modification region sets may be in large genomic blocks
#' and not show a peak, despite having variation across samples.
#'
#' @template signal
#' @template signalCoord
#' @template regionSet
#' @template signalCol
#' @template signalCoordType
#' @param binNum Number of bins to split each region into when
#' making the aggregate profile. More bins will
#' give a higher resolution but perhaps more noisy profile.
#' @template verbose
#' @templateVar usesAggrMethod TRUE
#' @template scoringMetric
#' @template absVal
#' @return A data.frame with the binned meta-region profile,
#' one row per bin. columns: binID and one column for each target variable
#' in signalCol. The function will return NULL if there
#' is no overlap between signalCoord and any of the bin groups that come 
#' from regionSet (e.g. none of the bin1's overlapped signalCoord, 
#' NULL returned).
#' 
#' @examples 
#' data("brcaATACCoord1")
#' data("brcaATACData1")
#' data("esr1_chr1")
#' featureContributionScores <- prcomp(t(brcaATACData1))$rotation
#' esr1_chr1_expanded <- resize(esr1_chr1, 12000, fix="center")
#' mrProfile <- getMetaRegionProfile(signal=featureContributionScores,
#'                                   signalCoord=brcaATACCoord1,
#'                                   regionSet=esr1_chr1_expanded,
#'                                   signalCol=c("PC1", "PC2"),
#'                                   binNum=21)
#' @export

getMetaRegionProfile <- function(signal, signalCoord, regionSet,
                    signalCol = c("PC1", "PC2"),
                    signalCoordType = "default", 
                    binNum = 21,
                    verbose=TRUE,  
                    aggrMethod = "default", absVal=TRUE) {
    
    ################### checking inputs  #################################
    
    ########## check that inputs are the correct class, convert
    checkConvertInputClasses(signal=signal,
                             signalCoord=signalCoord,
                             regionSet=regionSet,
                             signalCol = signalCol)
    
    ########## check that dimensions of inputs are consistent
    # length of signal coord = nrow of signal
    if (length(signalCoord) != nrow(signal)) {
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
    
    if (!(aggrMethod %in% getScoringMethods("metaRegionProfile"))) {
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
        # when signalCoord is a GRanges object
        if (any(start(signalCoord) != end(signalCoord))) {
            signalCoordType <- "multiBase"
        } else {
            signalCoordType <- "singleBase"
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
    
    GRDT <- grToDt(regionSet)
    
    loadProf <- BSBinAggregate(BSDT = loadingDT,
                               rangeDT = GRDT, 
                               binCount = binNum,
                               BSCoord = signalCoord,
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
# with the PCs to be annotated. One column for the signal of each 
# target variable
# and also has columns with the coordinates for the epigenetic features:
# chr (chromosome) and start column (also possibly end)
# @param rangeDT A data.table with the sets of regions to be binned, 
# with columns named start, end
# @param binCount Number of bins across the region
# @param BSCoord GRanges. Coordinates for BSDT. If NULL, then "chr" and "start"
# columns must be in BSDT.
# @param byRegionGroup Pass along to binCount (see ?binCount)
# @template signalCol
# @param verbose A "logical" object. Whether progress 
# of the function should be shown, one
# bar indicates the region set is completed.
# useful when using BSBinAggregate with 'apply' to do many 
# region sets at a time.
# @param aggrMethod see ?getMetaRegionProfile()
BSBinAggregate <- function(BSDT, rangeDT, binCount,
                           BSCoord=NULL,
                           byRegionGroup = TRUE,
                           splitFactor = NULL,
                           signalCol,
                           verbose = FALSE,
                           aggrMethod) {
    if (!is(rangeDT, "data.table")) {
        stop("rangeDT must be a data.table")
    }
    
    if (is.null(BSCoord)) {
        BSCoord <- BSdtToGRanges(list(BSDT))
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
                               FUN = function(x) regionOLWeightedMean(signalMat = BSDT, 
                                                                      # signalGR = dtToGr(BSDT[, .(chr, start, end)]), 
                                                                      signalGR = BSCoord,
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
                                                              # signalGR = dtToGr(BSDT[, .(chr, start, end)]),
                                                              signalGR = BSCoord,
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
                                  BSCoord = BSCoord, 
                                  jExpr=buildJ(signalCol,
                                               rep("mean", length(signalCol))),
                                  byRegionGroup = byRegionGroup,
                                  splitFactor = splitFactor)
        
        # if any bins had no data
        if (nrow(binnedBSDT) < binCount) {
            return(NULL)
        }
    } else if (aggrMethod == "regionMedian") {
        # what is output if a region set has no overlap?
        binnedBSDT <- BSAggregate(BSDT,
                                  regionsGRL=GRangesList(binnedGR),
                                  BSCoord = BSCoord, 
                                  jExpr=buildJ(signalCol,
                                               rep("median", length(signalCol))),
                                  byRegionGroup = byRegionGroup,
                                  splitFactor = splitFactor)
        
        # if any bins had no data
        if (nrow(binnedBSDT) < binCount) {
            return(NULL)
        }
    }
    # RGenomeUtils::BSAggregate

    # # If we aren't aggregating by bin, then don't restrict to min reads!
    # if (byRegionGroup) {
    #     binnedBSDT <- binnedBSDT[readCount > minReads,]
    # }
    if (verbose) {
        message(".", appendLF=FALSE)
    }
    
    return(binnedBSDT)
}

# modification of BSAggregate to just return mean per region
# 
# @template signal
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
# target variable in signalCol. 
# Regions are not in order along the rows of the data.table.
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
    
    ########## check that inputs are the correct class and converts
    checkConvertInputClasses(signal=signal,
                             signalCoord=signalCoord,
                             regionSet=regionSet,
                             signalCol = signalCol)
    
    ########## check that dimensions of inputs are consistent
    # length of signal coord = nrow of signal
    if (length(signalCoord) != nrow(signal)) {
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
    # any NAs that overlap the regionSet will cause the score to be NA
    if (is(signal, "data.table")) {
        naRows = apply(X = signal[, signalCol, with=FALSE, drop=FALSE], 
                       MARGIN = 1, FUN = function(x) any(is.na(x)))
    } else {
        naRows = apply(X = signal[, signalCol, drop=FALSE], 
                       MARGIN = 1, FUN = function(x) any(is.na(x)))    
    }
    
    if (any(naRows)) {
        signal <- signal[!naRows, ]
        signalCoord <- signalCoord[!naRows]
        warning("Removing rows with NA from `signal`")
    }
    
    #################################################################
    
    #################################################################

    # determine whether coordinates are single base or a range
    # if (!("end" %in% colnames(coordinateDT))) {
    #     dataCoordType <- "singleBase"
    # } else {
        if (any(start(signalCoord) != end(signalCoord))) {
            dataCoordType <- "multiBase"
        } else {
            dataCoordType <- "singleBase"
        }
    # }

    # take absolute value or not
    if (absVal) {
        signalDT <- as.data.table(abs(signal))
    } else {
        signalDT <- as.data.table(signal)
    }

    # use different function for single base data and for region data
    if (dataCoordType == "singleBase") {
        # linking coordinates to loading values, has columns chr start, signalCol
        BSDT  <- signalDT[, .SD, .SDcols = signalCol]
        jExpr <- buildJ(signalCol, rep("mean", length(signalCol)))

        avPerRegion <- BSAggregate(BSDT = BSDT,
                                   regionsGRL = regionSet,
                                   BSCoord = signalCoord,
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




#' Get regions that are most associated with target variable
#'
#' Get a GRanges with top regions from the region set based on 
#' average feature contribution scores
#' for the regions or the quantile of the region's average
#' feature contribution score based on the 
#' distribution of all feature contribution scores for the target variable.
#' Returns average feature contribution score or quantile as GRanges metadata.
#' 
#' @template signal
#' @template signalCoord
#' @template regionSet
#' @template signalCol
#' @param cutoff Numeric. Only regions with at least this value will be 
#' returned (either above this average `signal` value or above this quantile
#' if returnQuantile=TRUE).
#' @param returnQuantile Logical. If FALSE, return region averages. If TRUE,
#' for each region, return the quantile of that region's average value
#' based on the distribution of individual feature values in `signal` for
#' that `signalCol`.
#' @return A GRanges object with region coordinates for regions with
#' scores/quantiles above "cutoff" for any target variable in signalCol. 
#' The scores/quantiles
#' for signalCol are given as metadata in the GRanges.

# Are regions in order along the rows of the data.table?
#
#' @examples 
#' data("brcaATACCoord1")
#' data("brcaATACData1")
#' data("esr1_chr1")
#' featureContributionScores <- prcomp(t(brcaATACData1))$rotation
#' topRegions <- getTopRegions(signal=featureContributionScores,
#'                             signalCoord=brcaATACCoord1,
#'                             regionSet=esr1_chr1,
#'                             returnQuantile = TRUE)
#' @export

getTopRegions <- function(signal, 
                          signalCoord, 
                          regionSet, 
                          signalCol = c("PC1", "PC2"), 
                          cutoff = 0.8, 
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
# @param BSCoord GRanges. Coordinates for BSDT. If NULL, then "chr" and "start"
# columns must be in BSDT.
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
# @template rsOL
# 
BSAggregate <- function(BSDT, regionsGRL, BSCoord=NULL, excludeGR=NULL, 
                        regionsGRL.length = NULL, splitFactor=NULL, 
                        keepCols=NULL, sumCols=NULL, 
                        jExpr=NULL, byRegionGroup=FALSE, 
                        keep.na=FALSE, returnSD=FALSE, 
                        returnOLInfo=FALSE, meanPerRegion=FALSE,
                        returnQuantile=FALSE, rsOL=NULL) {

    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if(is(regionsGRL, "GRanges")) {
        regionsGRL <- GRangesList(regionsGRL)
    } else if (!is(regionsGRL, "GRangesList") && is.null(rsOL)) {
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
    if (is.null(BSCoord)) {
        bsgr <- BSdtToGRanges(list(BSDT))
    } else {
        bsgr <- BSCoord
    }
    

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
    
    if (is.null(rsOL)) {
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
        
        # one regionID per region of "regionsGR", withinGroupID: for when regionsGRL
        # has multiple GRanges, Gives the index of each region within its corresponding
        # GRanges object. regionGroupID: for when regionsGRL has multiple GRanges,
        # indicates which regionsGRL GRanges object the region has come from (
        # index of the GRanges object in regionsGRL)
        region2group <- data.table(
            regionID=seq_along(regionsGR), 
            chr=as.vector(seqnames(regionsGR)), 
            start=as.vector(start(regionsGR)), 
            end=as.vector(end(regionsGR)),
            # withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
            regionGroupID=rep(seq_along(regionsGRL), regionsGRL.length)) # repeats each number regionGRL.length times
        setkey(region2group, regionID)
        
        #TODO: if overlap method is single, but "end" is present, find center
        # and set that to start!
        # if (all(start(bsgr[[1]]) != end(bsgr[[length(unique(queryHits(hits)))1]]))) {
        if (all(start(bsgr) != end(bsgr))) {
            stop("BSDT start and end coordinates are not the same. Choose a different aggrMethod.")
        } else {
            # fo <- findOverlaps(query = bsgr[[1]], subject = regionsGR)
            fo <- findOverlaps(query = bsgr, subject = regionsGR)
        }
        
        if (length(subjectHits(fo)) < 1) {
            warning("Insufficient overlap between signalCoord and the region set.")
            return(NULL)
        }
        
    } else { 
        # Build a table to keep track of which regions belong to which group
        # BIOC note: sapply returns a list where each item is of different length
        # therefore, I'm not using vapply
        
        # region2group <- data.table(
        #     regionID=seq_len(regionsGR), 
        #     withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
        #     regionGroupID=rep(seq_along(regionsGRL), regionsGRL.length))
        
        if (length(subjectHits(rsOL)) < 1) {
            warning("Insufficient overlap between signalCoord and the region set.")
            return(NULL)
        }
        
        # only works for when one region set is given in regionsGRL (ie does
        # not work for metaregion profiles)
        region2group <- data.table(
            regionID=seq_len(max(subjectHits(rsOL))),
            regionGroupID=rep(1, max(subjectHits(rsOL)))) # assumes only 1 region set in regionsGRL
        setkey(region2group, regionID)
        
        fo <- rsOL
    }
    

    ### use info from findOverlaps to see how many individual
    # cytosines (or input regions) and how many regions (in region sets) 
    # overlap with one another
    if (returnOLInfo) {
        signalCoverage    <- length(unique(queryHits(fo)))
        regionSetCoverage <- length(unique(subjectHits(fo)))
    }

    # if("end" %in% colnames(BSDT)){
    #     setkey(BSDT, chr, start, end)
    # } else{
    #     setkey(BSDT, chr, start) 
    # }

    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))

    # message("Setting regionIDs...")
    BSDT <- BSDT[queryHits(fo),] #restrict the table to CpGs (or input region) in any region set.

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

    # Now actually do the aggregate: (aggregate within each region)
    # message("Combining...")
    bsCombined <- BSDT[,eval(parse(text=jExpr)), by=eval(parse(text=byString))]
    setkey(bsCombined, regionID)

    if (meanPerRegion) {
        setkey(region2group, regionID)

        avPerRegion <- merge(bsCombined, region2group)
        avPerRegion[, c("regionID", 
                        # "withinGroupID", 
                        "regionGroupID") := NULL]

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
    # byRegionGroup=TRUE means to aggregate within each individiual GRanges from
    # regions.GRL
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
        # aggregate by regionGroupID (region set/GRanges in regionsGRL)
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
#' For each target variable, get index of original region sets 
#' but ordered by rsScores
#' ranking for each target variable. 
#' The original index refers to that region set's position
#' in the `GRList` param given to `aggregateSignalGRList` which is 
#' also that region set's
#' row index in the COCOA output. The first number in a given column 
#' of this function's output will be the
#' original index of the region set ranked first for that target variable.
#' The second row for a
#' column will be the original index of the region set that ranked second
#' for that target variable, etc. You can use this function to make it easier 
#' when you want to select the top region sets for further analysis or
#' just for sorting the results. Region set scores are sorted in decreasing
#' or increasing order according to the `decreasing` parameter.
#' 
#' @template rsScores
#' @templateVar isRSRankingIndex TRUE
#' @template signalCol
#' @param decreasing Logical. Whether to sort rsScores in decreasing 
#' or increasing order. 
#' @param newColName Character. The names of the columns of the output data.frame.
#' The order should correspond to the order of the
#'  input columns given by signalCol.
#' @return A data.frame with one column for each `signalCol`. 
#' Column names are given by `signalCol` or `newColName` (if used). 
#' Each column has been 
#' sorted by score for region sets for that target variable 
#' (order given by `decreasing`
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
rsRankingIndex <- function(rsScores, signalCol, 
                           decreasing=TRUE, newColName = signalCol) {
    
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
# @param dataDT a data.table with
# columns to get metrics of eg (PC1, PC2). All columns
# will be considered 
# columns to get the metrics from so no unnecessary columns should be
# included.
# @param dataGR GRanges. Coordinates for dataDT.
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
                            dataGR,
                            regionSet,
                            signalCol = colnames(dataDT)[!(colnames(dataDT) %in% c("chr", "start", "end"))],
                            metrics=c("mean", "sd"),
                            alsoNonOLMet=TRUE, rsOL=NULL) {
    
    # convert DT to GR for finding overlaps
    # dataGR <- BSdtToGRanges(list(dataDT))[[1]]
    
    if (is.null(rsOL)) {
        OL <- findOverlaps(query = dataGR, subject = regionSet)   
        # region set info
        totalRegionNumber <- length(regionSet)
        meanRegionSize    <- round(mean(width(regionSet)), 1)
    } else {
        OL <- rsOL
    }
    
    # if no overlap, exit
    if (length(OL) == 0) {
        return(NULL)
    }
    
    # get indices for overlapping and non overlapping CpGs
    olCpG <- queryHits(OL)
    
    # get info on degree of overlap
    # number of CpGs that overlap
    signalCoverage <- length(unique(olCpG))
    # number of regions that overlap
    regionSetCoverage <- length(unique(subjectHits(OL)))
    
    # gets metrics for all columns except chr, start, end
    jExpr <- buildJ(cols=rep(signalCol, each=length(metrics)),
                    funcs=rep(metrics, length(signalCol)),
                    newColNames = paste0(rep(signalCol,
                                             each=length(metrics)),
                                         "_", metrics))
    
    # getting the metrics
    olMetrics <- as.data.frame(dataDT[olCpG, eval(parse(text=jExpr))])
    
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
    
    if (alsoNonOLMet) {
        nonOLCpG <- (seq_len(nrow(dataDT)))[-olCpG]
        # if no OL for region set, don't calculate for non region set 
        # TODO make conditional on region set having any overlap
        nonOLMetrics <- as.data.frame(dataDT[nonOLCpG, eval(parse(text=jExpr))])
        
        nonOLResults <- vapply(X = metrics,
                               FUN = function(x) as.numeric(nonOLMetrics[, grepl(pattern = x, colnames(nonOLMetrics))]),
                               as.numeric(seq_along(signalCol)))
        nonOLResults <- as.data.table(nonOLResults)
        setnames(nonOLResults, old = colnames(nonOLResults), new = paste0(colnames(nonOLResults), "_nonOL"))
        
        if (is.null(rsOL)) {
            metricDT <- cbind(data.table(testCol=signalCol),
                              olResults,
                              nonOLResults,
                              data.table(signalCoverage,
                                         regionSetCoverage,
                                         totalRegionNumber,
                                         meanRegionSize))
        } else {
            metricDT <- cbind(data.table(testCol=signalCol),
                              olResults,
                              nonOLResults,
                              data.table(signalCoverage,
                                         regionSetCoverage))
        }
    } else {
        if (is.null(rsOL)) {
            metricDT <- cbind(data.table(testCol=signalCol),
                              olResults,
                              data.table(signalCoverage,
                                         regionSetCoverage,
                                         totalRegionNumber,
                                         meanRegionSize))
        } else {
            metricDT <- cbind(data.table(testCol=signalCol),
                              olResults,
                              data.table(signalCoverage,
                                         regionSetCoverage))
        }

    }
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


# I optimized upstream so that a matrix would be given to this function
# if this function is rewritten and no longer requires a matrix input,
# then in order to prevent unnecessary object copying, 
# rewrite upstream code that converts signalDT to matrix class

# @param signalMat Data to be aggregated (e.g. raw data: ATAC-seq,
# region based DNA methylation or loading values)
# @param signalGR GRanges object with coordinates for signalMat
# @template regionSet 
# The region set to score.
# @param calcCols character object. Column names. A weighted sum will be done 
# for each of these columns (columns should be numeric).
# @template rsOL
# @param pOlap see "?aggregateSignal"
# @template returnCovInfo
# @value Returns data.frame with columns 'calcCols', signalCoverage col has
# number of signalGR regions that overlapped with any regionSet regions, 
# regionSetCoverage has the sum of all proportion overlaps of regions from 
# signalGR with regionSet (regionSet region is denominator)
# containing weighted mean for each col.
# Returns NULL if there is no overlap between signalGR and regionSet

regionOLWeightedMean <- function(signalMat, signalGR, 
                                 regionSet, calcCols, rsOL=NULL,
                                 pOlap=NULL, returnCovInfo=TRUE) {
    
    if (!is(signalMat, "matrix")) {
        signalMat <- as.matrix(signalMat)
    }
    
    if (is.null(rsOL)) {
        hits  <- findOverlaps(query = signalGR, subject = regionSet) 
    } else {
        hits <- rsOL
        
    }
    
    # if no overlap, return NULL
    if (length(hits) == 0) {
        return(NULL)
    }
    
    if (is.null(pOlap)) {
        olap  <- pintersect(signalGR[queryHits(hits)],
                            regionSet[subjectHits(hits)])
        pOlap <- width(olap) / width(regionSet[subjectHits(hits)])
    }

    
    # some rows may be duplicated if a signalMat region overlapped multiple
    # regions from signalGR but that is ok
    # signalMat <- signalMat[queryHits(hits), ] # done in next step to prevent extra copying
    
    # weight the signalMat values by the proportion overlap (weighted average)
    weightedSum <- t(pOlap) %*% signalMat[queryHits(hits), calcCols]
    
    # weighted average
    denom <- sum(pOlap)
    weightedAve <- as.data.frame(weightedSum / denom)
    colnames(weightedAve) <- calcCols
    
    # add columns for coverage info
    if (returnCovInfo) {
        weightedAve$signalCoverage = length(unique(queryHits(hits)))
        weightedAve$regionSetCoverage = length(unique(subjectHits(hits))) 
        weightedAve$sumProportionOverlap = denom
    }

    return(weightedAve)
}


# @param signalDT Data to be aggregated (e.g. raw data: ATAC-seq,
# region based DNA methylation or loading values)
# @param signalGR GRanges object with coordinates for signalDT
# @template regionSet 
# The region set to score.
# @param calcCols character object. Column names. A mean will be calculated for 
# each of these columns (columns should be numeric).
# @param metric character. "mean" or "median"
# @template rsOL
# @template returnCovInfo
# @value Returns data.frame with columns 'calcCols', signalCoverage col has
# number of signalGR regions that overlapped with any regionSet regions, 
# regionSetCoverage has the number of regions from 
# signalGR that overlapped with regionSet
# Returns NULL if there is no overlap between signalGR and regionSet

regionOLMean <- function(signalDT, signalGR, regionSet, 
                         calcCols, metric="mean", rsOL=NULL, returnCovInfo=TRUE) {

    if (is.null(rsOL)) {
        hits  <- findOverlaps(query = signalGR, subject = regionSet)
    } else {
        hits <- rsOL
    }
        # if no overlap, return NULL
    if (length(hits) == 0) {
        return(NULL)
    }

    # some rows may be duplicated if a signalDT region overlapped multiple
    # regions from signalGR but that is ok


    if (metric == "mean") {
        # mean of the overlapping signalDT values
        signalAve <- as.data.frame(t(colMeans(signalDT[queryHits(hits),..calcCols])))
    } else if (metric == "median") {
        # median of the overlapping signalDT values
        signalAve <- as.data.frame(t(apply(X = signalDT[queryHits(hits),..calcCols], 2, median)))
        
    } else {
        stop("Error in regionOLMean function. Invalid metric specified.")
    }
    
    if (returnCovInfo) {
        # add columns for coverage info
        signalAve$signalCoverage    <- length(unique(queryHits(hits)))
        signalAve$regionSetCoverage <- length(unique(subjectHits(hits)))
    }

    return(signalAve)
}

##########################################################################
# matrix scoring
# 1. make a region set matrix. The dimensions are nrows=nfeatures of 
# epigenetic data segmentation (e.g. ATAC consensus peaks), ncol=nregionsets. 
# It has a 1 for data regions that are overlapped by a given region set 
# and a zero for data regions that do not overlap the region set
# This will produce an unweighted mean.
# 2. Multiply the region set matrix times the loading/correlation matrix.
# 3. Divide by total covered regions for that region set to get the mean. This
# is the region set score.

# @template signalCoord
# @template GRList
# @template scoringMetric
# @value Returns a list with two items: 1. a weight matrix where each 
# column corresponds to one region set
# and rows are data regions. 2. A data.frame with size and coverage information
# about the region sets.

# Matrix weights depend on the scoringMetric
# All signalCoord regions and region set regions should be run through function
# at the same time so that there will be proper weighting if a region set 
# region overlaps multiple signalCoord regions and a scoringMetric other than
# "simpleMean" is being used.
# @examples data("brcaMCoord1")
# data("nrf1_chr1")
# data("esr1_chr1")
# myGRL <- GRangesList(esr1_chr1, nrf1_chr1)
# res = olToMat(signalCoord=, GRList=myGRL, scoringMetric="regionMean")
# olMatList <- res[[1]]
# coverageInfo <- res[[2]]
# 
olToMat = function(signalCoord, GRList, scoringMetric, 
                   maxRow=500000, minRSCov=0) {

    # input checks
    if (is.null(names(GRList))) {
        names(GRList) <- paste0("regionSet", seq_along(GRList))
    }
    ###########################################################################

    # calculate overlaps only once
    # region set must be subject to fit with scoring functions
    olList <- lapply(X = GRList, FUN = function(x) findOverlaps(query = signalCoord, 
                                                                subject = x))
    
    regionSetCoverage <- rep(0, length(GRList))
    signalCoverage <- rep(0, length(GRList))
    
    regionSetCoverage <- vapply(X = olList, 
                                FUN = function(x) length(unique(subjectHits(x))), 
                                FUN.VALUE = -1)
    keepInd <- which(regionSetCoverage >= minRSCov)
    rsNames <- names(GRList)[keepInd]
    regionSetCoverage <- regionSetCoverage[keepInd]
    signalCoverage <- vapply(X = olList[keepInd], 
                             FUN = function(x) length(unique(queryHits(x))), 
                             FUN.VALUE = -1)
    
    # # each column is a region set
    # rsMat <- matrix(data = rep(0, length(GRList) * length(signalCoord)), 
    #                nrow=length(signalCoord))
    # 
    nMat <- ceiling(length(signalCoord) / maxRow)
    finalMatRows <- length(signalCoord) %% maxRow
    if (finalMatRows == 0) {
        finalMatRows = maxRow
    }
    # rep(x, 0) is fine if nMat=1
    matRowNVec = c(rep(maxRow, nMat-1), finalMatRows)
    rsMatList <- list()
    if (nMat > 1) {
        for (matCount in 1:(nMat-1)) {
            rsMatList[[matCount]] = matrix(data=rep(0, maxRow * length(keepInd)), 
                                                    nrow=maxRow)
            colnames(rsMatList[[matCount]]) <- rsNames
        }
        # final matrix might not be the same size
        rsMatList[[nMat]] = matrix(data=rep(0, finalMatRows * length(keepInd)), 
                                                nrow=finalMatRows)
        colnames(rsMatList[[nMat]]) <- rsNames
  
    } else { # nMat = 1
        rsMatList[[nMat]] = matrix(data=rep(0, finalMatRows * length(keepInd)), 
                                            nrow=finalMatRows) 
        colnames(rsMatList[[nMat]]) <- rsNames
    }
    

    
    #######################################################################
    # calculate scores and assign to matrices in chunks
    if (scoringMetric == "simpleMean") {
        for (i in seq_along(keepInd)) {
            
            if (signalCoverage[i] != 0) {
                tmpCount = 1
                tmpVec = rep(0, length(signalCoord))
                # normalize by number so that matrix multiplication during COCOA scoring will produce mean
                tmpVec[unique(queryHits(olList[[keepInd[i]]]))] <- 1 / signalCoverage[i]  
                
                # for coordinate chunks
                for (j in seq_along(rsMatList)) {
                    rsMatList[[j]][, i] <- tmpVec[tmpCount:(tmpCount + matRowNVec[j]-1)]
                    tmpCount = tmpCount + matRowNVec[j]
                } 
            }
        }
        # put in correct orientation for matrix multiplication to prevent transposition
        rsMatList <- lapply(rsMatList, FUN = t) # transpose
    } else if (scoringMetric == "regionMean") {
        # instead of 1 give weight proportional to how many signalCoord are 
        # overlapped 
        
        for (i in seq_along(keepInd)) {
            if (signalCoverage[i] != 0) {
                tmpCount = 1
                tmpVec = rep(0, length(signalCoord))
                # aggregate number of overlaps by region set region
                tmp <- as.data.table(olList[[keepInd[i]]])
                # want the count per region set number
                tmp[, rCount := (1/.N), , by=subjectHits]
                normFactor <- sum(tmp$rCount)
                tmpVec[tmp$queryHits] <- tmp$rCount / normFactor
                
                # for coordinate chunks
                for (j in seq_along(rsMatList)) {
                    rsMatList[[j]][, i] <- tmpVec[tmpCount:(tmpCount + matRowNVec[j]-1)]
                    tmpCount = tmpCount + matRowNVec[j]
                }
            }
        }
        # put in correct orientation for matrix multiplication to prevent transposition
        rsMatList <- lapply(rsMatList, FUN = t) # transpose
        
    } else if (scoringMetric == "proportionWeightedMean") {
        sumProportionOverlap = rep(0, length(keepInd))
        for (i in seq_along(keepInd)) {
            if (signalCoverage[i] != 0) {
                tmpCount = 1
                tmpVec = rep(0, length(signalCoord))
                olap  <- pintersect(GRList[[keepInd[i]]][subjectHits(olList[[keepInd[i]]])],
                                    signalCoord[queryHits(olList[[keepInd[i]]])])
                pOlap <- width(olap) / width(GRList[[keepInd[i]]][subjectHits(olList[[keepInd[i]]])])
                
                # weighted average
                denom <- sum(pOlap)
                sumProportionOverlap[i] = denom
    
                # aggregate pOlap by signalCoord region
                olDT <- data.table(queryHits = queryHits(olList[[keepInd[i]]]), 
                                   pOlap=pOlap)
                normDT <- olDT[, .(coordSum = sum(pOlap)), by=queryHits]
                tmpVec[normDT$queryHits] <- normDT$coordSum / denom
                
                # for coordinate chunks
                for (j in seq_along(rsMatList)) {
                    rsMatList[[j]][, i] <- tmpVec[tmpCount:(tmpCount + matRowNVec[j]-1)]
                    tmpCount = tmpCount + matRowNVec[j]
                }
                
                
                # below comment outdated?
                # some rows may be duplicated if a signalMat region overlapped multiple
                # regions from signalGR but that is ok
                # signalMat <- signalMat[queryHits(hits), ] # done in next step to prevent extra copying
            }
        }
        # put in correct orientation for matrix multiplication to prevent transposition
        rsMatList <- lapply(rsMatList, FUN = t) # transpose
    } else {
        stop("The given scoringMetric cannot be used with matrix scoring.")
    }
    
    totalRegionNumber <- vapply(X = GRList[keepInd], FUN = length, FUN.VALUE = -1)
    meanRegionSize <- vapply(X = GRList[keepInd], FUN = function(x) mean(width(x)), 
                             FUN.VALUE = -1)
    # list item 2
    if (scoringMetric == "proportionWeightedMean") {
        rsInfo = data.frame(rsName=rsNames, signalCoverage, 
                            regionSetCoverage, sumProportionOverlap, 
                            totalRegionNumber, 
                            meanRegionSize)
    } else {
        rsInfo = data.frame(rsName=rsNames, signalCoverage, 
                            regionSetCoverage, totalRegionNumber, 
                            meanRegionSize)
    }

    # overlap matrix list and region set coverage info as data.frame 
    return(list(rsMatList, rsInfo))
    # totalRegionNumber = sapply(X = GRList, length)
    # meanRegionSize = sapply(X = GRList, function(x) round(mean(width(x))))
}

# multiply region set matrices with data matrices to get COCOA score
# rsMatList mats should have rows as region sets and cols as features
# ^ so data won't have to be copied during transpose
matScore <- function(rsMatList, signalMatList, rsInfo) {
    
    # rsMatList mats are features X region sets
    # signalMatList is features X target variables
    # multiply matrices
    scoreL <- mapply(FUN = function(x, y) x %*% y, 
                    x=rsMatList, y=signalMatList, SIMPLIFY = FALSE)
    # each item in scoreL is region sets X target variables
    # combine results. Already normalized to one so can just add to get mean.
    scoreDF <- Reduce("+", scoreL)
    
    resultsDF <- cbind(scoreDF, rsInfo)
    return(resultsDF)
}

splitSignal <- function(signal, maxRow=500000) {
    
    signalList <- list()
    nMat <- ceiling(nrow(signal) / maxRow)

    finalMatRows <- nrow(signal) %% maxRow
    if (finalMatRows == 0) {
        finalMatRows = maxRow
    }
    # rep(x, 0) is fine if nMat=1
    matRowNVec = c(rep(maxRow, nMat-1), finalMatRows)
    tmpCount = 1
    for (j in seq(nMat)) {
        signalList[[j]] <- signal[tmpCount:(tmpCount + matRowNVec[j]-1), ]
        tmpCount = tmpCount + matRowNVec[j]
    }
    
    return(signalList)
}


