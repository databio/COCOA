# PACKAGE DOCUMENTATION
#' Principal Component Region Set Analysis (PCRSA)
#' 
#'A method for annotating principal components of DNA methylation data 
#'with region sets. A region set is a set of genomic regions that share 
#'a biological annotation. PCRSA can identify biologically meaningful 
#'sources of variation along the PCs, increasing interpretability of the PCs 
#'and understanding of variation in the data. 
#'
#' @docType package
#' @name PCRSA
#' @author John Lawson
#' @author Nathan Sheffield
#'
#' @references \url{http://github.com/databio}
#' @importFrom GenomicRanges GRanges GRangesList elementMetadata strand
#'             seqnames granges
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_jitter geom_line
#'             theme_classic xlab ylab geom_hline ylim scale_color_discrete
#'             scale_x_discrete scale_fill_brewer scale_color_manual
#'             scale_color_brewer
#' @importFrom ComplexHeatmap Heatmap draw
#' @import BiocGenerics S4Vectors IRanges
#' @importFrom data.table ":=" setDT data.table setkey fread setnames 
#'             setcolorder rbindlist setattr setorder copy is.data.table 
#'             setorderv as.data.table
#' @importFrom Biobase sampleNames
#' @importFrom stats lm coefficients poly wilcox.test ecdf
#' @importFrom methods is
#' @importFrom MIRA binRegion
#' @importFrom tidyr gather
#' @importFrom grid grid.newpage grid.grabExpr grid.draw popViewport 
#'             pushViewport unit viewport
#' @importFrom grDevices dev.off
NULL


# Because of some issues with CRAN submission, 
# (see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation, 
# in order to pass some R CMD check NOTES.
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        ".", "bin", "binID", "chr", "element_blank", "featureID", 
        "geom_violin", "methylCount", "id", "meth", 
        "methylProp", "coverage", "sumCoverage", "regionGroupID", "regionID", 
        "sampleName", "sampleType", "theme", "ubinID", "V1", 
        "mean_region_size", "region_coverage", "rowIndex", "rsIndex",
        "total_region_number", "cytosine_coverage", ".SD")) 
}

#########################################################################


#' Function to aggregate PCA loading weights over a given region set
#' and then get p value for each PC based on a permutation
#' 
#' @param loadingMat matrix of loadings (the coefficients of 
#' the linear combination that defines each PC). One named column for each PC.
#' One row for each original dimension/variable (should be same order 
#' as original data/mCoord). The x$rotation output of prcomp().
#' @param mCoord a GRanges object or data frame with coordinates 
#' for the cytosines included in the PCA. Coordinates should be in the 
#' same order as the methylation data and loadings. If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.
#' @param regionSet A genomic ranges object with regions corresponding
#' to the same biological annotation.
#' @param PCsToAnnotate A character vector with principal components to 
#' include. eg c("PC1", "PC2")
# #UPDATE: make sure only aggregating PCsToAnnotate to save time
#' @param metric Scoring metric. "raw" is a simple
#' average of the loading values with no normalization (recommended). With
#' "raw" score, be cautious in interpretation for
#'  region sets with low number of regions.
#' Mean difference ("meanDiff") is also supported but is skewed toward 
#' ranking large region sets highly.
#' @param pcLoadAv The average absolute loading value for each PC. Will
#' significantly speed up computation if this is given.
#' @param verbose A boolean. Whether progress 
#' of the function should be shown, one
#' bar indicates the region set is completed. Useful when using 
#' aggregateLoadings with 'apply' to do many region sets at a time.

#' @return a data.table with one row and the following 
#' columns: one column for each item of PCsToAnnotate with name given
#' by the item which have scores for each PC, cytosine_coverage which
#' has number of cytosines that overlapped with regionSet,
#' region_coverage which has number of regions from regionSet
#' that overlapped any cytosines, total_region_number that has
#' number of regions in regionSet, mean_region_size that has average
#' size in base pairs of regions in regionSet, the average is based on
#' all regions in regionSet and not just ones that overlap. 

aggregateLoadings <- function(loadingMat, mCoord, regionSet, 
                  PCsToAnnotate = c("PC1", "PC2"), metric="raw", pcLoadAv=NULL,
                  verbose=FALSE) {
    # UPDATE: take out permutations
    permute <- FALSE
    
    if (is(mCoord, "GRanges")) {
        coordinateDT <- grToDt(mCoord)
    } else if (is(mCoord, "data.frame")) {
        coordinateDT <- mCoord
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    
    numOfRegions <- length(regionSet)
    totalCpGs <- nrow(loadingMat)
    
    # extreme positive or negative values both give important information
    loadingMat <- abs(loadingMat) # required for later code
    loadingDT <- as.data.table(loadingMat)
    
    # reformat into data.table with chromosome location and weight
    # also restricting to PCsToAnnotate so unnecessary computations 
    # are not done
    loadingDT <- data.table(coordinateDT, 
                            loadingDT[, PCsToAnnotate, with=FALSE])
    # naming does not work if only using one PC so add this line for that case
    setnames(loadingDT, c("chr", "start", PCsToAnnotate))
    
    # would rounding speed up aggregation?, potentially make a sparse matrix
    # if a lot of entries became 0
    
    # specify aggregation operation
    # will be done separately for each PC specified
    aggrCommand <- buildJ(PCsToAnnotate, 
                                rep("mean", length(PCsToAnnotate)))
    
    # do the actual aggregation
    if (metric == "raw") {
        # previously used BSAggregate from RGenomeUtils but now using local, 
        # modified copy
        loadAgMain <- BSAggregate(BSDT = loadingDT, 
                                  regionsGRL = GRangesList(regionSet),
                                 jExpr = aggrCommand,
                                 byRegionGroup = TRUE,
                                 splitFactor = NULL)
        results <- loadAgMain[, .SD, .SDcols = PCsToAnnotate]
        # if no cytosines from loadings were included in regionSet, result is NA
        if (is.null(loadAgMain)) {
            results <- as.data.table(t(rep(NA, length(PCsToAnnotate))))
            setnames(results, PCsToAnnotate)
            results[, cytosine_coverage := 0]
            results[, region_coverage := 0]
            results[, total_region_number := numOfRegions]
            results[, mean_region_size := round(mean(width(regionSet)), 1)]
        } else {
            results[, cytosine_coverage := loadAgMain[, .SD, .SDcols = "numCpGsOverlapping"]]
            results[, region_coverage := loadAgMain[, .SD, .SDcols = "numRegionsOverlapping"]]
            results[, total_region_number := numOfRegions]
            results[, mean_region_size := round(mean(width(regionSet)), 1)]
        }
    } else if (metric == "raw_CpG") {
        
        # average of loadings for all CpGs within region set
        loadMetrics <- cpgOLMetrics(dataDT=loadingDT, regionGR=regionSet, 
                                   metrics="mean", 
                                   alsoNonOLMet=FALSE)
        if (is.null(loadMetrics)) {
            results <- as.data.table(t(rep(NA, length(PCsToAnnotate))))
            setnames(results, PCsToAnnotate)
            results[, cytosine_coverage := 0]
            results[, region_coverage := 0]
            results[, total_region_number := numOfRegions]
            results[, mean_region_size := round(mean(width(regionSet)), 1)]
        } else {
            
            # simple mean 
            results <- as.data.table(t(loadMetrics$mean_OL))
            colnames(results) <- loadMetrics$testCol
            
            # add information about degree of overlap
            results <- cbind(results, 
                             loadMetrics[1, .SD, 
                                         .SDcols = c("cytosine_coverage", 
                                                     "region_coverage", 
                                                     "total_region_number", 
                                                     "mean_region_size")]) 
        }
        
        
    } else if (metric == "meanDiff") {
        # if (is.null(pcLoadAv)) {
        #     # calculate (should already be absolute)
        #     pcLoadAv <- apply(X = loadingDT[, PCsToAnnotate, with=FALSE], MARGIN = 2, FUN = mean)
        # }
        loadMetrics <- cpgOLMetrics(dataDT=loadingDT, regionGR=regionSet, 
                                   metrics=c("mean", "sd"), 
                                   alsoNonOLMet=TRUE)
        if (is.null(loadMetrics)) {
            results <- as.data.table(t(rep(NA, length(PCsToAnnotate))))
            setnames(results, PCsToAnnotate)
            results[, cytosine_coverage := 0]
            results[, region_coverage := 0]
            results[, total_region_number := numOfRegions]
            results[, mean_region_size := round(mean(width(regionSet)), 1)]
        } else {
            # calculate mean difference
            # pooled standard deviation
            sdPool <- sqrt((loadMetrics$sd_OL^2 + loadMetrics$sd_nonOL^2) / 2)
            
            # mean difference
            # error if numCpGsOverlapping > (1/2) * totalCpGs
            meanDiff <- (loadMetrics$mean_OL - loadMetrics$mean_nonOL) / 
                (sdPool * sqrt((1 / loadMetrics$cytosine_coverage) - (1 / (totalCpGs - loadMetrics$cytosine_coverage))))
            results <- as.data.table(t(meanDiff))
            colnames(results) <- loadMetrics$testCol
            
            # add information about degree of overlap
            results <- cbind(results, loadMetrics[1, .SD, .SDcols = c("cytosine_coverage", "region_coverage", 
                                                      "total_region_number", "mean_region_size")]) 
         }
        
        
    } else if (metric == "rankSum"){
        # one sided test since I took the absolute value of the loadings 
        wRes <- rsWilcox(dataDT = loadingDT, regionGR=regionSet, alternative="greater")
        
        if (is.null(wRes)) {
            results <- as.data.table(t(rep(NA, length(PCsToAnnotate))))
            setnames(results, PCsToAnnotate)
            results[, cytosine_coverage := 0]
            results[, region_coverage := 0]
            results[, total_region_number := numOfRegions]
            results[, mean_region_size := round(mean(width(regionSet)), 1)]
        } else {
            results <- as.data.table(wRes)
        }
    } else {
        stop(cleanws("metric was not recognized. 
                      Check spelling and available options."))
    }
    
    
    # cpgOLMetrics() # make sure it works with no overlap
    if (verbose) {
        message("|")
    }
        
        return(results)
    }



#' Do PCRSA for many region sets
#'
#' @param loadingMat matrix of loadings (the coefficients of 
#' the linear combination that defines each PC). One named column for each PC.
#' One row for each original dimension/variable (should be same order 
#' as original data/mCoord). The x$rotation output of prcomp().
#' @param mCoord a GRanges object or data frame with coordinates 
#' for the cytosines included in the PCA. Coordinates should be in the 
#' same order as the methylation data and loadings. If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.
#' @param GRList GRangesList object. Each list item is 
#' a distinct region set to test (region set: regions that correspond to 
#' the same biological annotation).
#' @param PCsToAnnotate A character vector with principal components to 
#' include. eg c("PC1", "PC2")
#' @param scoringMetric Scoring metric. "raw" is a simple
#' average of the loading values with no normalization (recommended). With
#' "raw" score, be cautious in interpretation for
#'  region sets with low number of regions.
#' Mean difference ("meanDiff") is also supported but is skewed toward 
#' ranking large region sets highly.
#' @param verbose A boolean. Whether progress 
#' of the function should be shown, one
#' bar indicates the region set is completed.
#' @return data.frame of results, one row for each region set. 
#' One column for each PC in PCsToAnnotate
#' with score for that PC for a given region set (specific score depends
#' on "scoringMetric" parameter). 
#' Rows will be in the same order as region sets in GRList 
#' and will be named if GRList was named
#' "cytosine_coverage" column has number of cytosines that 
#' overlapped with the given region set.
#' "region_coverage" column has number of regions 
#' that overlapped with any cytosines.
#' "total_region_number" column has total number of regions. 
#' "mean_region_size" has average region size (average of all regions,
#' not just those that overlap a cytosine).
#' 
#' @examples 
#' data("brcaCoord1")
#' data("brcaLoadings1")
#' data("esr1_chr1")
#' data("nrf1_chr1")
#' GRList = GenomicRanges::GRangesList(esr1_chr1, nrf1_chr1)
#' pcRegionSetEnrichment(loadingMat=brcaLoadings1, mCoord=brcaCoord1, 
#'                       GRList=GRList, PCsToAnnotate=c("PC1", "PC2"), 
#'                       scoringMetric="raw")

#' 
#' @export

pcRegionSetEnrichment <- function(loadingMat, mCoord, GRList, 
                  PCsToAnnotate = c("PC1", "PC2"), scoringMetric = "raw",
                  verbose=TRUE) {
    
    if (!(any(c(is(loadingMat, "matrix"), is(loadingMat, "data.frame"))))) {
        warning("loadingMat should be a matrix or data.frame.")
    }
 
    if (is(mCoord, "GRanges")) {
        coordinateDT <- grToDt(mCoord)
    } else if (is(mCoord, "data.frame")) {
        coordinateDT <- mCoord
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    # should be GRangesList
    if (is(GRList, "GRanges")) {
        GRList <- GRangesList(GRList)
    } else if (!is(GRList, "GRangesList")) {
        stop("GRList should be a GRangesList object.")
    }
    
    
    # apply over the list of region sets
    resultsList <- lapplyAlias(GRList, 
                              function(x) aggregateLoadings(loadingMat=loadingMat, 
                                                            mCoord=coordinateDT, 
                                                            regionSet=x, 
                                                            PCsToAnnotate = PCsToAnnotate,
                                                            metric = scoringMetric, 
                                                            verbose=verbose))
    resultsDT <- do.call(rbind, resultsList) 
    row.names(resultsDT) <- row.names(GRList)
    
    resultsDF <- as.data.frame(resultsDT)
    return(resultsDF)
}    



#' Function to create "meta-region" loading profile that indicates enrichment 
#' of cytosines with high loading values in region set but not in
#' surrounding genome. All regions in a given region set 
#' are combined into a single aggregate profile.
#'
#'
#' @param loadingMat matrix of loadings (the coefficients of 
#' the linear combination that defines each PC). One named column for each PC.
#' One row for each original dimension/variable (should be same order 
#' as original data/mCoord). The x$rotation output of prcomp().
#' @param mCoord a GRanges object or data frame with coordinates 
#' for the cytosines included in the PCA. Coordinates should be in the 
#' same order as the methylation data and loadings. If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.
#' @param GRList GRangesList object. Each list item is 
#' a distinct region set (regions that correspond to 
#' the same biological annotation).
#' @param PCsToAnnotate A character vector with principal components to 
#' include. eg c("PC1", "PC2")
#' @param binNum Number of bins to split the regions into when
#' making the aggregate loading profile. More bins will
#' give a higher resolution but perhaps more noisy profile.
#' @param verbose A boolean. Whether progress 
#' of the function should be shown, one
#' bar indicates the region set is completed.
#' @return A list of data.tables each data.table for
#' a separate region set. The data table has the binned loading profile,
#' one row per bin.
#' 
#' @examples 
#' data("brcaCoord1")
#' data("brcaLoadings1")
#' data("esr1_chr1")
#' data("nrf1_chr1")
#' GRList = GenomicRanges::GRangesList(esr1_chr1, nrf1_chr1)
#' pcEnrichmentProfile(loadingMat=brcaLoadings1, 
#'                     mCoord=brcaCoord1, 
#'                     GRList=GRList, 
#'                     PCsToAnnotate=c("PC1", "PC2"), 
#'                     binNum=25)
#' @export

pcEnrichmentProfile = function(loadingMat, mCoord, GRList,
                    PCsToAnnotate = c("PC1", "PC2"), binNum = 25,
                    verbose=TRUE) {
    
    if (!(any(c(is(loadingMat, "matrix"), is(loadingMat, "data.frame"))))) {
        warning("loadingMat should be a matrix or data.frame.")
    }
    
    # checks for correct input
    if (is(mCoord, "GRanges")) {
        coordinateDT <- grToDt(mCoord)
    } else if (is(mCoord, "data.frame")) {
        coordinateDT <- mCoord
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    # should be GRangesList
    if (is(GRList, "GRanges")) {
        GRList <- GRangesList(GRList)
    } else if (is(GRList, "list")) {
        #making sure each list item is a GRanges object
        if (all(sapply(X = GRList, function(x) is(x, "GRanges")))) {
            GRList <- GRangesList(GRList)
        } else {
            stop("GRList should be a GRanges or GRangesList object.")
        }
    } else if (!is(GRList, "GRangesList")) {
        stop("GRList should be a GRanges or GRangesList object.")
    } 
    
    
    loadingDT <- as.data.table(abs(loadingMat))
    loadingDT <- cbind(coordinateDT, loadingDT)
    
    GRDTList <- lapply(X = GRList, FUN = grToDt)
    profileList <- lapply(X = GRDTList, function(x) BSBinAggregate(BSDT = loadingDT, 
                                                    rangeDT = x, 
                                                    binCount = binNum, 
                                                    minReads = 0, 
                                                    byRegionGroup = TRUE, 
                                                    splitFactor = NULL,
                                                    PCsToAnnotate = PCsToAnnotate))
    
    profileList <- lapply(profileList, FUN = makeSymmetric)
    profileList <- lapply(profileList, function(x) x[, regionGroupID := 1:binNum])
    
    return(profileList)
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
# @param BSDT A data.table. For PCRSA, a data.table of loading values
# with the PCs to be annotated. One column for the loadings of each PC
# and also has columns with the coordinates for CpGs that the loadings
# are for: chr (chromosome) and start column
# @param rangeDT A data table with the sets of regions to be binned, 
# with columns named start, end
# @param binCount Number of bins across the region
# @param byRegionGroup Pass along to binCount (see ?binCount)
# @param minReads Filter out bins with fewer than X reads before returning.
# @param PCsToAnnotate A character vector with principal components to 
# analyze. eg c("PC1", "PC2")
# @param verbose A boolean. Whether progress 
# of the function should be shown, one
# bar indicates the region set is completed.
# useful when using BSBinAggregate with 'apply' to do many 
# region sets at a time.
BSBinAggregate <- function(BSDT, rangeDT, binCount, minReads = 500, 
                          byRegionGroup = TRUE, 
                          splitFactor = NULL,
                          PCsToAnnotate, verbose=FALSE) {
    if (! "data.table" %in% class(rangeDT)) {
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
    binnedDT <- rangeDT[, MIRA::binRegion(start, end, binCount, get(seqnamesColName))]
    binnedGR <- sapply(split(binnedDT, binnedDT$binID), dtToGr)
    # message("Aggregating...")
    
    # RGenomeUtils::BSAggregate
    binnedBSDT <- BSAggregate_RGenomeUtils(BSDT, 
                             regionsGRL=GRangesList(binnedGR), 
                             jExpr=buildJ(PCsToAnnotate, rep("mean", length(PCsToAnnotate))), 
                                          byRegionGroup=byRegionGroup, 
                                          splitFactor = splitFactor)
    # # If we aren't aggregating by bin, then don't restrict to min reads!
    # if (byRegionGroup) {
    #     binnedBSDT <- binnedBSDT[readCount > minReads,]
    # }
    if (verbose) {
        message("|")
    }
    
    return(binnedBSDT)
}

#' modification of BSAggregate to just return mean per region
#' 
#' @param loadingMat matrix of loadings (the coefficients of 
#' the linear combination that defines each PC). One named column for each PC.
#' One row for each original dimension/variable (should be same order 
#' as original data/mCoord). The x$rotation output of prcomp().
#' @param mCoord a GRanges object or data frame with coordinates 
#' for the cytosines included in the PCA. Coordinates should be in the 
#' same order as the methylation data and loadings. If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.
#' @param regionSet A GRanges object with regions corresponding
#' to the same biological annotation.
#' @param PCsToAnnotate A character vector with principal components to 
#' include. eg c("PC1", "PC2")
#' @param returnQuantile Boolean. If FALSE, return region averages. If TRUE,
#' for each region, return the quantile of that region's average value
#' based on the distribution of individual cytosine values
#' @return a data.table with region coordinates and average loading 
#' values for each region. Has columns chr, start, end, and a column for each
#' PC in PCsToAnnotate. Regions are not in order along the rows of the data.table.
#'
# @example averageByRegion(BSDT = BSDT, regionsGRL, 
#          jCommand = MIRA:::buildJ(cols = "methylProp", "mean")) 
# Devel note: I could add a column for how many cytosines are in each region 

averageByRegion <- function(loadingMat, 
                            mCoord, 
                            regionSet, 
                            PCsToAnnotate = c("PC1", "PC2"), 
                            returnQuantile=FALSE) {
    
    
    
    if (is(mCoord, "GRanges")) {
        coordinateDT <- grToDt(mCoord)
    } else if (is(mCoord, "data.frame")) {
        coordinateDT <- mCoord
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    
    # old parameters of BSAggregate
    excludeGR <- NULL
    splitFactor <- NULL
    regionsGRL.length <- NULL
    keep.na <- FALSE
    keepCols <- NULL
    sumCols <- NULL
    byRegionGroup <- FALSE
    
    # reformating parameters and assigning new names
    GRList <- regionSet
    regionsGRL <- GRList
    loadingDT <- as.data.table(abs(loadingMat))
    # linking coordinates to loading values, has columns chr start, PCsToAnnotate
    BSDT <- cbind(coordinateDT, loadingDT[, .SD, .SDcols = PCsToAnnotate])
    jCommand <- buildJ(PCsToAnnotate, rep("mean", length(PCsToAnnotate)))
    
    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if (is(regionsGRL, "GRanges")) {
        regionsGRL <- GRangesList(regionsGRL);
    } else if (!is(regionsGRL, "GRangesList")) {
        stop("regionsGRL is not a GRanges or GRangesList object");
    }
    
    if (! is.null(excludeGR)) {
        BSDT <- BSFilter(BSDT, minReads = 0, excludeGR)
    }
    
    bsgr <- BSdtToGRanges(list(BSDT));
    
    colModes <- sapply(BSDT, mode);
    if (is.null(sumCols)) {
        sumCols <- setdiff(colnames(BSDT), c("chr", "start", "end", 
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
    
    if (is.null(regionsGRL.length)) {
        # if (length(regionsGRL) > 100) {
        #     message(cleanws("BSAggregate: Calculating sizes. You can speed this
        #                     up by supplying a regionsGRL.length vector..."),
        #             appendLF = FALSE)
        # }
        regionsGRL.length <- sapply(regionsGRL, length)
        # message("Done counting regionsGRL lengths.");
    }
    
    # Build a table to keep track of which regions belong to which group
    region2group <- data.table(
        regionID = seq_along(regionsGR), 
        chr = as.vector(seqnames(regionsGR)), 
        start = as.vector(start(regionsGR)), 
        end = as.vector(end(regionsGR)), 
        withinGroupID = as.vector(unlist(sapply(regionsGRL.length, seq))), 
        regionGroupID = rep(seq_along(regionsGRL), regionsGRL.length))
    setkey(region2group, regionID)
    
    
    # message("Finding overlaps...");
    fo <- findOverlaps(bsgr[[1]], regionsGR)
    
    setkey(BSDT, chr, start)
    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$coverage == BSDT$coverage))
    
    # message("Setting regionIDs...");
    BSDT <- BSDT[queryHits(fo), ] # restrict the table to CpGs in any region.
    
    if (NROW(BSDT) < 1) {
        warning("No BSDT sites in the given region list. 
                Please expand your regionsGRL")
        return(NULL)
    }
    
    # record which region they overlapped, corresponds to withinGroupID
    BSDT[, regionID := subjectHits(fo)]
    # if (!keep.na) {
    # BSDT <- BSDT[queryHits(fo), ]
    #}
    
    if (is.null(jCommand)) {
        cols <- c(sumCols, keepCols)
        funcs <- c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
        jCommand <- buildJ(cols, funcs)
    }
    # message("jCommand: ", jCommand)
    
    # Build the by string
    if (is.null(splitFactor)) {
        byString <- paste0("list(regionID)");
    } else {
        byString <- paste0("list(", paste("regionID", paste0(splitFactor, ""), 
                                          collapse = ",", sep = ","), ")")
    }
    
    # Now actually do the aggregate:
    # message("Combining...");
    # for MIRA: average methylProp and sum coverage within each instance of 
    # each bin (ie average methylProp's in bin1 of first region, average 
    # methylProp's in bin1 of second region etc. for all bins and all regions
    # separately)
    bsCombined <- BSDT[, eval(parse(text = jCommand)), 
                       by = eval(parse(text = byString))]
    setkey(bsCombined, regionID)
    setkey(region2group, regionID)
    
    avPerRegion <- merge(bsCombined, region2group)
    avPerRegion[, c("regionID", "withinGroupID", "regionGroupID") := NULL]
    
    if (returnQuantile) {
       
        for (i in seq_along(PCsToAnnotate)) {
            # perhaps this could be more efficient with mapply
            avPerRegion[, c(PCsToAnnotate[i]) := ecdf(loadingDT[[PCsToAnnotate[i]]])(avPerRegion[[PCsToAnnotate[i]]])]
            # I tried set() to improve performance but it took about the same time
        }

    }
    
    
    return(avPerRegion)
    
    # Now aggregate across groups.
    # I do this in 2 steps to avoid assigning regions to groups, 
    # which takes awhile. I think this preserves memory and is faster.
    
    # # Define aggregation column. aggregate by region or by region group?
    # if (byRegionGroup) {
    #     # must set allow = TRUE here in case there are multiple IDs (splitCol)
    #     # adds regionGroupID column from region2group to bsCombined
    #     bsCombined[region2group, regionGroupID := regionGroupID, allow = TRUE]
    #     if (! is.null(splitFactor)) { 
    #         byStringGroup <- paste0("list(", 
    #                                 paste("regionGroupID", 
    #                                       paste0(splitFactor, collapse = ","), 
    #                                       sep = ","), 
    #                                 ")")
    #     } else {
    #         byStringGroup <- "list(regionGroupID)"
    #     }
    #     
    #     
    #     # actual aggregation operation
    #     # for normal MIRA use: averaging methylProp's and summing
    #     # coverage for all bins with the same number (ie all
    #     # bin1's together, all bin2's together, etc.)
    #     # NOTE: 2nd use of the jCommand so if the first jCommand use changed
    #     # column names that are required by the jCommand 
    #     # then this 2nd jCommand use will cause an error
    #     bsCombined <- bsCombined[, eval(parse(text = jCommand)), 
    #                              by = eval(parse(text = byStringGroup))]
    #     
    #     # if any strand information was not given, averaging the profiles 
    #     # about the center to account for unknown strand orientation, 
    #     # also averaging coverage about center
    #     # ie if any "*" are present then average
    #     if ("*" %in% unique(as.character(strand(regionsGR)))) {
    #         bsCombined[, methylProp := (methylProp + rev(methylProp)) / 2]
    #         if (hasCoverage) {
    #             bsCombined[, coverage := (coverage + rev(coverage)) / 2]
    #         }
    #     }
    #     
    #     # changing "regionGroupID" name to "bin" which is less confusing
    #     # for normal MIRA use cases
    #     setnames(bsCombined, old = "regionGroupID", new = "bin")
    #     return(bsCombined[]);
    # } else {
    #     warning(cleanws("Using byRegionGroup = FALSE may 
    #          result in missing functionalities such as symmetrical averaging"))
    #     e <- region2group[bsCombined, ]
    #     setkey(e, regionID);
    #     return(e);
    # }
    # # WARNING: There are now 2^2 ways to aggregate, sum vs mean
    # # at each level: across regions, then across region sets. THis
    # # doesn't give you a choice at this point. 
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
# GenomicRanges::findOverlaps(regionSet, dtToGr(coordinateDT)) 

###########################################################################

#  BSAggregate from RGenomeUtils
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
# @param excludeGR A GenomicRanges object with regions you want to 
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
# and also for rows that do not overlap with region set.
#
# 
BSAggregate <- function(BSDT, regionsGRL, excludeGR=NULL, regionsGRL.length = NULL, splitFactor=NULL, keepCols=NULL, 
                       sumCols=NULL, jExpr=NULL, byRegionGroup=FALSE, keep.na=FALSE, returnSD=FALSE) {
    
    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if(is(regionsGRL, "GRanges")) {
        regionsGRL <- GRangesList(regionsGRL)
    } else if (!is(regionsGRL, "GRangesList")) {
        stop("regionsGRL is not a GRanges or GRangesList object")
    }
    
    if(! is.null(excludeGR)) {
        BSDT <- BSFilter(BSDT, minReads=0, excludeGR)
    }
    
    bsgr <- BSdtToGRanges(list(BSDT))
    
    additionalColNames <- setdiff(colnames(BSDT), c("chr","start", "end","hitCount","readCount", splitFactor))
    
    colModes <- sapply(BSDT,mode)
    if (is.null(sumCols)) {
        sumCols <- setdiff(colnames(BSDT),c("chr", "start", "end", "strand", splitFactor, keepCols))
        # Restrict to numeric columns.		
        sumCols <- intersect(sumCols, names(colModes[which(colModes == "numeric")]))
        
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
        regionsGRL.length <- sapply(regionsGRL, length)
        # message("Done counting regionsGRL lengths.")
    }
    
    # Build a table to keep track of which regions belong to which group
    region2group <- data.table(
        regionID=1:length(regionsGR), 
        chr=as.vector(seqnames(regionsGR)), 
        start=as.vector(start(regionsGR)), 
        end=as.vector(end(regionsGR)),
        withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
        regionGroupID=rep(1:length(regionsGRL), regionsGRL.length))
    setkey(region2group, regionID)
    
    
    # message("Finding overlaps...")
    fo <- findOverlaps(query = bsgr[[1]], subject = regionsGR)
    
    ### use info from findOverlaps to see how many individual
    # cytosines and how many regions overlap with each other
    numCpGsOverlapping <- length(unique(queryHits(fo)))
    numRegionsOverlapping <- length(unique(subjectHits(fo)))
    
    
    setkey(BSDT, chr, start)
    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))
    
    # message("Setting regionIDs...")
    BSDT <- BSDT[queryHits(fo),] #restrict the table to CpGs in any region.
    
    if (NROW(BSDT) < 1) {
        warning("No BSDT sites in the given region list; please expand your regionsGRL")
        return(NULL)
    }
    
    BSDT[,regionID:=subjectHits(fo)] #record which region they overlapped.
    #BSDT[queryHits(fo),regionID:=subjectHits(fo)]
    #if (!keep.na) {
    #	BSDT = BSDT[queryHits(fo),]
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
        byString <- paste0("list(", paste("regionID", paste0(splitFactor, ""), collapse=", ", sep=", "), ")")
    }
    
    # Now actually do the aggregate:
    # message("Combining...")
    bsCombined <- BSDT[,eval(parse(text=jExpr)), by=eval(parse(text=byString))]
    setkey(bsCombined, regionID)
    # Now aggregate across groups.
    # I do this in 2 steps to avoid assigning regions to groups,
    # which takes awhile. I think this preserve memory and is faster.
    
    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        # must set allow=TRUE here in case there are multiple IDs (splitCol)
        bsCombined[region2group, regionGroupID:=regionGroupID, allow=TRUE]
        if (! is.null(splitFactor) ) { 
            byStringGroup <- paste0("list(", paste("regionGroupID", paste0(splitFactor, collapse=", "), sep=", "), ")")
        } else {
            byStringGroup <- "list(regionGroupID)"
        }
        bsCombined=bsCombined[,eval(parse(text=jExpr)), by=eval(parse(text=byStringGroup))]
        bsCombined[, numCpGsOverlapping := numCpGsOverlapping]
        bsCombined[, numRegionsOverlapping := numRegionsOverlapping]
        return(bsCombined)
    } else {
        e <- region2group[bsCombined,]
        setkey(e, regionID)
        return(e)
    }
    # WARNING: There are now 2^2 ways to aggregate, sum vs mean
    # at each level: across regions, then across region sets. THis
    # doesn't give you a choice at this point. 
}


 


#' Function to see how representative the cytosines in a region set are:
#' can cytosines from a region set reproduce PC ordering and how high is 
#' the correlation of original PCs with ordering from subset of cytosines?
#' @param regionSet The region set to subset by. Only loading values for
#' cytosines in these regions will be used to make a PC score.
#' @param mPCA The pca data. Output from prcomp() of DNA methylation data. 
#' Should have 
#' mPCA$rotation with PCA loadings, mPCA$x with PC scores and 
#' mPCA$center with values used to center DNA methylation.
#' @param methylData DNA methylation levels in matrix or data.frame. 
#' Rows are cytosines. Columns are 
#' samples.
#' @param mCoord a GRanges object or data frame with coordinates 
#' for the cytosines included in the PCA. Coordinates should be in the 
#' same order as the methylation data and loadings. If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.
#' @param pc A character object with PCs of interest (eg "PC1"). Loadings from
#' this PC will be used to generate new PC scores based on only loadings
#' that overlap with regionSet. Can be multiple PCs.
#' @param returnCor Option to return correlation between these scores and
#' original PC scores. 
#' 
#' @return a score for each patient from only loading values for cytosines in 
#' regionSet. If returnCor = TRUE, returns a correlation score for scores from
#' each "pc" with scores from a subset of loading values for that PC.
#' 

pcFromSubset <- function(regionSet, mPCA, methylData, mCoord, 
                         pc="PC1", returnCor=FALSE) {
    
    # test for appropriateness of inputs/right format
    if (is(mCoord, "GRanges")) {
        coordGR <- mCoord
    } else if (is(mCoord, "data.frame")) {
        # UPDATE: does the work on data.frames that are not data.tables?
        coordGR <- dtToGr(mCoord)
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    
    # get subset of loading values
   #  subsetInd = 
    newColNames <- paste0(PCofInterest, "_subset")
    
    
    olList <- findOverlaps(query = regionSet, subject = coordGR)
    # regionHitInd <- sort(unique(queryHits(olList)))
    cytosineHitInd <- sort(unique(subjectHits(olList)))
    thisRSMData <- t(methylData[cytosineHitInd, ])
    # subject_ID <- row.names(thisRSMData)
    # centeredPCAMeth <- t(apply(t(methylData), 1, function(x) x - pcaData$center)) # center first 
    # reducedValsPCA <- centeredPCAMeth %*% pcaData$rotation
    # reducedValsPCA <- pcaData$x
    # use subset of loading values on subset of DNA methylation values
    # getting PC scores manually so artificial PCs will be included (PC1m4 and PC1p3)
    centeringSubset <- mPCA$center[cytosineHitInd]
    centeredMeth <- t(apply(thisRSMData, 1, function(x) x - centeringSubset)) # center first 
    reducedValsPCA <- centeredMeth %*% mPCA$rotation[cytosineHitInd, PCofInterest]
    colnames(reducedValsPCA) <- newColNames
    # pcaValDF <- as.data.frame(reducedValsPCA)
    
    
    # optional calculate correlation
    if (returnCor) {
        # plot(reducedValsPCA[, PCofInterest],mPCA$x[, PCofInterest])
        corMat <- cor(reducedValsPCA[, newColNames],mPCA$x[, PCofInterest])
        # only correlation between subset and matching PC, not with other PCs
        mainCor <- diag(corMat)
        names(mainCor) <- PCofInterest
        return(mainCor)
        # origPCCor <- cor(pca$x[, PCofInterest])
        # subsetPCCor <- cor(reducedValsPCA[, PCofInterest])
        # Heatmap(corMat, cluster_columns = FALSE, cluster_rows = FALSE)
        # Heatmap(origPCCor, cluster_columns = FALSE, cluster_rows = FALSE)
        # Heatmap(subsetPCCor, cluster_columns = FALSE, cluster_rows = FALSE)
        # Heatmap(pca$rotation[cytosineHitInd, PCofInterest])
    } else {
        colnames(reducedValsPCA) <- PCofInterest
        return(reducedValsPCA)
    }
    
}

#' Determine what percent of cytosines in region sets overlap.
#' Slightly different than percent overlap of region sets.
#' @param mCoord a GRanges object or data frame with coordinates 
#' for the cytosines included in the PCA. Coordinates should be in the 
#' same order as the methylation data and loadings. If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.
#' @param GRList GRangesList object. Each list item is 
#' a distinct region set (regions that correspond to 
#' the same biological annotation).
#' @return Total number of cytosines from methylation data that overlapped with
#' each region set. Also the percent of cytosines that overlapped with
#' other region set.
#'

percentCOverlap <- function(mCoord, GRList) {
    
    # test for appropriateness of inputs/right format
    if (is(mCoord, "GRanges")) {
        coordGR <- mCoord
    } else if (is(mCoord, "data.frame")) {
        # UPDATE: does the work on data.frames that are not data.tables?
        coordGR <- dtToGr(mCoord)
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    
    # overlap between each region set and methylation data to see 
    # which cytosines were covered by each region set
    olList <- lapply(X = GRList, 
                    function(x) findOverlaps(query = x, subject = coordGR))
    
    # indices for cytosines that overlapped with each region set
    cHitInd <- lapply(X = olList, function(x) sort(unique(subjectHits(x))))

    # total number
    total_cytosines <- sapply(X = cHitInd, length)
    
    # make matrix of combinatorial intersection, will eventually have percent overlap
    pOL <- matrix(nrow=length(GRList), ncol=length(GRList))
    
    # intersect each iten in cHitInd with each item in cHitInd
    sharedC <- lapply(X = cHitInd, 
                     FUN = function(x) sapply(X = cHitInd, FUN = function(y) length(intersect(x, y))))
    
    # converting to proportion of cytosines covered by each region set
    prop_overlap <- lapply(X = sharedC, function(x) x / total_cytosines)
    prop_overlap <- do.call(rbind, prop_overlap)

    return(list(prop_overlap, total_cytosines))    
    
}

#' A convenience function to easily get indices for
#' the top region sets that were enriched for each PC. 
#' For each PC, get index of original region sets but ordered by rsScores
#' ranking for each PC. First number in a given column will be 
#' original index of the region set ranked first for that PC. Second row for a
#' column will be the original index of the region set that ranked second
#' for that PC, etc. 
#' 
#' Use this function when you want to look at top region sets to make it 
#' easier to get the original indices to select them from a list of region sets.
#' 
#' @param rsScores a data.table with scores for each 
#' region set from main PCRSA function. 
#' Each row is a region set. Columns are PCs and info on region set overlap
#' with DNA methylation data. Should be in the same order as GRList (the list of 
#' region sets used to create it.)
#' @param PCsToAnnotate PCs in rsScores for which you want
#' the indices of the original region sets
#' @return A data.frame with columns PCsToAnnotate. Each column has been 
#' ranked by enrichment score for region sets for that PC.
#' Original indices for region sets that were used to create rsScores
#' are given. 
#' @export
#' 
rsRankingIndex <- function(rsScores, PCsToAnnotate) {
    rsScores <- as.data.table(rsScores) # UPDATE to not be data.table
    # so by references changes will not be a problem
    rsScores <- copy(rsScores)
    rsScores[, rsIndex := 1:nrow(rsScores)]
    
    PCsToAnnotate <- PCsToAnnotate[PCsToAnnotate %in% colnames(rsScores)]
    
    rsEnSortedInd <- subset(rsScores, select= PCsToAnnotate)
    
    # then scores by each PC and make a column with the original index for sorted region sets
    # this object will be used to pull out region sets that were top hits for each PC
    for (i in seq_along(PCsToAnnotate)) {
        
        # -1 for decreasing order of scores
        setorderv(rsScores, cols = PCsToAnnotate[i], order=-1L)
        
        rsEnSortedInd[, PCsToAnnotate[i] := rsScores[, rsIndex]]
    }
    
    # reset order
    # setorderv(rsScores, cols = "rsIndex", order=1L)
    return(as.data.frame(rsEnSortedInd))
}

#################### Metric functions ########################################
# scores, metrics, or statistical tests

# Instead of averaging within regions first as BSAggregate does,
# this function does a simple average and standard deviation
#  for all CpGs that overlap
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
# @param regionGR GenomicRanges object. Metrics will be calculated on
# only coordinates within this region set (and optionally separately
# on those outside this region set with alsoNonOLMet parameter)
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
cpgOLMetrics <- function(dataDT, regionGR, metrics=c("mean", "sd"), 
                         alsoNonOLMet=TRUE) {
    
    # convert DT to GR for finding overlaps
    dataGR <- BSdtToGRanges(list(dataDT))[[1]]
    
    OL <- findOverlaps(query = regionGR, subject = dataGR)
    
    # if no overlap, exit
    if (length(OL) == 0) {
        return(NULL)
    }
    
    # get indices for overlapping and non overlapping CpGs
    olCpG <- subjectHits(OL)
    
    # region set info
    total_region_number <- length(regionGR)
    mean_region_size <- round(mean(width(regionGR)), 1)
    
    
    # get info on degree of overlap
    # number of CpGs that overlap
    cytosine_coverage <- length(unique(olCpG))
    # number of regions that overlap
    region_coverage <- length(unique(queryHits(OL)))
    
    
    
    
    nonOLCpG <- (1:nrow(dataDT))[-olCpG]
    
    # gets metrics for all columns except chr, start, end
    testCols <- colnames(dataDT)[!(colnames(dataDT) %in% 
                                       c("chr", "start", "end"))] 
    
    jExpr <- buildJ(cols=rep(testCols, each=length(metrics)), 
                    funcs=rep(metrics, length(testCols)),
                    newColNames = paste0(rep(testCols, 
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
    olResults <- sapply(X = metrics, FUN = function(x) as.numeric(olMetrics[, grepl(pattern = x, colnames(olMetrics))]))
    olResults <- as.data.table(olResults)
    setnames(olResults, old = colnames(olResults), new = paste0(colnames(olResults), "_OL"))
    
    
    nonOLResults <- sapply(X = metrics, FUN = function(x) as.numeric(nonOLMetrics[, grepl(pattern = x, colnames(nonOLMetrics))]))
    nonOLResults <- as.data.table(nonOLResults)
    setnames(nonOLResults, old = colnames(nonOLResults), new = paste0(colnames(nonOLResults), "_nonOL"))
    
    
    
    metricDT <- cbind(data.table(testCol=testCols), olResults, nonOLResults, data.table(cytosine_coverage, region_coverage, total_region_number, mean_region_size))
    
    return(metricDT)
}


#' Wilcoxon rank sum test for a region set
#' @param dataDT a data.table with chr, start, end columns as well
#' as columns to get metrics of eg (PC1, PC2). All columns
#' except chr, start, and end will be considered 
#' columns to get the metrics from so no unnecessary columns should be
#' included.
#' @param regionGR Region set, GenomicRanges object
#' @param ... Additional parameters of wilcox.test function. See ?wilcox.test.
#' For instance specify alternative hypothesis: alternative = "greater".
#' @return A vector with a p value for each column other than chr, start or end. 

# could make this a generic function then apply it across PCs of interest

rsWilcox <- function(dataDT, regionGR, ...) {
    
    # region set info
    total_region_number <- length(regionGR)
    mean_region_size <- round(mean(width(regionGR)), 1)
    

    # columns of interest
    testCols <- colnames(dataDT)[!(colnames(dataDT) %in% c("chr", "start", "end"))]
    
    
    # convert DT to GR for finding overlaps
    dataGR <- BSdtToGRanges(list(dataDT))[[1]]
    
    OL <- findOverlaps(query = regionGR, subject = dataGR)
    
    # if no overlap, exit
    if (length(OL) == 0) {
        return(NULL)
    }
    
    # get indices for overlapping and non overlapping CpGs
    olCpG <- subjectHits(OL)
    nonOLCpG <- (1:nrow(dataDT))[-olCpG]
    
    
    # get info on degree of overlap
    # number of CpGs that overlap
    cytosine_coverage <- length(unique(olCpG))
    # number of regions that overlap
    region_coverage <- length(unique(queryHits(OL)))
    
    
    # calculate Wilcoxon rank sum test for each column
    # additional parameters given with ...
    pVals <- sapply(X = testCols, FUN = function(x) wilcox.test(x = as.numeric(as.matrix(dataDT[olCpG, x, with=FALSE])), 
                                               y = as.numeric(as.matrix(dataDT[nonOLCpG, x, with=FALSE])), ...)$p.value)
     
    
    
    wRes <- data.frame(t(pVals), cytosine_coverage, region_coverage, total_region_number, mean_region_size)
    return(wRes)
}


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
# @param regionsGRL Regions across which you want to aggregate.
# @param excludeGR A GenomicRanges object with regions you want to 
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
#
BSAggregate_RGenomeUtils = function(BSDT, regionsGRL, excludeGR=NULL, regionsGRL.length = NULL, splitFactor=NULL, keepCols=NULL, 
                       sumCols=NULL, jExpr=NULL, byRegionGroup=FALSE, keep.na=FALSE) {
    
    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if (is(regionsGRL, "GRanges")) {
        regionsGRL <- GRangesList(regionsGRL);
    } else if (!is(regionsGRL, "GRangesList")) {
        stop("regionsGRL is not a GRanges or GRangesList object");
    }
    
    if(! is.null(excludeGR)) {
        BSDT <- BSFilter(BSDT, minReads=0, excludeGR)
    }
    
    bsgr <- BSdtToGRanges(list(BSDT))
    
    additionalColNames <- setdiff(colnames(BSDT), c("chr","start", "end","hitCount","readCount", splitFactor))
    
    colModes <- sapply(BSDT,mode)
    if (is.null(sumCols)) {
        sumCols <- setdiff(colnames(BSDT),c("chr", "start", "end", "strand", splitFactor, keepCols))
        # Restrict to numeric columns.		
        sumCols <- intersect(sumCols, names(colModes[which(colModes == "numeric")]))
        
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
        regionsGRL.length <- sapply(regionsGRL, length)
        # message("Done counting regionsGRL lengths.")
    }
    
    # Build a table to keep track of which regions belong to which group
    region2group <- data.table(
        regionID=1:length(regionsGR), 
        chr=as.vector(seqnames(regionsGR)), 
        start=as.vector(start(regionsGR)), 
        end=as.vector(end(regionsGR)),
        withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
        regionGroupID=rep(1:length(regionsGRL), regionsGRL.length))
    setkey(region2group, regionID)
    
    
    # message("Finding overlaps...")
    fo <- findOverlaps(bsgr[[1]], regionsGR)
    
    setkey(BSDT, chr, start)
    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))
    
    # message("Setting regionIDs...")
    BSDT <- BSDT[queryHits(fo),] #restrict the table to CpGs in any region.
    
    if (NROW(BSDT) < 1) {
        warning("No BSDT sites in the given region list; please expand your regionsGRL")
        return(NULL)
    }
    
    BSDT[,regionID:=subjectHits(fo)] #record which region they overlapped.
    #BSDT[queryHits(fo),regionID:=subjectHits(fo)]
    #if (!keep.na) {
    #	BSDT <- BSDT[queryHits(fo),]
    #}
    
    if (is.null(jExpr)) {
        cols <- c(sumCols, keepCols)
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
        byString <- paste0("list(", paste("regionID", paste0(splitFactor, ""), collapse=", ", sep=", "), ")")
    }
    
    # Now actually do the aggregate:
    # message("Combining...")
    bsCombined <- BSDT[,eval(parse(text=jExpr)), by=eval(parse(text=byString))]
    setkey(bsCombined, regionID)
    # Now aggregate across groups.
    # I do this in 2 steps to avoid assigning regions to groups,
    # which takes awhile. I think this preserve memory and is faster.
    
    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        # must set allow=TRUE here in case there are multiple IDs (splitCol)
        bsCombined[region2group, regionGroupID:=regionGroupID, allow=TRUE]
        if (! is.null(splitFactor) ) { 
            byStringGroup <- paste0("list(", paste("regionGroupID", paste0(splitFactor, collapse=", "), sep=", "), ")")
        } else {
            byStringGroup <- "list(regionGroupID)"
        }
        bsCombined=bsCombined[,eval(parse(text=jExpr)), by=eval(parse(text=byStringGroup))]
        return(bsCombined)
    } else {
        e <- region2group[bsCombined,]
        setkey(e, regionID)
        return(e)
    }
    # WARNING: There are now 2^2 ways to aggregate, sum vs mean
    # at each level: across regions, then across region sets. THis
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