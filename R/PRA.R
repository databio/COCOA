# package to annotate PCA components of DNA methylation data
# based on region set enrichment
library(RGenomeUtils)

#' Function to aggregate PCA loading weights over a given region set
#' and then get p value for each PC based on a permutation
#' 
#' @param loadingMat 
#' @param regionSet A genomic ranges object with regions corresponding
#' to the same biological annotation.
#' #UPDATE: make sure only aggregating PCsToAnnotate to save time

aggregateLoadings <- function(loadingMat, coordinateDT, regionSet, 
                  PCsToAnnotate = c("PC1", "PC2"), permute=TRUE) {
    # extreme positive or negative values both give important information
    loadingMat = abs(loadingMat) 
    
    # reformat into data.table with chromosome location and weight
    loadingDT = data.table(coordinateDT, loadingMat[, PCsToAnnotate, with=FALSE])
    # naming does not work if only using one PC so add this line for that case
    setnames(loadingDT, c("chr", "start", PCsToAnnotate)) 
    
    # would rounding speed up aggregation?, potentially make a sparse matrix
    # if a lot of entries became 0
    
    # specify aggregation operation
    # will be done separately for each PC specified
    aggrCommand = MIRA:::buildJ(PCsToAnnotate, 
                                rep("mean", length(PCsToAnnotate)))
    
    # do the actual aggregation
    loadAgMain = RGenomeUtils::BSAggregate(BSDT = loadingDT, regionsGRL = GRangesList(regionSet),
                jExpr = aggrCommand,
                byRegionGroup = TRUE,
                splitFactor = NULL)
    
    
    if (permute) {
        # if no cytosines from loadings overlapped with regionSet, no need for
        # permutations
        if (is.null(loadAgMain)) {
            results = as.data.table(t(rep(NA, length(PCsToAnnotate))))
            setnames(results, PCsToAnnotate)
            return(results)
        } else {
            # permutation test?
            # shuffle chromosomal coordinate labels and rerun 1000? times
            set.seed(100) # will this cause the user problems if they independently have set
            # ...the seed for other functions?
            loadAgPerm = list()
            for (i in 1:100) {
                permInd = sample(1:nrow(loadingDT), replace = FALSE)
                # reformat into data.table with chromosome location and weight
                loadingDT = data.table(coordinateDT[permInd, ], loadingMat[, PCsToAnnotate, with=FALSE])
                
                # naming does not work if only using one PC so add this line for that case
                setnames(loadingDT, c("chr", "start", PCsToAnnotate)) 
                
                loadAgPerm[[i]] = RGenomeUtils::BSAggregate(BSDT = loadingDT, 
                                                            regionsGRL = GRangesList(regionSet),
                                                            jExpr = aggrCommand,
                                                            byRegionGroup = TRUE,
                                                            splitFactor = NULL)
                
            }
            loadAgPerm = rbindlist(loadAgPerm)
            
            testCDF = lapply(X = loadAgPerm[, .SD, .SDcols = PCsToAnnotate], ecdf)
            pVals = mapply(function(x, y) y(x), loadAgMain[, .SD, .SDcols = PCsToAnnotate],
                           testCDF)
            pVals2 = (0.5 - abs(pVals - 0.5)) * 2 # two sided?    
            
            # UPDATE to return p value from permutation test?
            return(pVals2)
            #return(loadAgMain[, PCsToAnnotate])
        }
    } else {
        # if no cytosines from loadings were included in regionSet, result is NA
        if (is.null(loadAgMain)) {
            results = as.data.table(t(rep(NA, length(PCsToAnnotate))))
            setnames(results, PCsToAnnotate)
        } else {
            results = loadAgMain[, .SD, .SDcols = PCsToAnnotate]
        }
        return(results)
    }
}


#' Wrapper function to do PCA region set enrichment 
#' analysis for many region sets
#' 
#' For parallel processing, region sets are split up between the cores

pcRegionSetEnrichment <- function(loadingMat, coordinateDT, GRList, 
                  PCsToAnnotate = c("PC1", "PC2"), permute=TRUE) {
 
    # apply over the list of region sets
    resultsList = MIRA:::lapplyAlias(GRList, 
                              function(x) aggregateLoadings(loadingMat=loadingMat, 
                                                            coordinateDT=coordinateDT, 
                                                            regionSet=x, 
                                                            PCsToAnnotate = PCsToAnnotate,
                                                            permute=permute))
    resultsDT = do.call(rbind, resultsList) 
    row.names(resultsDT) = row.names(GRList)
    
    return(resultsDT)
}    



#' function to create profile that indicates enrichment 
#' of cytosines with high loading values in region set but not in
#' surrounding genome
#'
#' A wrapper of BSAggregate that first bins regions and then aggregates
#' each bin across a set of regions, individually.
#'
pcEnrichmentProfile = function(loadingMat, coordinateDT, GRList,
                    PCsToAnnotate = c("PC1", "PC2"), binNum = 25) {
    
    
    loadingDT = as.data.table(abs(loadingMat))
    loadingDT = cbind(coordinateDT, loadingDT)
    
    GRDTList = lapply(X = GRList, FUN = MIRA:::grToDt)
    profileList = lapply(X = GRDTList, function(x) BSBinAggregate(BSDT = loadingDT, 
                                                    rangeDT = x, 
                                                    binCount = binNum, 
                                                    minReads = 0, 
                                                    byRegionGroup = TRUE, 
                                                    splitFactor = NULL,
                                                    PCsToAnnotate = PCsToAnnotate))
    
    profileList = lapply(profileList, FUN = makeSymmetric)
    profileList = lapply(profileList, function(x) x[, regionGroupID := 1:binNum])
    
    return(profileList)
}

makeSymmetric = function(prof) {
    symProf = apply(prof, 2, .makeSymmetric)
    symProf = as.data.table(symProf)
    return(symProf)
}

.makeSymmetric = function(vec) {
    symVec = (vec + rev(vec)) / 2
    return(symVec)
}

#' Produced originally for binning Ewing RRBS data across various region sets
#'
#' @param rangeDT A data table with the sets of regions to be binned, 
#' with columns named start, end
#' @param binCount Number of bins across the region
#' @param byRegionGroup Pass along to binCount (see ?binCount)
#' @param minReads Filter out bins with fewer than X reads before returning.
BSBinAggregate = function(BSDT, rangeDT, binCount, minReads = 500, 
                          byRegionGroup = TRUE, 
                          splitFactor = NULL,
                          PCsToAnnotate=PCsToAnnotate) {
    if (! "data.table" %in% class(rangeDT)) {
    stop("rangeDT must be a data.table")
}
    seqnamesColName = "seqnames"  # default column name
    if (! "seqnames" %in% colnames(rangeDT)) {
        if ("chr" %in% colnames(rangeDT)) {
            message("seqnames column name set to: chr")
            seqnamesColName = "chr"
        } else {
            # Got neither.
            stop("rangeDT must have a seqnames column")
        }
    }
    
    message("Binning...")
    binnedDT = rangeDT[, MIRA:::binRegion(start, end, binCount, get(seqnamesColName))]
    binnedGR = sapply(split(binnedDT, binnedDT$binID), dtToGr)
    message("Aggregating...")
    binnedBSDT = RGenomeUtils::BSAggregate(BSDT, 
                             regionsGRL=GRangesList(binnedGR), 
                             jExpr=buildJ(PCsToAnnotate, rep("mean", length(PCsToAnnotate))), 
                                          byRegionGroup=byRegionGroup, 
                                          splitFactor = splitFactor)
    # # If we aren't aggregating by bin, then don't restrict to min reads!
    # if (byRegionGroup) {
    #     binnedBSDT = binnedBSDT[readCount > minReads,]
    # }
    return(binnedBSDT)
}


# modification of BSAggregate to just return mean per region
# averageByRegion(BSDT = BSDT, regionsGRL, jCommand = MIRA:::buildJ(cols = "methylProp", "mean"))
averageByRegion <- function(loadingMat, coordinateDT, GRList, PCsToAnnotate = c("PC1", "PC2")) {
    
    # old parameters
    excludeGR = NULL
    splitFactor = NULL
    regionsGRL.length = NULL
    keep.na = FALSE
    keepCols = NULL
    sumCols = NULL
    byRegionGroup = FALSE
    
    # reformating parameters and assigning new names
    regionsGRL = GRList
    loadingDT = as.data.table(abs(loadingMat))
    BSDT = cbind(coordinateDT, loadingDT[, .SD, .SDcols = PCsToAnnotate])
    jCommand = MIRA:::buildJ(PCsToAnnotate, rep("mean", length(PCsToAnnotate)))
    
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
    
    bsgr <- MIRA:::BSdtToGRanges(list(BSDT));
    
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
        if (length(regionsGRL) > 100) {
            message(cleanws("BSAggregate: Calculating sizes. You can speed this
                            up by supplying a regionsGRL.length vector..."),
                    appendLF = FALSE)
        }
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
    
    avPerRegion = merge(bsCombined, region2group)
    avPerRegion[, c("regionID", "withinGroupID", "regionGroupID") := NULL]
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


#' different scoring metrics
#' remniscent of LOLA: 
#' support (number of regions, for us regions that have at least 1 cytosine and can be scored)
#' mean loading value (or could do ratio of peak to surrounding regions)
#' signal to noise ratio, how big is peak compared to noise of surrounding area
#' with SNR, even a small peak could have a high SNR




