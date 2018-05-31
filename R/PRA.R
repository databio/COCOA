# package to annotate PCA components of DNA methylation data
# based on region set enrichment
library(RGenomeUtils)
library(GenomicRanges)
library(data.table)

#' Function to aggregate PCA loading weights over a given region set
#' and then get p value for each PC based on a permutation
#' 
#' @param loadingMat 
#' @param regionSet A genomic ranges object with regions corresponding
#' to the same biological annotation.
#' #UPDATE: make sure only aggregating PCsToAnnotate to save time
#' # permute=TRUE is deprecated and that code chunk is not up to date

aggregateLoadings <- function(loadingMat, coordinateDT, regionSet, 
                  PCsToAnnotate = c("PC1", "PC2"), permute=FALSE) {
    
    numOfRegions = length(regionSet)
    
    # extreme positive or negative values both give important information
    loadingMat = abs(loadingMat) 
    loadingMat = as.data.table(loadingMat)
    
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
    # previously used BSAggregate from RGenomeUtils but now using local. 
    # modified copy
    loadAgMain = BSAggregate(BSDT = loadingDT, regionsGRL = GRangesList(regionSet),
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
                
                loadAgPerm[[i]] = BSAggregate(BSDT = loadingDT, 
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
            results[, cytosine_coverage := 0]
            results[, region_coverage := 0]
            results[, total_region_number := numOfRegions]
            results[, mean_region_size := round(mean(width(regionSet)), 1)]
        } else {
            results = loadAgMain[, .SD, .SDcols = PCsToAnnotate]
            results[, cytosine_coverage := loadAgMain[, .SD, .SDcols = "numCpGsOverlapping"]]
            results[, region_coverage := loadAgMain[, .SD, .SDcols = "numRegionsOverlapping"]]
            results[, total_region_number := numOfRegions]
            results[, mean_region_size := round(mean(width(regionSet)), 1)]
        }
        return(results)
    }
}


#' Wrapper function to do PCA region set enrichment 
#' analysis for many region sets
#'
#' For parallel processing, region sets are split up between the cores
#' @param loadingMat
#' @param 
#' @return data.table of results, one row for each region set. 
#' One column for each PC in PCsToAnnotate
#' with average loading score for that PC for a given region set. 
#' rsName column has region set name.
#' rsDescription has a description for each region set. 
#' column has number of cytosines that overlapped with the given region set.
#' column has number of regions that overlapped with any cytosines.
#' column has total number of regions. 

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


#' modification of BSAggregate to just return mean per region
#' averageByRegion(BSDT = BSDT, regionsGRL, jCommand = MIRA:::buildJ(cols = "methylProp", "mean"))
#' @param loadingMat matrix of loading values. One named column for each PC and
#' rows are the original dimensions/variables. The $rotation output of prcomp.
#' @param GRList a single region set as a GRangesList. If regionSet is a GRanges
#' object, this can be made with the command: GRangesList(regionSet)
#' @return a data.table with region coordinates and average loading 
#' values for each region. Has columns chr, start, end, and a column for each
#' PC in PCsToAnnotate. Regions are not in order along the rows of the data.table.
#' 
# Devel note: I could add a column for how many cytosines are in each region 
averageByRegion <- function(loadingMat, coordinateDT, GRList, PCsToAnnotate = c("PC1", "PC2"), returnQuantile=FALSE) {
    
    # old parameters of BSAggregate
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
    # linking coordinates to loading values, has columns chr start, PCsToAnnotate
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


#' different scoring metrics
#' remniscent of LOLA: 
#' support (number of regions, for us regions that have at least 1 cytosine and can be scored)
#' mean loading value (or could do ratio of peak to surrounding regions)
#' signal to noise ratio, how big is peak compared to noise of surrounding area
#' with SNR, even a small peak could have a high SNR

#' @param loadingProf
loadingProfileSNR = function() {
    # magnitude of the peak divided by standard deviation of 
    # noise (signal in surrounding areas)
}


# support is just number of regions from each input region set that overlap at all with
# the cytosines that have loading values 
# GenomicRanges::findOverlaps(regionSet, MIRA:::dtToGr(coordinateDT)) 

###########################################################################

# BSAggregate from RGenomeUtils
#' BSaggregate -- Aggregate a BSDT across regions or region groups,
#' for multiple samples at a time.
#' This function is as BScombineByRegion, but can handle not only multiple
#' samples in BSDT, but also simultaneously multiple region sets by passing
#' a regionsGRL (GRangesList object).
#' you can use jExpr to do other functions.

#' Given a bisulfite data table as input, with an identifier column for
#' different samples; plus a GRanges objects with regions to aggregate.
#'
#' @param BSDT The bisulfite data.table (output from one of the parsing
#' functions for methylation calls) that you wish to aggregate. It can
#' be a combined table, with individual samples identified by column passed
#' to splitFactor.
#' @param regionsGRL Regions across which you want to aggregate. Should be 
#' from a single region set. eg GRangesList(regionSet)
#' @param excludeGR A GenomicRanges object with regions you want to 
#' exclude from the aggregation function. These regions will be eliminated
#' from the input table and not counted.
#' @param jExpr You can pass a custom command in the j slot to data.table
#' specifying which columns to aggregate, and which functions to use. You
#' can use buildJ() to build a jExpr argument easily.
#' @param byRegionGroup You can aggregate by regionID or by regionGroupID; 
#' this reflects the regionsGRL that you pass; by default, BSAggregate will
#' aggregate each region individually -- scores will then be contiguous, and
#' the output is 1 row per region.
#' Turn on this flag to aggregate across all region groups, making the result
#' uncontiguous, and resulting in 1 row per *region group*.
#'
#' @export
BSAggregate = function(BSDT, regionsGRL, excludeGR=NULL, regionsGRL.length = NULL, splitFactor=NULL, keepCols=NULL, 
                       sumCols=NULL, jExpr=NULL, byRegionGroup=FALSE, keep.na=FALSE) {
    
    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if( "GRanges" %in% class(regionsGRL)) {
        regionsGRL = GRangesList(regionsGRL)
    } else if (! "GRangesList" %in% class(regionsGRL)) {
        stop("regionsGRL is not a GRanges or GRangesList object")
    }
    
    if(! is.null(excludeGR)) {
        BSDT = BSFilter(BSDT, minReads=0, excludeGR)
    }
    
    bsgr = MIRA:::BSdtToGRanges(list(BSDT))
    
    additionalColNames = setdiff(colnames(BSDT), c("chr","start", "end","hitCount","readCount", splitFactor))
    
    colModes = sapply(BSDT,mode)
    if (is.null(sumCols)) {
        sumCols = setdiff(colnames(BSDT),c("chr", "start", "end", "strand", splitFactor, keepCols))
        # Restrict to numeric columns.		
        sumCols = intersect(sumCols, names(colModes[which(colModes == "numeric")]))
        
    }
    # It's required to do a findoverlaps on each region individually,
    # Not on a GRL, because of the way overlaps with GRLs work. So,
    # we must convert the GRL to a GR, but we must keep track of which
    # regions came from which group.
    regionsGR = unlist(regionsGRL)
    
    if(is.null(regionsGRL.length)) {
        if (length(regionsGRL) > 100) {
            message("BSAggregate: Calculating sizes. You can speed this up by supplying a regionsGRL.length vector...", appendLF=FALSE)
        }
        regionsGRL.length = sapply(regionsGRL, length)
        message("Done counting regionsGRL lengths.")
    }
    
    # Build a table to keep track of which regions belong to which group
    region2group = data.table(
        regionID=1:length(regionsGR), 
        chr=as.vector(seqnames(regionsGR)), 
        start=as.vector(start(regionsGR)), 
        end=as.vector(end(regionsGR)),
        withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
        regionGroupID=rep(1:length(regionsGRL), regionsGRL.length))
    setkey(region2group, regionID)
    
    
    message("Finding overlaps...")
    fo = findOverlaps(query = bsgr[[1]], subject = regionsGR)
    
    ### use info from findOverlaps to see how many individual
    # cytosines and how many regions overlap with each other
    numCpGsOverlapping = length(unique(queryHits(fo)))
    numRegionsOverlapping = length(unique(subjectHits(fo)))
    
    
    setkey(BSDT, chr, start)
    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))
    
    message("Setting regionIDs...")
    BSDT = BSDT[queryHits(fo),] #restrict the table to CpGs in any region.
    
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
        funcs = c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
        jExpr = buildJ(cols, funcs)
    }
    message("jExpr: ", jExpr)
    
    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        agCol = "regionGroupID"
    } else {
        agCol = "regionID" # Default
    }
    
    # Build the by string
    if (is.null(splitFactor)) {
        byString = paste0("list(regionID)")
    } else {
        byString = paste0("list(", paste("regionID", paste0(splitFactor, ""), collapse=", ", sep=", "), ")")
    }
    
    # Now actually do the aggregate:
    message("Combining...")
    bsCombined = BSDT[,eval(parse(text=jExpr)), by=eval(parse(text=byString))]
    setkey(bsCombined, regionID)
    # Now aggregate across groups.
    # I do this in 2 steps to avoid assigning regions to groups,
    # which takes awhile. I think this preserve memory and is faster.
    
    # Define aggregation column. aggregate by region or by region group?
    if (byRegionGroup) {
        # must set allow=TRUE here in case there are multiple IDs (splitCol)
        bsCombined[region2group, regionGroupID:=regionGroupID, allow=TRUE]
        if (! is.null(splitFactor) ) { 
            byStringGroup = paste0("list(", paste("regionGroupID", paste0(splitFactor, collapse=", "), sep=", "), ")")
        } else {
            byStringGroup = "list(regionGroupID)"
        }
        bsCombined=bsCombined[,eval(parse(text=jExpr)), by=eval(parse(text=byStringGroup))]
        bsCombined[, numCpGsOverlapping := numCpGsOverlapping]
        bsCombined[, numRegionsOverlapping := numRegionsOverlapping]
        return(bsCombined)
    } else {
        e = region2group[bsCombined,]
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
#' @param pca The pca data. Output from prcomp().
#' @param methylData DNA methylation levels in matrix or data.frame. 
#' Rows are cytosines. Columns are 
#' samples.
#' @param coordinateDT One row per cytosine coordinate, corresponding to
#' a row of methylMat. Columns are chr, start, end.   
#' @param returnCor Option to return correlation between these scores and
#' original PC scores. 
#' 
#' @return a score for each patient from only loading values for cytosines in 
#' regionSet. If returnCor = TRUE, returns a correlation score for scores from
#' each PCofInterest with scores from a subset of loading values for that PC.

pcFromSubset <- function(regionSet, pca, methylData, coordinateDT, PCofInterest="PC1", returnCor=FALSE) {
    
    # test for appropriateness of inputs/right format
    
    # get subset of loading values
   #  subsetInd = 
    newColNames <- paste0(PCofInterest, "_subset")
    
    coordGR = MIRA:::dtToGr(coordinateDT)
    olList = findOverlaps(query = regionSet, subject = coordGR)
    # regionHitInd = sort(unique(queryHits(olList)))
    cytosineHitInd = sort(unique(subjectHits(olList)))
    thisRSMData = t(methylData[cytosineHitInd, ])
    # subject_ID = row.names(thisRSMData)
    # centeredPCAMeth = t(apply(t(methylData), 1, function(x) x - pcaData$center)) # center first 
    # reducedValsPCA = centeredPCAMeth %*% pcaData$rotation
    # reducedValsPCA = pcaData$x
    # use subset of loading values on subset of DNA methylation values
    # getting PC scores manually so artificial PCs will be included (PC1m4 and PC1p3)
    centeringSubset = pca$center[cytosineHitInd]
    centeredMeth = t(apply(thisRSMData, 1, function(x) x - centeringSubset)) # center first 
    reducedValsPCA = centeredMeth %*% pca$rotation[cytosineHitInd, PCofInterest]
    colnames(reducedValsPCA) <- newColNames
    # pcaValDF = as.data.frame(reducedValsPCA)
    
    
    # optional calculate correlation
    if (returnCor) {
        # plot(reducedValsPCA[, PCofInterest], pca$x[, PCofInterest])
        corMat = cor(reducedValsPCA[, newColNames], pca$x[, PCofInterest])
        # only correlation between subset and matching PC, not with other PCs
        mainCor = diag(corMat)
        names(mainCor) <- PCofInterest
        return(mainCor)
        # origPCCor = cor(pca$x[, PCofInterest])
        # subsetPCCor = cor(reducedValsPCA[, PCofInterest])
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
#' @param coordGR GenomicRanges object with coordinates for cytosines 
#' included in the PCA.
#' @param GRList
#' @return Total number of cytosines from methylation data that overlapped with
#' each region set. Also the percent of cytosines that overlapped with
#' other region set.

percentCOverlap <- function(coordGR, GRList) {
    
    if (is.data.table(coordGR)) {
        coordGR <- MIRA:::dtToGr(coordGR)
    }
    
    # overlap between each region set and methylation data to see 
    # which cytosines were covered by each region set
    olList = lapply(X = GRList, 
                    function(x) findOverlaps(query = x, subject = coordGR))
    
    # indices for cytosines that overlapped with each region set
    cHitInd = lapply(X = olList, function(x) sort(unique(subjectHits(x))))

    # total number
    total_cytosines = sapply(X = cHitInd, length)
    
    # make matrix of combinatorial intersection, will eventually have percent overlap
    pOL = matrix(nrow=length(GRList), ncol=length(GRList))
    
    # intersect each iten in cHitInd with each item in cHitInd
    sharedC = lapply(X = cHitInd, 
                     FUN = function(x) sapply(X = cHitInd, FUN = function(y) length(intersect(x, y))))
    
    # converting to proportion of cytosines covered by each region set
    prop_overlap = lapply(X = sharedC, function(x) x / total_cytosines)
    prop_overlap = do.call(rbind, prop_overlap)

    return(list(prop_overlap, total_cytosines))    
    
}