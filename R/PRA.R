# package to annotate PCA components of DNA methylation data
# based on region set enrichment

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
    } else {
        return(loadAgMain[, .SD, .SDcols = PCsToAnnotate])
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


