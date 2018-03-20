# package to annotate PCA components of DNA methylation data
# based on region set enrichment

#' Function to aggregate PCA loading weights over a given region set
#' and then get p value for each PC based on a permutation
#' 
#' @param loadingMat 
#' @param regionSet A genomic ranges object with regions corresponding
#' to the same biological annotation.

aggregateLoadings(loadingMat, coordinateDT, regionSet, 
                  PCsToAnnotate = c("PC1", "PC2")) {
    # extreme positive or negative values both give important information
    loadingMat = abs(loadingMat) 
    
    # reformat into data.table with chromosome location and weight
    loadingDT = data.table(coordinateDT, loadingMat[, PCsToAnnotate])
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
                jCommand = aggrCommand,
                byRegionGroup = TRUE,
                splitFactor = NULL)
    
    # permutation test?
    # shuffle chromosomal coordinate labels and rerun 1000? times
    set.seed(100) # will this cause the user problems if they independently have set
    # ...the seed for other functions?
    loadAgPerm = list()
    for (i in 1:1000) {
        permInd = sample(1:nrow(loadingDT), replace = FALSE)
        # reformat into data.table with chromosome location and weight
        loadingDT = data.table(coordinateDT[permInd, ], loadingMat[, PCsToAnnotate])
        
        # naming does not work if only using one PC so add this line for that case
        setnames(loadingDT, c("chr", "start", PCsToAnnotate)) 
        
        loadAgPerm[[i]] = RGenomeUtils::BSAggregate(BSDT = loadingDT, 
                                                    regionsGRL = GRangesList(regionSet),
                                                    jCommand = aggrCommand,
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


#' Wrapper function to do PCA region set enrichment 
#' analysis for many region sets
#' 
#' For parallel processing, region sets are split up between the cores

pcRegionSetEnrichment(loadingMat, coordinateDT, GRList, 
                  PCsToAnnotate = c("PC1", "PC2")) {
 
    # apply over the list of region sets
    resultsList = lapplyAlias(GRList, 
                              function(x) aggregateLoadings(loadingMat=loadingMat, 
                                                            coordinateDT=coordinateDT, 
                                                            GRList=x, 
                                                            PCsToAnnotate = PCsToAnnotate))
    resultsDT = do.call(rbind, resultsList) 
    row.names(resultsDT) = row.names(GRList)
    
    return(resultsDT)
}    