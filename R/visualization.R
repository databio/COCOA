# functions to visualize results of PCRSA, relevant regions, and variation in the dataset


#' Function to look at methylation in regions of interest across patients, 
#' with patients ordered according to PC of interest 
# https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/
simpleCache(allMPCAString, assignToVariable = "mPCAALL")
simpleCache(top10MPCAString, assignToVariable = "mPCATop10")
simpleCache("rsEnrichmentTop10")
simpleCache("rsEnrichment")

rsEnrichmentTop10[, rsIndex := 1:nrow(rsEnrichmentTop10)]
rsEnrichment[, rsIndex := 1:nrow(rsEnrichment)]

# getting appropriate region
pc2 = rsEnrichment[order(PC2, decreasing = TRUE), ]
# GRList[[pc1_top10$rsIndex[1]]] # most enriched region set
rsMethylHeatmap(methylData = bigSharedC$methylProp, 
                coordGR =  MIRA:::dtToGr(bigSharedC$coordinates), 
                regionSet = GRList[[pc2$rsIndex[2]]], 
                pcaData = mPCAallC, 
                pc = "PC2")




# getting PC scores for all data (including NBM which were not used in PCA)

hist(stat5MData[, 30])
#' @param methylData
#' @param coordGR coordinates for methylation loci in methylData. (chromosome, start, end)
#' @param regionSet 
#' @param pcaData The principal component scores for the samples 
#' (ie transformed methylation data).
#' @param pc The principal component used to order the samples in the heatmap
rsMethylHeatmap <- function(methylData, coordGR, regionSet, pcaData, pc="PC1") {
    
    # coordGR =
    olList = findOverlaps(regionSet, coordGR)
    # regionHitInd = sort(unique(queryHits(olList)))
    cytosineHitInd = sort(unique(subjectHits(olList)))
    thisRSMData = t(methylData[cytosineHitInd, ])
    subject_ID = row.names(thisRSMData)
    # centeredPCAMeth = t(apply(t(methylData), 1, function(x) x - pcaData$center)) # center first 
    # reducedValsPCA = centeredPCAMeth %*% pcaData$rotation
    reducedValsPCA = pcaData$x
    thisRSMData = thisRSMData[names(sort(pcaData$x[, "PC1"])), ]
    message(paste0("Number of cytosines: ", ncol(thisRSMData)))
    message(paste0("Number of regions: ", length(unique(queryHits(olList)))))
    heatmap(thisRSMData, Rowv = NA, Colv = NA, col=cm.colors(256), scale = "none")
}
