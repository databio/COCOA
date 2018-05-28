# functions to visualize results of PCRSA, relevant regions, and variation in the dataset
# data.table: copy, :=, setorder
# MIRA: dtToGr
# ComplexHeatmap

#' # functions to visualize results of PCRSA, relevant regions, and variation in the dataset
#' 
#' 
#' #' Function to look at methylation in regions of interest across patients, 
#' #' with patients ordered according to PC of interest 
#' # https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/
#' simpleCache(allMPCAString, assignToVariable = "mPCAALL")
#' simpleCache(top10MPCAString, assignToVariable = "mPCATop10")
#' simpleCache("rsEnrichmentTop10")
#' simpleCache("rsEnrichment")
#' 
#' rsEnrichmentTop10[, rsIndex := 1:nrow(rsEnrichmentTop10)]
#' rsEnrichment[, rsIndex := 1:nrow(rsEnrichment)]
#' 
#' # getting appropriate region
#' pc2 = rsEnrichment[order(PC2, decreasing = TRUE), ]
#' # GRList[[pc1_top10$rsIndex[1]]] # most enriched region set
#' rsMethylHeatmap(methylData = bigSharedC$methylProp, 
#'                 coordGR =  MIRA:::dtToGr(bigSharedC$coordinates), 
#'                 regionSet = GRList[[pc2$rsIndex[2]]], 
#'                 pcaData = mPCAallC, 
#'                 pc = "PC2")
#' 
#' 
#' 
#' 
#' # getting PC scores for all data (including NBM which were not used in PCA)
#' 
#' hist(stat5MData[, 30])

#' @param methylData
#' @param coordGR coordinates for methylation loci in methylData. (chromosome, start, end)
#' @param regionSet 
#' @param pcaData The principal component scores for the samples 
#' (ie transformed methylation data).
#' @param pc The principal component used to order the samples in the heatmap
# library(ComplexHeatmap)
rsMethylHeatmap <- function(methylData, coordGR, regionSet, pcaData, pc="PC1", ...) {
    
    # coordGR =
    olList = findOverlaps(regionSet, coordGR)
    # regionHitInd = sort(unique(queryHits(olList)))
    cytosineHitInd = sort(unique(subjectHits(olList)))
    thisRSMData = t(methylData[cytosineHitInd, ])
    subject_ID = row.names(thisRSMData)
    # centeredPCAMeth = t(apply(t(methylData), 1, function(x) x - pcaData$center)) # center first 
    # reducedValsPCA = centeredPCAMeth %*% pcaData$rotation
    # reducedValsPCA = pcaData$x
    thisRSMData = thisRSMData[names(sort(pcaData[, "PC1"])), ]
    message(paste0("Number of cytosines: ", ncol(thisRSMData)))
    message(paste0("Number of regions: ", length(unique(queryHits(olList)))))
    ComplexHeatmap::Heatmap(thisRSMData, cluster_rows = FALSE, cluster_columns = FALSE, ...)# ,
                            # use_raster=TRUE, raster_device = "png")
}



#' Heatmap of enrichment scores across PCs
#' A visualization of the enrichment data.table.
#' 
#' @param rsEnrichment
#' @param PCsToAnnotate A character vector with principal components to 
#' include in the plot, eg c("PC1", "PC2")
#' @param orderByPC PC to order by (decreasing order) in heatmap
#
#' @example enrichmentHeatmap = rsEnrichHeatmap(rsEnrichment, PCsToAnnotate=paste0("PC", 1:10), orderByPC = "PC2")
rsEnrichHeatmap <- function(rsEnrichment, PCsToAnnotate=paste0("PC", 1:5),
                            orderByPC="PC1", rsNameCol = "rsNames", topX = 20) {
    
    # so by reference operations will not affect original object
    rsEn = data.table::copy(rsEnrichment)
    
    # only ones you have data for
    PCsToAnnotate = PCsToAnnotate[PCsToAnnotate %in% colnames(rsEn)]
    if (length(PCsToAnnotate) == 0) {
        stop("Please check format of PC names in PCsToAnnotate.")
    }
     
    rsEn = rsEn[, c(PCsToAnnotate, rsNameCol), with = FALSE] # apparently erases row names
    
    
    # how to deal with NA?
    
    # number of region sets tested
    rsNum = nrow(rsEn)
    
    for (i in seq_along(PCsToAnnotate)) {
        # first convert to rank
        setorderv(rsEn, PCsToAnnotate[i], order = -1L) # descending order
        rsEn[, PCsToAnnotate[i] := 1:rsNum]
        
        # center around zero
        rsEn[, PCsToAnnotate[i] := ((rsNum + 1) / 2) - get(PCsToAnnotate[i])]
    }
    
    # heatmap of the centered ranks
    setorderv(rsEn, orderByPC, order = -1L) # back to first order
    rowNames =  rsEn[, get(rsNameCol)] # redefined/reordered later
    row.names(rsEn) <- rowNames
    rsEn[, c(rsNameCol) := NULL]
    rsEn = as.matrix(rsEn)
    row.names(rsEn) <- rowNames
    Heatmap(rsEn[1:topX, ], cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, row_names_max_width = unit(100000, "mm"))
    
}



#' create pdf with multiple heatmap plots (number = length(PCsToRankBy)). 
#' Plot i will be ranked by PCsToRankBy[i]. 
#' # see https://github.com/jokergoo/ComplexHeatmap/issues/110
comparePCHeatmap <- function(rsEnrichment, PCsToRankBy=paste0("PC", 1:5), PCsToInclude=paste0("PC", 1:10), fileName=NULL) {
    if (!is.null(fileName)) {
        grDevices::pdf(file = fileName, width = 11, height = 8.5 * length(PCsToRankBy))
        grid.newpage()
        for (i in seq_along(PCsToRankBy)) {
            multiHM = grid.grabExpr(draw(rsEnrichHeatmap(rsEnrichment = rsEnrichment, PCsToAnnotate = PCsToInclude,
                                                         orderByPC = PCsToRankBy[i], rsNameCol = "rsName", topX = 40)))
            
            pushViewport(viewport(y = unit((8.5*length(PCsToRankBy))-(i-1)*8.5, "in"), height = unit(8, "in"), just = "top"))
            grid.draw(multiHM)
            popViewport()
        }
        dev.off()
    }
    
    # multiColPlots = marrangeGrob(grobs = plotList, ncol = 1, nrow = 1)
    # ggsave(filename = paste0(Sys.getenv("PLOTS"), "rsEnrichHeatmap.pdf"), plot = multiColPlots, device = "pdf")
    # gridextra
    # multiColPlots = marrangeGrob(plotList, ncol = 2, nrow = 2)
}

