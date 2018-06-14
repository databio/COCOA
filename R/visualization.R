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
#' (ie transformed methylation data). Must have subject_ID as row names,
#' These same subject_IDs must be column names of methylData
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
    # pcaData must have subject_ID as row name
    thisRSMData = thisRSMData[names(sort(pcaData[, pc])), ]
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

#' Heatmap of average loading score for individual regions across PCs.
#' For a single region set.
#' @param

regionPCHeatmap <- function() {
    
}



# looking at methylation level data at individual cytosines ordered by PC 
# only looking at regions with high average loading scores
# still individual cytosine methylation

#' raw methylation at top enriched regions for a single region set and single PC.
#' Patients are ordered by PC score for given PC
#' #' TODO: deal with bug when nrow(highVariable) = 1
#' @param loadingMat
#' @param loadingThreshold Only select regions with average 
#' loading at least this high. Based on loading values from orderByPC.
#' @param pcScores
#' @param coordinateDT
#' @param methylData
#' @param GRList
#' @param orderByPC PC to order patients by (order rows of heatmap by PC score)
#' REMOVE: rsInd 
#' @param topXRegions max number of regions to plot, avoids excessively large 
#' plots which can be hard to load. Number of regions on plot will be less
#' than or equal to topXRegions (less than if there are not that many regions
#' total). 50 is arbitrary 
#' 



methylAlongPC <- function (loadingMat, loadingThreshold, 
                           pcScores, coordinateDT, methylData, 
                           GRList, orderByPC, 
                           topXRegions=50) {
    
    # number of region sets to plot
    nRSToPlot = length(GRList)
    rsNames = names(GRList)
  
        
    # once for each plot/region set
    for (i in seq_along(GRList)) { # loop through top region sets
        regionSet = GRList[[i]] 
        regionSetName = rsNames[i]
        
        length(regionSet)
        # if no overlap, will return NULL
        regionLoadAv = averageByRegion(loadingMat = mPCA$rotation, 
                                       coordinateDT = bigSharedC$coordinates, 
                                       GRList = regionSet, 
                                       PCsToAnnotate = orderByPC)
        
        #if (!is.null(regionLoadAv)) {
            
            # finding a suitable threshold for "high" average loading score
            # loadingMeans = apply(X = abs(mPCATop10$rotation[, PCsToAnnotate]), 2, mean)
            # getting 95th percentile for each PC
            # loadingXPerc = apply(abs(mPCA$rotation[, PCsToAnnotate]), 2, function(x) quantile(x, 0.95))
            loadingXPerc = quantile(abs(mPCA$rotation[, orderByPC]), loadingThreshold)
            
            # hist(mPCATop10$rotation[, "PC1"])
            highVariable = regionLoadAv[get(orderByPC) > loadingXPerc, .(chr, start, end, score=get(orderByPC))]
            
            # reducing to top X regions so plot won't be too large
            if (nrow(highVariable) > topXRegions) {
                highVariable[, rowIndex :=  1:nrow(highVariable)]
                tmp = highVariable[order(score, decreasing = TRUE), ]
                tmp = tmp[1:50,]
                tmp = tmp[order(rowIndex, decreasing = FALSE), ]
                highVariable = highVariable[tmp$rowIndex, ]
            }
            
            nrow(regionLoadAv)
            if (nrow(highVariable) > 1) {
                regionSet = MIRA:::dtToGr(highVariable)
                
                
                # text(paste0(pc1$rsDescription[i], ":", pc1$rsName[i]))
                # gives error if nrow(highVariable = 1) (happened for PC2)
                multiHM = grid.grabExpr(draw(rsMethylHeatmap(methylData = methylData, 
                                                             coordGR = MIRA:::dtToGr(coordinateDT), 
                                                             regionSet = regionSet, 
                                                             pcaData = pcScores, 
                                                             pc = orderByPC, column_title= regionSetName))) # use_raster=TRUE, raster_device="jpeg")
                pushViewport(viewport(y = unit((8.5*nRSToPlot)-(i-1)*8.5, "in"), height = unit(8, "in"), just = "top"))
                grid.draw(multiHM)
                popViewport()
                # 
                #                 name = paste0(pc1$rsDescription[i], " : ", pc1$rsName[i]), 
                #                 column_title = paste0(pc1$rsDescription[i], " : ", pc1$rsName[i]),
                #                 column_title_side = "top",
                #                 column_title_gp = gpar(fontsize = 14))
                #                 # column_title = paste0(pc1$rsDescription[i], " : ", pc1$rsName[i]))
            
            
            
        }
        
        #}
        
    }
}



#' plot individual region scores/percentiles across PCs for a single region set
#' One plot for each region set
#' @param loadingMat
#' @param coordinateDT
#' @param GRList
#' @param rsNames
#' @param PCsToAnnotate

regionQuantileByPC <- function(loadingMat, coordinateDT, GRList, 
                               rsNames, PCsToAnnotate=paste0("PC", 1:5)) {
    
    for (i in seq_along(GRList)) { # loop through top region sets
        
        # rsIndex = sort.int(rsEnrichment$PC2, index.return = TRUE, decreasing = TRUE)$ix[6]
        regionSet = GRList[[i]]
        rsRegionAverage = averageByRegion(loadingMat = mPCA$rotation, coordinateDT = bigSharedC$coordinates, 
                                          GRList = regionSet, PCsToAnnotate = PCsToAnnotate,
                                          returnQuantile = TRUE)
        # ranking in terms of percentiles in case there were different distributions of loading scores for each PC
        
        
        multiHM = grid.grabExpr(draw(Heatmap(matrix = as.matrix(rsRegionAverage[, PCsToAnnotate, with=FALSE]), column_title = rsNames[i], 
                                             cluster_columns = FALSE, name = "Percentile of Loading Scores in PC"))) # use_raster=TRUE, raster_device="jpeg")
        pushViewport(viewport(y = unit((8.5*length(GRList))-(i-1) * 8.5, "in"), height = unit(8, "in"), just = "top"))
        grid.draw(multiHM)
        popViewport()
        
    }
    
}
