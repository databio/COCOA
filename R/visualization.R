################################################################################
# functions to visualize results of PCRSA, relevant regions, 
# and variation in the dataset
###############################################################################
# Some imported functions:
# data.table: copy, :=, setorder
# MIRA: dtToGr
# ComplexHeatmap
###
# plotting functions to be exported:
# rsScoreHeatmap (heatmap of top region sets (rows) by PCs (cols), ordered by one PC), uses ComplexHeatmap
# rsMethylHeatmap (raw methylation along PC), uses ComplexHeatmap
# regionQuantileByPC?

# functions to visualize results of PCRSA, relevant regions, and variation in the dataset
# 

#' Look at methylation in regions of interest across samples, 
#' with samples ordered according to PC of interest 
#'
#' @param methylData DNA methylation levels (0 to 1) in matrix or data.frame. 
#' Rows are cytosines. Columns are samples.
#' @param mCoord a GRanges object or data frame with coordinates 
#' for the cytosines included in the PCA. Coordinates should be in the 
#' same order as the methylation data and loadings. If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.
#' @param regionSet A genomic ranges object with regions corresponding
#' to the same biological annotation.
#' @param pcScores The principal component scores for the samples 
#' (ie transformed methylation data). Must have subject_ID as row names,
#' These same subject_IDs must be column names of methylData
#' @param orderByPC PC to order samples by (order rows of heatmap by PC score, 
#' from high to low score)
#' @param ... optional parameters for ComplexHeatmap::Heatmap()
#  library(ComplexHeatmap)
#' @return A heatmap of DNA methylation levels in regions of interest (regionSet).
#'
#' @examples data("brcaMethylData1")
#' data("brcaCoord1")
#' data("esr1_chr1")
#' data("brcaPCScores")
#' PCRSA:::rsMethylHeatmap(methylData=brcaMethylData1,
#'                 mCoord=brcaCoord1,
#'                 regionSet=esr1_chr1,
#'                 pcScores=brcaPCScores,
#'                 orderByPC="PC1", cluster_columns=TRUE)
rsMethylHeatmap <- function(methylData, mCoord, regionSet, 
                            pcScores, orderByPC="PC1", ...) {
    

    
    # test for appropriateness of inputs/right format
    if (is(mCoord, "GRanges")) {
        coordGR <- mCoord
    } else if (is(mCoord, "data.frame")) {
        # UPDATE: does the work on data.frames that are not data.tables?
        coordGR <- dtToGr(mCoord)
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    # PCA object must have subject_ID as row.names (corresponding 
    # to column names of methylData)
    if (sum(row.names(pcScores) %in% colnames(methylData)) < 2) {
        stop(cleanws("Sample names on pca data (row names) 
                      must match sample names on methylation
                             (column names)"))
    }
    
    
    
    # coordGR =
    olList <- findOverlaps(regionSet, coordGR)
    # regionHitInd <- sort(unique(queryHits(olList)))
    cytosineHitInd <- sort(unique(subjectHits(olList)))
    thisRSMData <- t(methylData[cytosineHitInd, ])
    subject_ID <- row.names(thisRSMData)
    # centeredPCAMeth <- t(apply(t(methylData), 1, 
    #                            function(x) x - pcaData$center)) #center first 
    # reducedValsPCA <- centeredPCAMeth %*% pcaData$rotation
    # reducedValsPCA <- pcaData$x
    # pcaData must have subject_ID as row name
    thisRSMData <- thisRSMData[names(sort(pcScores[, orderByPC], 
                                          decreasing = TRUE)), ]
    message(paste0("Number of cytosines: ", ncol(thisRSMData)))
    message(paste0("Number of regions: ", length(unique(queryHits(olList)))))
    if (hasArg("cluster_columns")) {
        ComplexHeatmap::Heatmap(thisRSMData, cluster_rows = FALSE, ...)  
    } else {
        ComplexHeatmap::Heatmap(thisRSMData, cluster_rows = FALSE, 
                               cluster_columns = FALSE, ...)# ,
        # use_raster=TRUE, raster_device = "png")    
    }
   
}



#' Heatmap of enrichment scores across PCs
#' A visualization of the enrichment data.frame.
#' 
#' @param rsScores a data.table with scores for each 
#' region set from main PCRSA function. 
#' Each row is a region set. Columns are PCs and info on region set overlap
#' with DNA methylation data. Should be in the same order as GRList (the list of 
#' region sets used to create it.)
#' @param PCsToAnnotate A character vector with principal components to 
#' include. eg c("PC1", "PC2")
#' @param orderByPC PC to order by (decreasing order) in heatmap
#' @param rsNameCol character. Name of the column in rsScores that has the
#' names/identifiers for the region sets so this information can be included 
#' in the plot.
#' @param topX Number of top region sets to include in the heatmap
#' @return A heatmap of region set scores across. Each row is a region set,
#' each column is a PC. The color corresponds to a region set's relative
#' rank for a given PC out of all tested region sets.
#
# @examples scoreHeatmap <- rsScoreHeatmap(rsScores, 
#           PCsToAnnotate=paste0("PC", 1:10), orderByPC = "PC2")

rsScoreHeatmap <- function(rsScores, PCsToAnnotate=paste0("PC", 1:5),
                            orderByPC="PC1", rsNameCol = "rsName", topX = 20) {
    
    rsEnrichment <- rsScores
    # prevent indexing out of bounds later
    if (nrow(rsEnrichment) < topX) {
        topX = nrow(rsEnrichment)
    }
    # so by reference operations will not affect original object
    rsEn <- as.data.table(data.table::copy(rsEnrichment))
    
    # only ones you have data for
    PCsToAnnotate <- PCsToAnnotate[PCsToAnnotate %in% colnames(rsEn)]
    if (length(PCsToAnnotate) == 0) {
        stop("Please check format of PC names in PCsToAnnotate.")
    }
     
    # apparently erases row names
    rsEn <- rsEn[, c(PCsToAnnotate, rsNameCol), with=FALSE] 

    
    # how to deal with NA?
    
    # number of region sets tested
    rsNum <- nrow(rsEn)
    
    # convert to data.table to do some data.table operations
    rsEn = as.data.table(rsEn)
    
    for (i in seq_along(PCsToAnnotate)) {
        # first convert to rank
        setorderv(rsEn, PCsToAnnotate[i], order = -1L) # descending order
        rsEn[, PCsToAnnotate[i] := 1:rsNum]
        
        # center around zero
        rsEn[, PCsToAnnotate[i] := ((rsNum + 1) / 2) - get(PCsToAnnotate[i])]
    }
    
    # heatmap of the centered ranks
    setorderv(rsEn, orderByPC, order = -1L) # back to first order
    rowNames <-  rsEn[, get(rsNameCol)] # redefined/reordered later
    row.names(rsEn) <- rowNames
    rsEn[, c(rsNameCol) := NULL]
    rsEn <- as.matrix(rsEn)
    row.names(rsEn) <- rowNames
    Heatmap(rsEn[1:topX, ], cluster_rows = FALSE, cluster_columns = FALSE, 
            show_row_names = TRUE, row_names_max_width = unit(100000, "mm"))
    
}



#' create pdf with multiple heatmap plots (number = length(PCsToRankBy)). 
#' Plot i will be ranked by PCsToRankBy[i]. A wrapper for rsScoreHeatmap
#' 
#' @param rsScores a data.table with scores for each 
#' region set from main PCRSA function. 
#' Each row is a region set. Columns are PCs and info on region set overlap
#' with DNA methylation data. Should be in the same order as GRList (the list of 
#' region sets used to create it.)
#' @param PCsToRankBy PC to order by (decreasing order) in heatmap. One 
#' heatmap for each PC in PCsToRankBy. 
#' @param PCsToInclude A character vector with names of PCs that 
#' should be present in the heatmap.
#' @param fileName A character vector. All plots from this function will
#' be saved to a single pdf. fileName should give the name of that file.
#' By default, it will be saved in the working directory but filename
#' can also include the a file path to save the plot in another directory
#' @param topX Number of top region sets to include in the heatmap
#' 
#' # see https://github.com/jokergoo/ComplexHeatmap/issues/110
#' 
comparePCHeatmap <- function(rsScores, PCsToRankBy=paste0("PC", 1:5), 
                             PCsToInclude=paste0("PC", 1:10), fileName=NULL,
                             topX=40) {
    rsEnrichment <- rsScores
    
    if (!is.null(fileName)) {
        grDevices::pdf(file = fileName, 
                       width = 11, 
                       height = 8.5 * length(PCsToRankBy))
        grid.newpage()
        for (i in seq_along(PCsToRankBy)) {
            multiHM <- grid.grabExpr(draw(rsScoreHeatmap(rsScores = rsEnrichment, 
                                                         PCsToAnnotate = PCsToInclude,
                                                         orderByPC = PCsToRankBy[i], 
                                                         rsNameCol = "rsName", 
                                                         topX = topX)))
            
            pushViewport(viewport(y = unit((8.5 * length(PCsToRankBy)) -
                                               (i - 1) * 8.5, "in"), 
                                  height = unit(8, "in"), just = "top"))
            grid.draw(multiHM)
            popViewport()
        }
        dev.off()
    }
    
    # multiColPlots <- marrangeGrob(grobs = plotList, ncol = 1, nrow = 1)
    # ggsave(filename <- paste0(Sys.getenv("PLOTS"), "rsScoreHeatmap.pdf"), plot = multiColPlots, device = "pdf")
    # gridextra
    # multiColPlots <- marrangeGrob(plotList, ncol = 2, nrow = 2)
}


# looking at methylation level data at individual cytosines ordered by PC 
# only looking at regions with high average loading scores
# still individual cytosine methylation

# raw methylation at top enriched regions for a single region set and single PC.
# Patients are ordered by PC score for given PC
# #' TODO: deal with bug when nrow(highVariable) = 1
# @param loadingMat matrix of loadings (the coefficients of 
# the linear combination that defines each PC). One named column for each PC.
# One row for each original dimension/variable (should be same order 
# as original data/mCoord). The x$rotation output of prcomp().
# @param loadingThreshold Only select regions with average 
# loading at least this high. Based on loading values from orderByPC.
# @param pcScores The principal component scores for the samples
#  (ie transformed methylation data). The $x output of prcomp() but must 
# have subject_ID as row names.
# These same subject_IDs must be column names of methylData.
# @param mCoord a GRanges object or data frame with coordinates 
# for the cytosines included in the PCA. Coordinates should be in the 
# same order as the methylation data and loadings. If a data.frame, 
# must have chr and start columns. If end is included, start 
# and end should be the same. Start coordinate will be used for calculations.
# @param methylData DNA methylation levels (0 to 1) in matrix or data.frame. 
# Rows are cytosines. Columns are samples.
# @param GRList GRangesList object. Each list item is 
# a distinct region set (regions that correspond to 
# the same biological annotation).
# @param orderByPC PC to order samples by (order rows of heatmap by PC score, 
# from high to low score)
# @param topXRegions max number of regions to plot, avoids excessively large 
# plots which can be hard to load. Number of regions on plot will be less
# than or equal to topXRegions (less than if there are not that many regions
# total). 50 is arbitrary 

methylAlongPC <- function (loadingMat, loadingThreshold, 
                           pcScores, mCoord, methylData, 
                           GRList, orderByPC, 
                           topXRegions=50) {
    
    if (is(mCoord, "GRanges")) {
        coordinateDT <- grToDt(mCoord)
    } else if (is(mCoord, "data.frame")) {
        coordinateDT <- mCoord
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    
    # PCA object must have subject_ID as row.names (corresponding 
    # to column names of methylData)
    if (sum(row.names(pcScores) %in% colnames(methylData)) < 2) {
        stop(cleanws("Sample names on pca data (row names) 
                             must match sample names on methylation
                             (column names)"))
    }
    
    
    
    # number of region sets to plot
    nRSToPlot <- length(GRList)
    rsNames <- names(GRList)
  
        
    # once for each plot/region set
    for (i in seq_along(GRList)) { # loop through top region sets
        regionSet <- GRList[[i]] 
        regionSetName <- rsNames[i]
        
        length(regionSet)
        # if no overlap, will return NULL
        regionLoadAv <- averageByRegion(loadingMat = loadingMat, 
                                       mCoord = coordinateDT, 
                                       regionSet = regionSet, 
                                       PCsToAnnotate = orderByPC)
        
        #if (!is.null(regionLoadAv)) {
            
            # find threshold for these loadings
            loadingXPerc <- quantile(abs(loadingMat[, orderByPC]), 
                                     loadingThreshold)
            
            # hist(loadingMat[, "PC1"])
            highVariable <- regionLoadAv[get(orderByPC) > loadingXPerc, 
                                         .(chr, 
                                           start, 
                                           end, 
                                           score=get(orderByPC))]
            
            # reducing to top X regions so plot won't be too large
            if (nrow(highVariable) > topXRegions) {
                highVariable[, rowIndex :=  1:nrow(highVariable)]
                tmp <- highVariable[order(score, decreasing = TRUE), ]
                tmp <- tmp[1:50,]
                tmp <- tmp[order(rowIndex, decreasing = FALSE), ]
                highVariable <- highVariable[tmp$rowIndex, ]
            }
            
            nrow(regionLoadAv)
            if (nrow(highVariable) > 1) {
                regionSet <- dtToGr(highVariable)
                
                
                # text(paste0(pc1$rsDescription[i], ":", pc1$rsName[i]))
                # gives error if nrow(highVariable = 1) (happened for PC2)
                multiHM <- grid.grabExpr(draw(rsMethylHeatmap(methylData = methylData, 
                                                             mCoord = dtToGr(coordinateDT), 
                                                             regionSet = regionSet, 
                                                             pcScores = pcScores, 
                                                             orderByPC = orderByPC, 
                                                             column_title= regionSetName))) # use_raster=TRUE, raster_device="jpeg")
                pushViewport(viewport(y = unit((8.5 * nRSToPlot) - 
                                                   (i - 1) * 8.5, "in"), 
                                      height = unit(8, "in"), 
                                      just = "top"))
                grid.draw(multiHM)
                popViewport()
                # column_title = paste0(pc1$rsDescription[i], " : ", pc1$rsName[i]),
                # column_title_side = "top",
                # column_title_gp = gpar(fontsize = 14))
            
            
        }
        
        #}
        
    }
}



#' plot individual region scores/percentiles across PCs for a single region set
#' One plot for each region set
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
#' @param rsNames character vector. Names of the region sets in the same
#' order as GRList. For use as a title for each heatmap.
#' @param PCsToAnnotate A character vector with principal components to 
#' include. eg c("PC1", "PC2")
#' @param maxRegionsToPlot how many top regions from region set to include
#' in heatmap. Including too many may slow down computation and increase memory
#' use.
#' @param cluster_rows Boolean, whether to cluster rows or not (may 
#' increase computation time significantly for large number of rows)
#' 
regionQuantileByPC <- function(loadingMat, mCoord, GRList, 
                               rsNames, PCsToAnnotate=paste0("PC", 1:5),
                               maxRegionsToPlot = 8000, cluster_rows = TRUE) {
    
    if (is(mCoord, "GRanges")) {
        coordinateDT <- grToDt(mCoord)
    } else if (is(mCoord, "data.frame")) {
        coordinateDT <- mCoord
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    
    for (i in seq_along(GRList)) { # loop through top region sets
        
        regionSet <- GRList[[i]]
        
    
        rsRegionAverage <- averageByRegion(loadingMat = loadingMat, 
                                           mCoord =coordinateDT, 
                                          regionSet = regionSet, 
                                          PCsToAnnotate = PCsToAnnotate,
                                          returnQuantile = TRUE)
        # ranking in terms of percentiles in case there were different 
        # distributions of loading scores for each PC
        
        # if there are too many regions, don't plot because the attempt to cluster
        # will cause a memory error
        if (nrow(rsRegionAverage) <= maxRegionsToPlot) {
            multiHM <- grid.grabExpr(draw(Heatmap(matrix = as.matrix(rsRegionAverage[, PCsToAnnotate, with=FALSE]), 
                                                 column_title = rsNames[i], 
                                                 cluster_rows = cluster_rows,
                                                 cluster_columns = FALSE, 
                                                 name = "Percentile of Loading Scores in PC"))) # use_raster=TRUE, raster_device="jpeg")
            pushViewport(viewport(y = unit((8.5*length(GRList))-(i-1) * 8.5, "in"), 
                                  height = unit(8, "in"), just = "top"))
            grid.draw(multiHM)
            popViewport()
        }

        
    }
    
}
