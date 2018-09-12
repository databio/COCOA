################################################################################
# functions to visualize results of PCRSA, relevant regions, 
# and variation in the dataset
###############################################################################
# Some imported functions:
# data.table: copy, :=, setorder
# MIRA: dtToGr
# ComplexHeatmap
###

# functions to visualize results of PCRSA, relevant regions, and variation in the dataset
# 

# color schemes: red/blue, yellow/red, red/grey, skyblue/coral, skyblue/yellow

#' Look at features (eg, DNA methylation values) in regions of 
#' interest across samples, 
#' with samples ordered according to PC of interest. 
#' The ComplexHeatmap package
#' is used and additional parameters for the ComplexHeatmap::Heatmap function
#' may be passed to this function to modify the heatmap.   
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
#' @param pcScores A matrix. The principal component scores for the samples 
#' (ie transformed methylation data). Must have subject_ID as row names,
#' These same subject_IDs must be column names of methylData
#' @param orderByPC PC to order samples by (order rows of heatmap by PC score, 
#' from high to low score)
#' @param cluster_columns boolean, whether to cluster columns (the features,
#' eg DNA methylation values for each CpG).
#' @param cluster_rows boolean, whether rows should be clustered. 
#' This should be kept as FALSE to keep the correct ranking of 
#' samples/observations according to their PC score.
#' @param name character object, legend title
#' @param ... optional parameters for ComplexHeatmap::Heatmap() (eg change
#' heatmap colors with "col" parameter)
#' @return A heatmap of feature values (eg DNA methylation levels) 
#' in regions of interest (regionSet).
#' Each row is a patient/sample and each column is an individual feature. 
#' Rows are ordered by PC score (orderByPC), high scores at top and low at 
#' the bottom.
#'
#' @examples data("brcaMethylData1")
#' data("brcaCoord1")
#' data("esr1_chr1")
#' data("brcaPCScores")
#' featuresAlongPC(methylData=brcaMethylData1,
#'                 mCoord=brcaCoord1,
#'                 regionSet=esr1_chr1,
#'                 pcScores=brcaPCScores,
#'                 orderByPC="PC1", cluster_columns=TRUE)
#' 
#' @export
# previously called rsMethylHeatmap
featuresAlongPC <- function(methylData, mCoord, regionSet, 
                            pcScores, orderByPC="PC1", cluster_columns = FALSE, 
                            cluster_rows = FALSE, name = "Feature Value", ...) {
    

    if (!(is(methylData, "matrix") || is(methylData, "data.frame"))) {
        stop("methylData should be a matrix or data.frame. Check object class.")
    }
    
    # test for appropriateness of inputs/right format
    if (is(mCoord, "GRanges")) {
        coordGR <- mCoord
    } else if (is(mCoord, "data.frame")) {
        # UPDATE: does the work on data.frames that are not data.tables?
        coordGR <- dtToGr(mCoord)
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    if (!(is(pcScores, "matrix") || is(pcScores, "data.frame"))) {
        stop("pcScores should be a matrix or data.frame.")
    }
    
    
    # PCA object must have subject_ID as row.names (corresponding 
    # to column names of methylData)
    if (sum(row.names(pcScores) %in% colnames(methylData)) < 2) {
        stop(cleanws("Sample names on pca data (row names) 
                      must match sample names on methylation
                             (column names)"))
    }
    
    
    if (!is(regionSet, "GRanges")) {
        stop("regionSet should be a GRanges object. Check object class.")
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
    ComplexHeatmap::Heatmap(thisRSMData, 
                            cluster_rows = cluster_rows, 
                            cluster_columns = cluster_columns, 
                            name = name, ...)
}


#' Heatmap of the ranking of region set scores across PCs
#' A visualization of rank of region sets in each PC, allowing the
#' user to see if a region set is ranked highly in all PCs or only a subset.
#' The ComplexHeatmap package
#' is used and additional parameters for the ComplexHeatmap::Heatmap function
#' may be passed to this function to modify the heatmap.  
#' 
#' @param rsScores a data.frame with scores for each 
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
#' @param cluster_rows boolean, whether rows should be clustered. 
#' This should be kept as FALSE to keep the correct ranking of region sets.
#' @param cluster_columns boolean, whether to cluster columns. It is recommended
#' to keep this as FALSE so it will be easier to compare PCs 
#' (with cluster_columns = FALSE, they will be in the same specified
#' order in different heatmaps)
#' @param show_row_names boolean, display row names (ie region set names)
#' @param row_names_max_width "unit" object. The amount of room to 
#' allocate for row names. See ?grid::unit for object type.
#' @param name character object, legend title
# @param col a vector of colors or a color mapping function which
# will be passed to the ComplexHeatmap::Heatmap() function. See ?Heatmap
# (the "col" parameter) for more details.
#' @param ... optional parameters for ComplexHeatmap::Heatmap(), 
#' including "col" to change heatmap color scheme.
#' @return A heatmap of region set scores across. Each row is a region set,
#' each column is a PC. The color corresponds to a region set's relative
#' rank for a given PC out of all tested region sets.
#'
#' @examples data("rsScores")
#' scoreHeatmap <- rsScoreHeatmap(rsScores, 
#'           PCsToAnnotate=paste0("PC", 1:2), orderByPC = "PC2")
#' @export

rsScoreHeatmap <- function(rsScores, PCsToAnnotate=paste0("PC", 1:5),
                            orderByPC="PC1", rsNameCol = "rsName", topX = 20, 
                           cluster_rows = FALSE, cluster_columns = FALSE, 
                           show_row_names = TRUE, 
                           row_names_max_width = unit(100000, "mm"), 
                           name="Rank within PC", ...) {
    
    
    if (!is(PCsToAnnotate, "character")) {
        stop("PCsToAnnotate should be a character object (eg 'PC1').")
    }
    
    # function is not meant for making a heatmap of a single PC but should work
    if (is(rsScores, "numeric")) {
        warning("rsScores should be a data.frame")
        rsScores <- data.frame(rsScores)
        colnames(rsScores) <- orderByPC
    }
    if (!is(rsScores, "data.frame")) {
        stop("rsScores should be a data.frame.")
    }
    
    if (!(orderByPC %in% colnames(rsScores))) {
        stop(cleanws(paste0("orderByPC parameter:", orderByPC, 
                            ", was not a column of rsScores")))
    }
    if (!(rsNameCol %in% colnames(rsScores))) {
        stop(cleanws(paste0(rsNameCol, "was not a column of rsScores. 
                            Specify the column containing region set names 
                            with the rsNameCol parameter.")))
    }
    
    # so by reference operations will not affect original object
    rsScores <- as.data.table(data.table::copy(rsScores))
    
    
    # prevent indexing out of bounds later
    if (nrow(rsScores) < topX) {
        topX = nrow(rsScores)
    }
    
    
    # only ones you have data for
    PCsToAnnotate <- PCsToAnnotate[PCsToAnnotate %in% colnames(rsScores)]
    if (length(PCsToAnnotate) == 0) {
        stop("Please check format of PC names in PCsToAnnotate.")
    }
     
    # apparently erases row names
    rsScores <- rsScores[, c(PCsToAnnotate, rsNameCol), with=FALSE] 

    
    # how to deal with NA?
    
    # number of region sets tested
    rsNum <- nrow(rsScores)
    
    # convert to data.table to do some data.table operations
    rsScores = as.data.table(rsScores)
    
    for (i in seq_along(PCsToAnnotate)) {
        # first convert to rank
        setorderv(rsScores, PCsToAnnotate[i], order = -1L) # descending order
        rsScores[, PCsToAnnotate[i] := seq_len(rsNum)]
        
        # center around zero
        rsScores[, PCsToAnnotate[i] := ((rsNum + 1) / 2) - get(PCsToAnnotate[i])]
    }
    
    # heatmap of the centered ranks
    setorderv(rsScores, orderByPC, order = -1L) # back to first order
    rowNames <-  rsScores[, get(rsNameCol)] # redefined/reordered later
    row.names(rsScores) <- rowNames
    rsScores[, c(rsNameCol) := NULL]
    rsScores <- as.matrix(rsScores)
    row.names(rsScores) <- rowNames
    Heatmap(rsScores[seq_len(topX), ], cluster_rows = cluster_rows, 
            cluster_columns = cluster_columns, 
            show_row_names = show_row_names, 
            row_names_max_width = row_names_max_width, 
            name = name, ...)
    
}


#' Plot individual region scores/percentiles across PCs for a single region set
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
#' @param regionSet A genomic ranges object with regions corresponding
#' to the same biological annotation.
#' @param rsName character vector. Names of the region sets in the same
#' order as GRList. For use as a title for each heatmap.
#' @param PCsToAnnotate A character vector with principal components to 
#' include. eg c("PC1", "PC2")
#' @param maxRegionsToPlot how many top regions from region set to include
#' in heatmap. Including too many may slow down computation and increase memory
#' use. If regionSet has more regions than maxRegionsToPlot, a number of regions 
#' equal to maxRegionsToPlot will be randomly sampled from the region set and
#' these regions will be plotted.
#' @param cluster_rows Boolean, whether to cluster rows or not (may 
#' increase computation time significantly for large number of rows)
#' @param cluster_columns boolean, whether to cluster columns. It is recommended
#' to keep this as FALSE so it will be easier to compare PCs 
#' (with cluster_columns = FALSE, they will be in the same specified
#' order in different heatmaps)
#' @param column_title character object, column title
#' @param name character object, legend title
#' @param col a vector of colors or a color mapping function which
#' will be passed to the ComplexHeatmap::Heatmap() function. See ?Heatmap
#' (the "col" parameter) for more details.
#' @param ... optional parameters for ComplexHeatmap::Heatmap()
#' @return a heatmap. This heatmap allows you to see if some regions are 
#' associated with certain PCs but not others. Also, you can see if a subset of 
#' regions in the region set are associated with PCs while another subset
#' are not associated with any PCs 
#' Columns are PCs, rows are regions. To color
#' each region, first the absolute loading values within that region are
#' averaged. Then this average is compared to the distribution of absolute
#' loading values for all individual features to get a quantile/percentile 
#' for that region. Colors are based on this quantile/percentile. 
#' The output is a Heatmap object (ComplexHeatmap package).
#' 
#' @examples data("brcaLoadings1")
#' data("brcaCoord1")
#' data("esr1_chr1")
#' data("brcaPCScores")
#' regionByPCHM <- regionQuantileByPC(loadingMat = brcaLoadings1, 
#'                                    mCoord = brcaCoord1, 
#'                                    regionSet = esr1_chr1, 
#'                                    rsName = "Estrogen Receptor Chr1", 
#'                                    PCsToAnnotate=paste0("PC", 1:2),
#'                                    maxRegionsToPlot = 8000, 
#'                                    cluster_rows = TRUE, 
#'                                    cluster_columns = FALSE, 
#'                                    column_title = rsName, 
#'                                    name = "Percentile of Loading Scores in PC")
#'                                    
#' 
#' @export
regionQuantileByPC <- function(loadingMat, mCoord, regionSet, 
                               rsName = "", PCsToAnnotate=paste0("PC", 1:5),
                               maxRegionsToPlot = 8000, cluster_rows = TRUE, 
                               cluster_columns = FALSE, column_title = rsName, 
                               name = "Percentile of Loading Scores in PC", 
                               col = c("skyblue", "yellow"), ...) {
    
    if (is(loadingMat, "matrix")) {
        loadingMat = as.data.frame(loadingMat)
    } else if (!is(loadingMat, "data.frame")) {
        stop("loadingMat should be a matrix or data.frame. Check object class.")
    }
    
    if (is(mCoord, "GRanges")) {
        coordinateDT <- grToDt(mCoord)
    } else if (is(mCoord, "data.frame")) {
        coordinateDT <- mCoord
    } else {
        stop("mCoord should be a data.frame or GRanges object.")
    }
    
    if (!is(regionSet, "GRanges")) {
        stop("regionSet should be a GRanges object. Check object class.")
    }
    
    
    # if too many regions to plot, randomly subsample regions
    if (length(regionSet) > maxRegionsToPlot) {
        regionInd <- sample(x = seq_along(regionSet), size = maxRegionsToPlot, 
                            replace = FALSE)
        # get subset
        regionSet <- regionSet[regionInd]
    }
    
    
    
    rsRegionAverage <- averageByRegion(loadingMat = loadingMat, 
                                       mCoord =coordinateDT, 
                                       regionSet = regionSet, 
                                       PCsToAnnotate = PCsToAnnotate,
                                       returnQuantile = TRUE)
    # ranking in terms of percentiles in case there were different 
    # distributions of loading scores for each PC
    
    # the heatmap
    Heatmap(matrix = as.matrix(rsRegionAverage[, PCsToAnnotate, with=FALSE]), 
            column_title = rsName, 
            cluster_rows = cluster_rows,
            cluster_columns = cluster_columns, 
            name = name, 
            col = col, ...)
}

