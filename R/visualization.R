################################################################################
# functions to visualize results of COCOA, relevant regions, 
# and variation in the dataset
###############################################################################
# Some imported functions:
# data.table: copy, :=, setorder
# MIRA: dtToGr
# ComplexHeatmap
###

# functions to visualize results of COCOA, relevant regions, and variation in the dataset
# 

# color schemes: red/blue, yellow/red, red/grey, skyblue/coral, skyblue/yellow

#' Visualize how genomic signal in a region set changes 
#' along principal component axis
#' 
#' Look at genomic signal (eg, DNA methylation values) in regions of 
#' interest across samples, 
#' with samples ordered according to score for PC of interest. 
#' The ComplexHeatmap package
#' is used and additional parameters for the ComplexHeatmap::Heatmap function
#' may be passed to this function to modify the heatmap.   
#'
#' @param genomicSignal The genomic signal (eg DNA methylation levels) 
#' in matrix or data.frame. 
#' Rows are individual signal/feature values. Columns are samples.
#' Must have sample names/IDs as column names,
#' These same sample names must be row names of pcScores.
#' @param signalCoord a GRanges object or data frame with coordinates 
#' for the genomic signal/original data (eg DNA methylation) 
#' included in the PCA. Coordinates should be in the 
#' same order as the original data and the loadings 
#' (each item/row in signalCoord
#' corresponds to a row in loadingMat). If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.
#' @param regionSet A genomic ranges object with regions corresponding
#' to the same biological annotation. The regions where you will visualize
#' the genomic signal. Must be from the same reference genome
#' as the coordinates for the actual data (signalCoord).
#' @param pcScores A matrix. The principal component scores for the samples 
#' (ie transformed methylation data). Must have sample names/IDs as row names,
#' These same sample names must be column names of genomicSignal
#' @param orderByPC a character object. PC to order samples by 
#' (order rows of heatmap by PC score, from high to low score). 
#' Must be the name of a column in pcScores.
#' @param cluster_rows "logical" object, whether rows should be clustered. 
#' This should be kept as FALSE to keep the correct ranking of 
#' samples/observations according to their PC score.
#' @param cluster_columns "logical" object, whether to cluster columns 
#' (the genomic signal, eg DNA methylation values for each CpG).
#' @param row_title character object, row title
#' @param column_title character object, column title
#' @param column_title_side character object, where to put the column title:
#' "top" or "bottom"
#' @param name character object, legend title
#' @param col a vector of colors or a color mapping function which
#' will be passed to the ComplexHeatmap::Heatmap() function. See ?Heatmap
#' (the "col" parameter) for more details. "#EEEEEE" is the code for a
#' color similar to white.
#' @param ... optional parameters for ComplexHeatmap::Heatmap()
#' @return A heatmap of genomic signal values (eg DNA methylation levels) 
#' in regions of interest (regionSet), with rows ordered by PC score.
#' Each row is a patient/sample and each column is an individual genomic signal value. 
#' Rows are ordered by PC score for `orderByPC`, high scores at top and low at 
#' the bottom.
#'
#' @examples data("brcaMethylData1")
#' data("brcaMCoord1")
#' data("esr1_chr1")
#' data("brcaPCScores")
#' signalHM <- signalAlongPC(genomicSignal=brcaMethylData1,
#'                              signalCoord=brcaMCoord1,
#'                              regionSet=esr1_chr1,
#'                              pcScores=brcaPCScores,
#'                              orderByPC="PC1", cluster_columns=TRUE)
#' 
#' @export
# previously called rsMethylHeatmap
signalAlongPC <- function(genomicSignal, signalCoord, regionSet, 
                            pcScores, orderByPC="PC1", cluster_columns = FALSE, 
                            cluster_rows = FALSE, 
                            row_title = "Sample",
                            column_title = "Genomic Signal", 
                            column_title_side = "bottom",
                            name = "Genomic Signal Value",
                            col = c("blue", "#EEEEEE", "red"), ...) {
    

    if (!(is(genomicSignal, "matrix") || is(genomicSignal, "data.frame"))) {
        stop("genomicSignal should be a matrix or data.frame. Check object class.")
    }
    
    # test for appropriateness of inputs/right format
    if (is(signalCoord, "GRanges")) {
        coordGR <- signalCoord
    } else if (is(signalCoord, "data.frame")) {
        # UPDATE: does the work on data.frames that are not data.tables?
        coordGR <- dtToGr(signalCoord)
    } else {
        stop("signalCoord should be a data.frame or GRanges object.")
    }
    
    if (!(is(pcScores, "matrix") || is(pcScores, "data.frame"))) {
        stop("pcScores should be a matrix or data.frame.")
    }
    
    
    # PCA object must have subject_ID as row.names (corresponding 
    # to column names of genomicSignal)
    if (sum(row.names(pcScores) %in% colnames(genomicSignal)) < 2) {
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
    thisRSMData <- t(genomicSignal[cytosineHitInd, ])
    subject_ID <- row.names(thisRSMData)
    # centeredPCAMeth <- t(apply(t(genomicSignal), 1, 
    #                            function(x) x - pcaData$center)) #center first 
    # reducedValsPCA <- centeredPCAMeth %*% pcaData$rotation
    # reducedValsPCA <- pcaData$x
    # pcaData must have subject_ID as row name
    thisRSMData <- thisRSMData[names(sort(pcScores[, orderByPC], 
                                          decreasing = TRUE)), ]
    message(paste0("Number of cytosines: ", ncol(thisRSMData)))
    message(paste0("Number of regions: ", length(unique(queryHits(olList)))))
    ComplexHeatmap::Heatmap(thisRSMData, 
                            col = col,
                            row_title = row_title,
                            column_title = column_title,
                            column_title_side = column_title_side,
                            cluster_rows = cluster_rows, 
                            cluster_columns = cluster_columns, 
                            name = name, ...)
}


#' Heatmap of region set scores
#' 
#' Heatmap of the ranking of region set scores across PCs
#' A visualization of the rank of region sets in each PC, allowing the
#' user to see if a region set is ranked highly in all PCs or only a subset.
#' Region sets will be ranked from highest scoring to lowest based on 
#' their score for `orderByPC`.
#' The ComplexHeatmap package
#' is used and additional parameters for the ComplexHeatmap::Heatmap function
#' may be passed to this function to modify the heatmap.  
#' 
#' @param rsScores a data.frame with scores for each 
#' region set from main COCOA function `runCOCOA`. 
#' Each row is a region set. Columns are scores, one column for each PCs 
#' Also can have columns with info on region set overlap
#' with the original data. Should be in the same order as GRList (the list of 
#' region sets used to create it.)
#' @param PCsToAnnotate A character vector with principal components to 
#' include. eg c("PC1", "PC2"). Must be column names of rsScores.
#' @param orderByPC a character object. PC to order by in heatmap 
#' (arranged in decreasing order for scores so p values should 
#' be -log transformed). Must be the name of a column in rsScores.
#' @param rsNameCol character. Name of the column in rsScores that has the
#' names/identifiers for the region sets so these can be included 
#' in the plot as row names.
#' @param topX Number of top region sets to include in the heatmap
#' @param row_title character object, row title
#' @param column_title character object, column title
#' @param column_title_side character object, where to put the column title:
#' "top" or "bottom"
#' @param cluster_rows "logical" object, whether rows should be clustered. 
#' This should be kept as FALSE to keep the correct ranking of region sets.
#' @param cluster_columns "logical" object, whether to cluster columns. It is recommended
#' to keep this as FALSE so it will be easier to compare PCs 
#' (with cluster_columns = FALSE, they will be in the same specified
#' order in different heatmaps)
#' @param show_row_names "logical" object, display row names (ie region set names)
#' @param row_names_max_width "unit" object. The amount of room to 
#' allocate for row names. See ?grid::unit for object type.
#' @param name character object, legend title
#' @param col a vector of colors or a color mapping function which
#' will be passed to the ComplexHeatmap::Heatmap() function. See ?Heatmap
#' (the "col" parameter) for more details. "#EEEEEE" is the code for a
#' color similar to white.
#' @param ... optional parameters for ComplexHeatmap::Heatmap().
#' @return A heatmap of region set scores across. Each row is a region set,
#' each column is a PC. The color corresponds to the relative rank of a 
#' region set's score for a given PC out of all tested region sets.
#'
#' @examples data("rsScores")
#' scoreHeatmap <- rsScoreHeatmap(rsScores, 
#'           PCsToAnnotate=paste0("PC", 1:2), orderByPC = "PC2")
#' @export

rsScoreHeatmap <- function(rsScores, PCsToAnnotate=paste0("PC", 1:5),
                            orderByPC="PC1", rsNameCol = "rsName", topX = 20, 
                           col=c("red", "#EEEEEE", "blue"),
                           row_title = "Region Set", column_title = "Principal Component",
                           column_title_side = "bottom",
                           cluster_rows = FALSE, cluster_columns = FALSE,
                           show_row_names = TRUE, 
                           row_names_max_width = unit(10000, "mm"),
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
        setorderv(rsScores, PCsToAnnotate[i], order = -1L, na.last=TRUE) # descending order
        rsScores[, PCsToAnnotate[i] := seq_len(rsNum)]
        
        # center around zero
        # rsScores[, PCsToAnnotate[i] := ((rsNum + 1) / 2) - get(PCsToAnnotate[i])]
    }
    
    # heatmap of the centered ranks
    # back to first order, -1L means decreasing order
    setorderv(rsScores, orderByPC, order = 1L) 
    rowNames <-  rsScores[, get(rsNameCol)] # redefined/reordered later
    row.names(rsScores) <- rowNames
    rsScores[, c(rsNameCol) := NULL]
    rsScores <- as.matrix(rsScores)
    row.names(rsScores) <- rowNames
    Heatmap(rsScores[seq_len(topX), ], 
            col = col,
            row_title = row_title, column_title = column_title, 
            cluster_rows = cluster_rows, 
            cluster_columns = cluster_columns, 
            column_title_side = column_title_side,
            show_row_names = show_row_names, 
            row_names_max_width = row_names_max_width, 
            name = name,
            heatmap_legend_param=, ...)
    
}


#' Visualize how individual regions are associated with principal components
#' 
#' Visualize how much each region in a region set is associated with each PC.
#' For each PC, the average absolute loading is calculated for 
#' each region in the region set. Then for a given PC, 
#' the average loading is converted to a percentile/quantile based 
#' on the distribution of all loadings for that PC. These values are
#' plotted in a heatmap.
#' 
#' @param loadingMat matrix of loadings (the coefficients of 
#' the linear combination that defines each PC). One named column for each PC.
#' One row for each original dimension/variable (should be same order 
#' as original data/signalCoord). The x$rotation output of prcomp().
#' @param signalCoord a GRanges object or data frame with coordinates 
#' for the genomic signal/original data (eg DNA methylation) 
#' included in the PCA. Coordinates should be in the 
#' same order as the original data and the loadings 
#' (each item/row in signalCoord
#' corresponds to a row in loadingMat). If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.
#' @param regionSet A genomic ranges object with regions corresponding
#' to the same biological annotation. These are the regions that will be 
#' visualized. Must be from the same reference genome
#' as the coordinates for the actual data (signalCoord).
#' @param rsName character vector. Names of the region sets in the same
#' order as GRList. For use as a title for each heatmap.
#' @param PCsToAnnotate A character vector with principal components to 
#' include. eg c("PC1", "PC2") These should be column names of loadingMat.
#' @param maxRegionsToPlot how many top regions from region set to include
#' in heatmap. Including too many may slow down computation and increase memory
#' use. If regionSet has more regions than maxRegionsToPlot, a number of regions 
#' equal to maxRegionsToPlot will be randomly sampled from the region set and
#' these regions will be plotted. Clustering rows is a major limiting factor
#' on how long it takes to plot the regions so if you want to plot many regions, 
#' you can also set cluster_rows to FALSE.
#' @param row_title character object, row title
#' @param column_title character object, column title
#' @param column_title_side character object, where to put the column title:
#' "top" or "bottom"
#' @param cluster_rows "logical" object, whether to cluster rows or not (may 
#' increase computation time significantly for large number of rows)
#' @param cluster_columns "logical" object, whether to cluster columns. 
#' It is recommended
#' to keep this as FALSE so it will be easier to compare PCs 
#' (with cluster_columns = FALSE, they will be in the same specified
#' order in different heatmaps)
#' @param name character object, legend title
#' @param col a vector of colors or a color mapping function which
#' will be passed to the ComplexHeatmap::Heatmap() function. See ?Heatmap
#' (the "col" parameter) for more details.
#' @param ... optional parameters for ComplexHeatmap::Heatmap()
#' @return a heatmap. Columns are PCs, rows are regions. 
#' This heatmap allows you to see if some regions are 
#' associated with certain PCs but not others. Also, you can see if a subset of 
#' regions in the region set are associated with PCs while another subset
#' are not associated with any PCs 
#' To color each region, first the absolute loading 
#' values within that region are
#' averaged. Then this average is compared to the distribution of absolute
#' loading values for all individual genomic signal values to get 
#' a quantile/percentile 
#' for that region. Colors are based on this quantile/percentile. 
#' The output is a Heatmap object (ComplexHeatmap package).
#' 
#' @examples data("brcaLoadings1")
#' data("brcaMCoord1")
#' data("esr1_chr1")
#' data("brcaPCScores")
#' regionByPCHM <- regionQuantileByPC(loadingMat = brcaLoadings1, 
#'                                    signalCoord = brcaMCoord1, 
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
regionQuantileByPC <- function(loadingMat, signalCoord, regionSet, 
                               rsName = "", PCsToAnnotate=paste0("PC", 1:5),
                               maxRegionsToPlot = 8000, cluster_rows = TRUE, 
                               row_title = "Region", column_title = rsName,
                               column_title_side = "top",
                               cluster_columns = FALSE, 
                               name = "Percentile of Loading Scores in PC", 
                               col = c("skyblue", "yellow"), ...) {
    
    ################### checking inputs  #################################
    
    ########## check that inputs are the correct class
    # exports coordinateDT to this environment (converts signalCoord)
    checkConvertInputClasses(loadingMat=loadingMat,
                             signalCoord=signalCoord,
                             regionSet=regionSet,
                             PCsToAnnotate = PCsToAnnotate)
    
    ########## check that dimensions of inputs are consistent
    # length of signal coord = nrow of loadingMat
    if (nrow(coordinateDT) != nrow(loadingMat)) {
        stop(cleanws("The number of coordinates in 
            signalCoord (length(signalCoord)) does not equal the number of 
                     rows in loadingMat"))
    } 
    
    ######### check that appropriate columns are present
    # PCsToAnnotate are column names of loadingMat
    if (!all(PCsToAnnotate %in% colnames(loadingMat))) {
        missingCols = PCsToAnnotate[!(PCsToAnnotate %in% colnames(loadingMat))]
        stop(cleanws(paste0("Some PCsToAnnotate are not 
                            columns of loadingMat: ", missingCols)))
    }
    
    #######
    # what happens if there are NAs or Inf in loadingMat?
    
    #################################################################

    
    # if too many regions to plot, randomly subsample regions
    if (length(regionSet) > maxRegionsToPlot) {
        regionInd <- sample(x = seq_along(regionSet), size = maxRegionsToPlot, 
                            replace = FALSE)
        # get subset
        regionSet <- regionSet[regionInd]
    }
    
    
    
    rsRegionAverage <- averagePerRegion(loadingMat = loadingMat, 
                                       signalCoord =coordinateDT, 
                                       regionSet = regionSet, 
                                       PCsToAnnotate = PCsToAnnotate,
                                       returnQuantile = TRUE)
    # ranking in terms of percentiles in case there were different 
    # distributions of loading scores for each PC
    
    # the heatmap
    Heatmap(matrix = as.matrix(rsRegionAverage[, PCsToAnnotate, with=FALSE]), 
            row_title = row_title,
            column_title = rsName,
            column_title_side = column_title_side,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_columns, 
            name = name, 
            col = col, ...)
}

