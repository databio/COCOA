################################################################################
# functions to visualize results of COCOA, relevant regions, 
# and variation in the dataset
###############################################################################
# Some imported functions:
# data.table: copy, :=, setorder
# MIRA: dtToGr
# ComplexHeatmap
###

# functions to visualize results of COCOA, relevant regions, 
# and variation in the dataset
# 

# color schemes: red/blue, yellow/red, red/grey, skyblue/coral, skyblue/yellow

#' Visualize how genomic signal in a region set changes 
#' along a given axis
#' 
#' Look at genomic signal (e.g., DNA methylation values) in regions of 
#' interest across samples, 
#' with samples ordered according to a variable of interest (e.g. PC score). 
#' The ComplexHeatmap package
#' is used and additional parameters for the ComplexHeatmap::Heatmap function
#' may be passed to this function to modify the heatmap.   
#' 
#' @templateVar requireSampleNames TRUE
#' @template genomicSignal
#' @template signalCoord
#' @templateVar refGenomeWarning TRUE
#' @templateVar rsVisualization TRUE
#' @template regionSet
#' @param sampleScores A matrix. Must contain a column for the 
#' variable of interest/target variable.
#' E.g. The variable of interest could be 
#' the principal component scores for the samples. 
#' `sampleScores` must have sample names/IDs as row names,
#' These same sample names must be column names of genomicSignal.
#' @param orderByCol A character object. A variable to order samples by
#' (order rows of heatmap by variable, from high to low value).
#' Must be the name of a column in sampleScores. For instance, if doing
#' unsupervised COCOA with PCA, orderByCol might be the name of one of the PCs
#' (e.g. "PC1"). If doing supervised COCOA, orderByCol might be the name
#' of the target variable of the supervised analysis. 
#' @param topXVariables Numeric. The number of variables from genomicSignal
#' to plot. The variables with the highest scores according to 
#' variableScores will be plotted. Can help to reduce the size of the plot.
#' @param variableScores Numeric. A vector that has a numeric score for
#' each variable in genomicSignal (length(variableScores) should equal
#' nrow(genomicSignal)). Only used if topXVariables is given. The highest
#' `topXVariables` will be plotted.
#' @param decreasing Logical. Whether samples should be sorted in 
#' decreasing order of `orderByCol` or not (FALSE is increasing order).
#' @param cluster_rows Logical. Whether rows should be clustered. 
#' This should be kept as FALSE to keep the correct ranking of 
#' samples/observations according to their PC score.
#' @param cluster_columns Logical. Whether to cluster columns 
#' (the genomic signal, e.g. DNA methylation values for each CpG).
#' @param row_title Character object, row title
#' @param column_title Character object, column title
#' @param column_title_side Character object, where to put the column title:
#' "top" or "bottom"
#' @param name Character object, legend title
#' @param col A vector of colors or a color mapping function which
#' will be passed to the ComplexHeatmap::Heatmap() function. See ?Heatmap
#' (the "col" parameter) for more details. "#EEEEEE" is the code for a
#' color similar to white.
#' @param ... Optional parameters for ComplexHeatmap::Heatmap()
#' @return A heatmap of genomic signal values (eg DNA methylation levels) 
#' in regions of interest (regionSet), with rows ordered by the
#' column of sampleScores given with `orderByCol`.
#' Each row is a patient/sample and 
#' each column is an individual genomic signal value. 
#'
#' @examples data("brcaMethylData1")
#' data("brcaMCoord1")
#' data("esr1_chr1")
#' data("brcaPCScores")
#' signalHM <- signalAlongAxis(genomicSignal=brcaMethylData1,
#'                              signalCoord=brcaMCoord1,
#'                              regionSet=esr1_chr1,
#'                              sampleScores=brcaPCScores,
#'                              orderByCol="PC1", cluster_columns=TRUE)
#' 
#' @export
# previously called rsMethylHeatmap
signalAlongAxis <- function(genomicSignal, signalCoord, regionSet, 
                          sampleScores, orderByCol="PC1", topXVariables=NULL, 
                          variableScores=NULL,
                          decreasing=TRUE,
                          cluster_columns = FALSE, 
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
    
    if (!(is(sampleScores, "matrix") || is(sampleScores, "data.frame"))) {
        stop("sampleScores should be a matrix or data.frame.")
    }
    
    
    # PCA object must have subject_ID as row.names (corresponding 
    # to column names of genomicSignal)
    if (sum(row.names(sampleScores) %in% colnames(genomicSignal)) < 2) {
        stop(cleanws("Sample names on pca data (row names) 
                     must match sample names on methylation
                     (column names)"))
    }
    
    
    if (!is(regionSet, "GRanges")) {
        stop("regionSet should be a GRanges object. Check object class.")
    }
    
    if (!is.null(topXVariables)) {
        if (is.null(variableScores)) {
            stop("To plot the topXVariables, variableScores must be given.")
        }
        if (!(length(variableScores) == nrow(genomicSignal))) {
            stop("length(variableScores) should equal nrow(genomicSignal)")
        }
        
    }
    
    
    
    # coordGR =
    olList <- findOverlaps(regionSet, coordGR)
    # regionHitInd <- sort(unique(queryHits(olList)))
    cytosineHitInd <- sort(unique(subjectHits(olList)))
    thisRSMData <- t(genomicSignal[cytosineHitInd, ])
    nRegion = length(unique(queryHits(olList)))
    # get top variables
    if (!is.null(topXVariables)) {
        if (nRegion > topXVariables) {
            # select top variables
            thisRSMData <- thisRSMData[, order(variableScores[cytosineHitInd], decreasing = TRUE)][, 1:topXVariables]
        }
    }
    subject_ID <- row.names(thisRSMData)
    # centeredPCAMeth <- t(apply(t(genomicSignal), 1, 
    #                            function(x) x - pcaData$center)) #center first 
    # reducedValsPCA <- centeredPCAMeth %*% pcaData$rotation
    # reducedValsPCA <- pcaData$x
    # pcaData must have subject_ID as row name
    thisRSMData <- thisRSMData[names(sort(sampleScores[, orderByCol], 
                                          decreasing = decreasing)), ]
    
    # message(paste0("Number of epigenetic features: ", length(cytosineHitInd)))
    # message(paste0("Number of region set regions: ", nRegion))
    
    
    
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
#' Heatmap of the ranking of region set scores across target variables.
#' A visualization of the rank of region sets in each target variable, 
#' allowing the
#' user to see if a region set is ranked highly for all target variables 
#' or only a subset.
#' Region sets will be ranked from highest scoring to lowest based on 
#' their score for `orderByCol`.
#' The ComplexHeatmap package
#' is used and additional parameters for the ComplexHeatmap::Heatmap function
#' may be passed to this function to modify the heatmap.  
#' 
#' @template rsScores
#' @templateVar usesRSScores TRUE
#' @template signalCol
#' @param orderByCol A character object. Target variable to order by in heatmap 
#' (arranged in decreasing order for scores so p values should 
#' be -log transformed). Must be the name of a column in rsScores.
#' @param rsNameCol Character. Name of the column in rsScores that has the
#' names/identifiers for the region sets so these can be included 
#' in the plot as row names.
#' @param topX Number of top region sets to include in the heatmap
#' @param row_title Character object, row title
#' @param column_title Character object, column title
#' @param column_title_side Character object, where to put the column title:
#' "top" or "bottom"
#' @param cluster_rows Logical object, whether rows should be clustered. 
#' This should be kept as FALSE to keep the correct ranking of region sets.
#' @param cluster_columns Logical object, whether to cluster columns. 
#' It is recommended
#' to keep this as FALSE so it will be easier to compare target variables
#' that are ordered (such as principal components). 
#' With cluster_columns = FALSE, they will be in the same specified
#' order in different heatmaps.
#' @param show_row_names Logical object, display row names (ie region set names)
#' @param row_names_max_width "unit" object. The amount of room to 
#' allocate for row names. See ?grid::unit for object type.
#' @param name Character object, legend title
#' @param col A vector of colors or a color mapping function which
#' will be passed to the ComplexHeatmap::Heatmap() function. See ?Heatmap
#' (the "col" parameter) for more details. "#EEEEEE" is the code for a
#' color similar to white.
#' @param ... Optional parameters for ComplexHeatmap::Heatmap().
#' @return A heatmap of region set scores across. Each row is a region set,
#' each column is a target variable. 
#' The color corresponds to the relative rank of a 
#' region set's score for a given target variable out of all tested region sets.
#'
#' @examples data("rsScores")
#' scoreHeatmap <- rsScoreHeatmap(rsScores, 
#'           signalCol=paste0("PC", 1:2), orderByCol = "PC2")
#' @export

rsScoreHeatmap <- function(rsScores, signalCol=paste0("PC", 1:5),
                            orderByCol="PC1", rsNameCol = "rsName", topX = 20, 
                           col=c("red", "#EEEEEE", "blue"),
                           row_title = "Region Set", column_title = "Principal Component",
                           column_title_side = "bottom",
                           cluster_rows = FALSE, cluster_columns = FALSE,
                           show_row_names = TRUE, 
                           row_names_max_width = unit(10000, "mm"),
                           name="Rank within PC", ...) {
    
    
    if (!is(signalCol, "character")) {
        stop("signalCol should be a character object (eg 'PC1').")
    }
    
    # function is not meant for making a heatmap of a single PC but should work
    if (is(rsScores, "numeric")) {
        warning("rsScores should be a data.frame")
        rsScores <- data.frame(rsScores)
        colnames(rsScores) <- orderByCol
    }
    if (!is(rsScores, "data.frame")) {
        stop("rsScores should be a data.frame.")
    }
    
    if (!(orderByCol %in% colnames(rsScores))) {
        stop(cleanws(paste0("orderByCol parameter:", orderByCol, 
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
    signalCol <- signalCol[signalCol %in% colnames(rsScores)]
    if (length(signalCol) == 0) {
        stop("Please check format of PC names in signalCol.")
    }
     
    # apparently erases row names
    rsScores <- rsScores[, c(signalCol, rsNameCol), with=FALSE] 

    
    # how to deal with NA?
    
    # number of region sets tested
    rsNum <- nrow(rsScores)
    
    # convert to data.table to do some data.table operations
    rsScores = as.data.table(rsScores)
    
    for (i in seq_along(signalCol)) {
        # first convert to rank
        setorderv(rsScores, signalCol[i], order = -1L, na.last=TRUE) # descending order
        rsScores[, signalCol[i] := seq_len(rsNum)]
        
        # center around zero
        # rsScores[, signalCol[i] := ((rsNum + 1) / 2) - get(signalCol[i])]
    }
    
    # heatmap of the centered ranks
    # back to first order, -1L means decreasing order
    setorderv(rsScores, orderByCol, order = 1L) 
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


#' Visualize how individual regions are associated with target variable
#' 
#' Visualize how much each region in a region set 
#' is associated with each target variable. 
#' For each target variable (`signalCol`), the average (absolute) 
#' signal value is calculated for 
#' each region in the region set. Then for a given target variable, 
#' the average signal is converted to a percentile/quantile based 
#' on the distribution of all signal values 
#' for that target variable. These values are
#' plotted in a heatmap.
#' 
#' @template signal
#' @template signalCoord
#' @templateVar refGenomeWarning TRUE
#' @templateVar rsVisualization TRUE
#' @template regionSet
#' @param rsName Character vector. Names of the region sets in the same
#' order as GRList. For use as a title for each heatmap.
#' @template signalCol
#' @param maxRegionsToPlot How many top regions from region set to include
#' in heatmap. Including too many may slow down computation and increase memory
#' use. If regionSet has more regions than maxRegionsToPlot, a number of regions 
#' equal to maxRegionsToPlot will be randomly sampled from the region set and
#' these regions will be plotted. Clustering rows is a major limiting factor
#' on how long it takes to plot the regions so if you want to plot many regions, 
#' you can also set cluster_rows to FALSE.
#' @param row_title Character object, row title
#' @param column_title Character object, column title
#' @param column_title_side Character object, where to put the column title:
#' "top" or "bottom"
#' @param cluster_rows Logical object, whether to cluster rows or not (may 
#' increase computation time significantly for large number of rows)
#' @param cluster_columns Logical object, whether to cluster columns. 
#' It is recommended
#' to keep this as FALSE so it will be easier to compare PCs 
#' (with cluster_columns = FALSE, they will be in the same specified
#' order in different heatmaps)
#' @param name Character object, legend title
#' @param col A vector of colors or a color mapping function which
#' will be passed to the ComplexHeatmap::Heatmap() function. See ?Heatmap
#' (the "col" parameter) for more details.
#' @template absVal
#' @param ... Optional parameters for ComplexHeatmap::Heatmap()
#' @return A heatmap. Columns are signalCol's, rows are regions. 
#' This heatmap allows you to see if some regions are 
#' associated with certain target variables but not others. 
#' Also, you can see if a subset of 
#' regions in the region set are associated with 
#' target variables while another subset
#' are not associated with any target variables 
#' To color each region, first the (absolute) signal 
#' values within that region are
#' averaged. Then this average is compared to the distribution of all (absolute)
#' individual signal values for the given target variable to get 
#' a quantile/percentile 
#' for that region. Colors are based on this quantile/percentile. 
#' The output is a Heatmap object (ComplexHeatmap package).
#' 
#' @examples
#' data("brcaATACCoord1")
#' data("brcaATACData1")
#' data("esr1_chr1")
#' featureContributionScores <- prcomp(t(brcaATACData1))$rotation
#' regionByPCHM <- regionQuantileByTargetVar(signal = featureContributionScores, 
#'                                    signalCoord = brcaATACCoord1, 
#'                                    regionSet = esr1_chr1, 
#'                                    rsName = "Estrogen Receptor Chr1", 
#'                                    signalCol=paste0("PC", 1:2),
#'                                    maxRegionsToPlot = 8000, 
#'                                    cluster_rows = TRUE, 
#'                                    cluster_columns = FALSE, 
#'                                    column_title = rsName, 
#'                                    name = "Percentile of Loading Scores in PC")
#'                                    
#' 
#' @export
regionQuantileByTargetVar <- function(signal, signalCoord, regionSet, 
                               rsName = "", signalCol=paste0("PC", 1:5),
                               maxRegionsToPlot = 8000, cluster_rows = TRUE, 
                               row_title = "Region", column_title = rsName,
                               column_title_side = "top",
                               cluster_columns = FALSE, 
                               name = "Percentile of Loading Scores in PC", 
                               col = c("skyblue", "yellow"),
                               absVal=TRUE, ...) {
    
    ################### checking inputs  #################################
    
    ########## check that inputs are the correct class
    checkConvertInputClasses(signal=signal,
                             signalCoord=signalCoord,
                             regionSet=regionSet,
                             signalCol = signalCol)
    
    ########## check that dimensions of inputs are consistent
    # length of signal coord = nrow of `signal`
    if (length(signalCoord) != nrow(signal)) {
        stop(cleanws("The number of coordinates in 
            signalCoord (length(signalCoord)) does not equal the number of 
                     rows in `signal`"))
    } 
    
    ######### check that appropriate columns are present
    # signalCol are column names of `signal`
    if (!all(signalCol %in% colnames(signal))) {
        missingCols = signalCol[!(signalCol %in% colnames(signal))]
        stop(cleanws(paste0("Some signalCol are not 
                            columns of `signal`: ", missingCols)))
    }
    
    #######
    # what happens if there are NAs or Inf in `signal`?
    
    #################################################################

    
    # if too many regions to plot, randomly subsample regions
    if (length(regionSet) > maxRegionsToPlot) {
        regionInd <- sample(x = seq_along(regionSet), size = maxRegionsToPlot, 
                            replace = FALSE)
        # get subset
        regionSet <- regionSet[regionInd]
    }
    
    
    
    rsRegionAverage <- averagePerRegion(signal = signal, 
                                       signalCoord =signalCoord, 
                                       regionSet = regionSet, 
                                       signalCol = signalCol,
                                       returnQuantile = TRUE,
                                       absVal=absVal)
    # ranking in terms of percentiles in case there were different 
    # distributions of loading scores for each PC
    
    # the heatmap
    Heatmap(matrix = as.matrix(rsRegionAverage[, signalCol, with=FALSE]), 
            row_title = row_title,
            column_title = rsName,
            column_title_side = column_title_side,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_columns, 
            name = name, 
            col = col, ...)
}

