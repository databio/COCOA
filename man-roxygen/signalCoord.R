#' @param signalCoord A GRanges object or data frame with coordinates 
#' for the genomic signal/original epigenetic data. 
#' Coordinates should be in the 
#' same order as the original data and the feature contribution scores 
#' (each item/row in signalCoord
#' corresponds to a row in signal). If a data.frame, 
#' must have chr and start columns (optionally can have end column, 
#' depending on the epigenetic data type).