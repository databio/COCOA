#' @param signalCoord a GRanges object or data.frame with coordinates 
#' for the genomic signal/original data (e.g. DNA methylation) 
#' included in the PCA. Coordinates should be in the 
#' same order as the original data and the loadings 
#' (each item/row in signalCoord
#' corresponds to a row in `signal`). If a data.frame, 
#' must have chr and start columns. If end is included, start 
#' and end should be the same. Start coordinate will be used for calculations.