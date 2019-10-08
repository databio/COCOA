#' @param absVal logical. If TRUE, take the absolute value of values in
#' signal. Choose TRUE if you think there may be some 
#' genomic loci in a region set that will increase and others
#' will decrease (if there may be anticorrelation between
#' regions in a region set). Choose FALSE if you expect regions in a 
#' given region set to all change in the same direction (all be positively
#' correlated with each other).