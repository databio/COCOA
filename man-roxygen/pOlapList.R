#' @param pOlapList list. This parameter is only used if the scoring metric is
#' "proportionWeightedMean" and olList is also provided as an argument. Each
#' item of the list should be a vector that contains the proportion overlap 
#' between signalCoord and regions from one region set (one item of GRList).
#' Specifically, each value should be the proportion of the region set region 
#' that is overlapped
#' by a signalCoord region.
#' The proportion overlap values should be in the same order as the overlaps
#' given by olList for the corresponding region set.  