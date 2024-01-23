#' @param rsOL "SortedByQueryHits" object (output of findOverlaps function) 
#' containing overlap information between signalCoord and one region set (GRList).
#' The region set must be the "subject" in findOverlaps, and signalCoord must be the "query".
#' Providing this information significantly improves permutation speed since overlaps do not 
#' needs to be recalculated for each permutation.
#' Make sure to maintain the original order of signalCoord, genomicSignal, and the region set 
#' when using this parameter, as altering the order may lead to incorrect genomic loci references.
