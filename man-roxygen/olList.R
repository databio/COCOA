#' @param olList list of "SortedByQueryHits" objects (output of findOverlaps function) 
#' containing overlap information between signalCoord and GRList.
#' Providing this information can significantly speed up permutation calculations in 
#' the "runCOCOAPerm" function. Make sure to maintain the original order of 
#' signalCoord, genomicSignal, and region sets in olList.
