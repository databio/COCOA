#' @param olList list. Each list item should be a "SortedByQueryHits" object 
#' (output of findOverlaps function). Each hits object should have the overlap
#' information between signalCoord and one item of GRList (one unique region set).
#' The region sets from GRList must be the "subject" in findOverlaps 
#' and signalCoord must be the "query". E.g. findOverlaps(subject=regionSet,
#' query=signalCoord).
#' Providing this information can greatly improve permutation speed since the 
#' overlaps will not have to be calculated for each permutation. 
#' The "runCOCOAPerm" function calculates this information only once, internally,
#' so this does not have to be provided when using that function. When using 
#' this parameter, signalCoord, 
#' genomicSignal, and each region set must be in the same order as they were
#' when olList was created. Otherwise, the wrong genomic loci will be referenced
#' (e.g. if epigenetic features were filtered out of genomicSignal after olList
#' was created.) 