#' @param rsScores data.frame. A data.frame with region set
#' scores. The output of the 'aggregateSignalGRList' function.
#' Each row is a region set. One column for each sample
#' variable of interest (e.g. PC or sample phenotype).
#' Also can have columns with info on the overlap between the 
#' region set and the epigenetic data. 
#' Rows should be in the same order as the region sets in GRList
#' (the list of region sets used to create rsScores.)
# (runCOCOAPerm) Must include columns with names given by 'colsToAnnotate'.
# (getPermStat)
# (rsRankingIndex)
# (rsScoreHeatmap)
# (getGammaPVal)

