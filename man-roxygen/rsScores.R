#' @param rsScores data.frame containing region set scores from 
#' the 'aggregateSignalGRList' function. Each row represents a region set, 
#' and columns contain sample variables of interest (e.g., PC or sample phenotype).
#' Additionally, it can have columns with overlap information between the region set 
#' and epigenetic data.
#' Make sure that the rows are in the same order as the region sets in GRList.
#' (the list of region sets used to create rsScores.)
# (runCOCOAPerm) Must include columns with names given by 'colsToAnnotate'.
# (getPermStat)
# (rsRankingIndex)
# (rsScoreHeatmap)
# (getGammaPVal)

