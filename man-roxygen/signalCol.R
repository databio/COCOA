#' @param signalCol A character vector with the sample variables
#' of interest (e.g. PCs or sample phenotypes). 
#' <%=ifelse(exists("isRSRankingIndex") && isRSRankingIndex, "The columns in rsScores for which you want", "") %>
#' <%=ifelse(exists("isRSRankingIndex") && isRSRankingIndex, "the indices of the original region sets.", "") %>
#' <%=ifelse(exists("usesRSScores") && usesRSScores, "Must be column names of rsScores.", "") %>
#' e.g. c("PC1", "PC2") 

# (aggregateSignal)
# (runCOCOA)
# (getMetaRegionProfile)
# (getTopRegions)
# (rsScoreHeatmap) usesRSScores
# (rsRankingIndex) usesRSScores isRSRankingIndex
# (regionQuantileByVOI)