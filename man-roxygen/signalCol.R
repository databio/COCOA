#' @param signalCol A character vector with the names of the sample variables
#' of interest/target variables (e.g. PCs or sample phenotypes). 

#' <%=ifelse(exists("usesRSScores") && usesRSScores, "Must be column names of rsScores.", "") %>

#' <%=ifelse(exists("usesTargetVar") && usesTargetVar, "The columns in `sampleLabels` for which to calculate", "") %>
#' <%=ifelse(exists("usesTargetVar") && usesTargetVar, "the variation related to the epigenetic data", "") %>
#' <%=ifelse(exists("usesTargetVar") && usesTargetVar, "(e.g. correlation) and then to run COCOA on.", "") %>

#' <%=ifelse(exists("isRSRankingIndex") && isRSRankingIndex, "The columns in rsScores for which you want", "") %>
#' <%=ifelse(exists("isRSRankingIndex") && isRSRankingIndex, "the indices of the original region sets.", "") %>





# (aggregateSignal)
# (aggregateSignalGRList)
# (getMetaRegionProfile)
# (getTopRegions)
# (rsScoreHeatmap) usesRSScores
# (rsRankingIndex) usesRSScores isRSRankingIndex
# (regionQuantileByTargetVar)
# (runCOCOAPerm) usesTargetVar
# (corPerm) usesTargetVar
# (getGammaPVal) usesRSScores
# (getPermStat) usesRSScores