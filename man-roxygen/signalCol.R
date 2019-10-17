#' @param signalCol A character vector with the names of the sample variables
#' of interest/target variables (e.g. PCs or sample phenotypes). 

#' <%=ifelse(exists("usesRSScores") && usesRSScores, "Must be column names of rsScores.", "") %>

#' <%=ifelse(exists("usesSampleLabels") && usesSampleLabels, "The columns in `sampleLabels` for which to calculate", "") %>
#' <%=ifelse(exists("usesSampleLabels") && usesSampleLabels, "the variation related to the epigenetic data", "") %>
#' <%=ifelse(exists("usesSampleLabels") && usesSampleLabels, "(e.g. correlation) and then to run COCOA on.", "") %>

#' <%=ifelse(exists("isRSRankingIndex") && isRSRankingIndex, "The columns in rsScores for which you want", "") %>
#' <%=ifelse(exists("isRSRankingIndex") && isRSRankingIndex, "the indices of the original region sets.", "") %>





# (aggregateSignal)
# (runCOCOA)
# (getMetaRegionProfile)
# (getTopRegions)
# (rsScoreHeatmap) usesRSScores
# (rsRankingIndex) usesRSScores isRSRankingIndex
# (regionQuantileByTargetVar)
# (runCOCOAPerm) usesSampleLabels
# (corPerm) usesSampleLabels
# (getGammaPVal) usesRSScores
# (getPermStat) usesRSScores