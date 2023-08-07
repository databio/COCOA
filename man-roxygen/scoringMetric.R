#' <%=ifelse(exists("usesAggrMethod") && usesAggrMethod, "@param aggrMethod character. A character object with the aggregation method.", "") %>
#' <%=ifelse(exists("usesAggrMethod") && usesAggrMethod, "Similar to aggregateSignalGRList `scoringMetric` parameter.", "@param scoringMetric A character object with the scoring metric.") %>
#' There are different methods available for 
#' signalCoordType="singleBase" vs  signalCoordType="multiBase".
#' For "singleBase", the available methods are "regionMean", 
#' "regionMedian", "simpleMean", and "simpleMedian". 
#' The default method is "regionMean".
#' For "multiBase", the methods are "proportionWeightedMean", 
#' "simpleMean", and "simpleMedian". The default is "proportionWeightedMean".
#' "regionMean" is a weighted
#' average of the signal, weighted by region (absolute value of signal 
#' if absVal=TRUE). The methods for averaging signal values within region sets are:
#' "regionMean": Average the signal within each regionSet region, then average all regions.
#' "regionMedian": Take the median within each regionSet region, then median across all regions.
#' "simpleMean": Unweighted average of all overlapping signal values for each regionSet region.
#' "simpleMedian": Take the median of all overlapping signal values for each regionSet region.
#' "proportionWeightedMean": Weighted average of signal values, considering the proportion 
#' of each signalCoord region that overlaps with a regionSet region. The weights are 
#' based on the proportion overlaps.
