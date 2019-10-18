#' <%=ifelse(exists("usesAggrMethod") && usesAggrMethod, "@param aggrMethod character. A character object with the aggregation method.", "") %>
#' <%=ifelse(exists("usesAggrMethod") && usesAggrMethod, "Similar to runCOCOA `scoringMetric` parameter.", "@param scoringMetric A character object with the scoring metric.") %>
#' There are different methods available for 
#' signalCoordType="singleBase" vs  signalCoordType="multiBase".
#' For "singleBase", the available methods are "regionMean" and 
#' "simpleMean". The default method is "regionMean".
#' For "multiBase", the methods are "proportionWeightedMean" and 
#' "simpleMean". The default is "proportionWeightedMean".
#' "regionMean" is a weighted
#' average of the signal, weighted by region (absolute value of signal 
#' if absVal=TRUE). First the signal is
#' averaged within each regionSet region, 
#' then all the regions are averaged. With
#' "regionMean" method, be cautious in interpretation for
#' region sets with low number of regions that overlap signalCoord. 
#' The "simpleMean"
#' method is just the unweighted average of all (absolute) signal values that
#' overlap the given region set. For multiBase data, this includes
#' signal regions that overlap a regionSet region at all (1 base
#' overlap or more) and the signal for each overlapping region is
#' given the same weight for the average regardless of how much it overlaps. 
#' "proportionWeightedMean" is a weighted average of all signalCoord 
#' regions that overlap with regionSet regions. For each signalCoord region
#' that overlaps with a regionSet region, we calculate what proportion
#' of the regionSet region is covered. Then this proportion is used to
#' weight the signal value when calculating the mean. 
#' The denominator of the mean
#' is the sum of all the proportion overlaps.