#' @param scoringMetric A character object with the scoring metric.
#' There are different scoring metrics available for 
#' signalCoordType="singleBase" vs  signalCoordType="multiBase".
#' For "singleBase", the available scoring methods are "regionMean", 
#' "simpleMean", and "rankSum". The default method is "regionMean".
#' For "multiBase", the scoring methods are "proportionWeightedMean" and 
#' "simpleMean". The default is "proportionWeightedMean".
#' "regionMean" is a weighted
#' average of the signal, weighted by region (absolute value of signal 
#' if absVal=TRUE). First the signal is
#' averaged within each regionSet region, 
#' then all the regions are averaged. With
#' "regionMean" score, be cautious in interpretation for
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
#' Wilcoxon rank sum test ("rankSum") is also supported but is
#' skewed toward ranking large region sets highly and is
#' significantly slower than the "regionMean" method. 
#' For the ranksum method, the absolute loadings for loadings that
#' overlap the given region set are taken as a group and all the
#' loadings that do not overlap the region set are taken as
#' the other group. Then p value is then given as the score.
#' It is a one sided test, with the alternative hypothesis
#' that the loadings in the region set will be greater than
#' the loadings not in the region set.