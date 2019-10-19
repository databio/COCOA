#' @param targetVar Matrix or data.frame. Rows should be samples. 
#' Columns should be the target variables 
#' (whatever variable you want to test for association with
#' the epigenetic signal: e.g. PC scores),
# all columns in featureMat will be used (subset when passing to function
# in order to not use all columns).

# (corPerm)
# (runCOCOAPerm) values will be shuffled for the permutation test.
# signalAlongAxis has @param sampleScores
