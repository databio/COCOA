#' @param sampleLabels Matrix or data.frame. Rows should be samples, 
#' columns should be the target variables 
#' (whatever variable you want to test for association with
#' the epigenetic signal: e.g. PC scores),
#' all columns in featureMat will be used (subset when passing to function
#' in order to not use all columns).

#' @param sampleLabels data.frame/matrix. Sample labels/values that 
#' you are running COCOA to find region sets associated with. These 
#' values will be shuffled for the permutation test. Rows are samples.
#' Each column is a sample label.


# (corPerm)
# (runCOCOAPerm)