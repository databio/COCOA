#' @param signal Matrix of feature contribution scores (the contribution of 
#' each epigenetic feature to each target variable). One named column for each 
#' target variable.
#' One row for each original epigenetic feature (should be same order 
#' as original data/signalCoord). For (an unsupervised) example, if PCA was
#' done on epigenetic data and the
#' goal was to find region sets associated with the principal components, you 
#' could use the x$rotation output of prcomp(epigenetic data) as the
#' feature contribution scores/`signal` parameter.