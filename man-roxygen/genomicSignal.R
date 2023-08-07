#' @param genomicSignal Matrix/data.frame. 
#' The genomic signal.
#' Columns of genomicSignal should be samples/patients. 
#' Rows should be individual signal/features
#' <%=ifelse(exists("requireSampleNames") && requireSampleNames, "Must have sample names/IDs as column names,", "") %>
#' <%=ifelse(exists("requireSampleNames") && requireSampleNames, "These same sample names must be row names of sampleScores.", "") %>
# (runCOCOAPerm)
# (corPerm)
# (signalAlongAxis) requireSampleNames
