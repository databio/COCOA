#' @param regionSet A genomic ranges (GRanges) object with regions corresponding
#' to the same biological annotation.
#' <%=ifelse(exists("refGenomeWarning") && refGenomeWarning, "Must be from the same reference genome as the coordinates for the actual data/samples (signalCoord).", "") %>
#' <%=ifelse(exists("rsVisualization") && rsVisualization, "The regions that will be visualized.", "") %>