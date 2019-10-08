#' @param signalCoordType character. Can be "default", "singleBase", or 
#' "multiBase". This describes whether the coordinates for `signal` 
#' (`signalCoord`) are each a single base (e.g. as for DNA methylation)
#' or a region/multiple bases (e.g. as for ATAC-seq). Different scoring
#' options are available for each type of data. If "default" is given,
#' the type of coordinates will be detected automatically. For "default", if each
#' coordinate start value equals the coordinate end value 
#' (all(start(signalCoord) == end(signalCoord))), "singleBase"
#' will be used. Otherwise, "multiBase" will be used.