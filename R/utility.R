# Utility functions such as format conversion

# Check whether function inputs are the right class
# if they are wrong class but close enough, convert them.
# if they are wrong class but not convertible, give an error.
# Having a function to do this means, I don't need to write the 
# same code within each exported function and can change it in one place.
# default NULL means do not check a given parameter

# @return Does not teturn anything but optionally can assign
# objects to the environment that called it. Can assign `signalCoord`,
# and `GRList` if they need to be converted to a different object class.

checkConvertInputClasses <- function(signal=NULL,
                                     signalCoord=NULL,
                                     regionSet=NULL,
                                     signalCol = NULL,
                                     GRList = NULL,
                                     .env=.enclosingEnv) {
    
    # default value for .env, I put this inside the function to clarify
    # that parent.frame() is being called inside the function 
    # (it is still called while inside the function when it is a 
    # default parameter but not when it is explicitly assigned to .env
    # in the function call)
    # I can define .enclosingEnv inside the function because of the 
    # lazy evaluation of .env
    .enclosingEnv = parent.frame(n=1)
    
    if (!is.null(signal)) {
        # preferred as matrix, data.frame works
        if (!(is(signal, "matrix") || is(signal, "data.frame"))) {
            stop("`signal` should be a matrix. Check object class.")
        }
    }
    if (!is.null(signalCoord)) {
        # if (is(signalCoord, "GRanges")) {
            # coordinateDT <- grToDt(signalCoord)
            # assign("coordinateDT", coordinateDT, envir=.env)
        # } 
        # signalCoord should be a GRanges
        if (is(signalCoord, "data.frame")) {
            signalCoord <- dtToGr(signalCoord)
            assign("signalCoord", signalCoord, envir=.env)
        } else if (!is(signalCoord, "GRanges")) {
            stop("signalCoord should be a GRanges object.")
        }
    }
    if (!is.null(regionSet)) {
        if (!is(regionSet, "GRanges")) {
            stop("regionSet should be a GRanges object. Check object class.")
        }
    }
    if (!is.null(signalCol)) {
        if (!is(signalCol, "character")) {
            stop("signalCol should be a character object (eg 'PC1').")
        }
    }
    if (!is.null(GRList)) {
        # should be GRangesList
        if (is(GRList, "GRanges")) {
            GRList <- GRangesList(GRList)
            assign("GRList", GRList, envir=.env)
        } else if (!is(GRList, "GRangesList")) {
            stop("GRList should be a GRangesList object.")
        }
    }
    
}

# convenience function to help in maintaining scoring methods
# @param scoringContext character. "singleBase", "multiBase", or "both"
# also can be "metaRegionProfile"
# "metaRegionProfile is a separate category in case some methods do not work
# for binned region sets
getScoringMethods <- (scoringContext = "both") {
    
    if (scoringContext == "singleBase") {
        sMethods <- c("default", "regionMean", "simpleMean", 
                      "regionMedian")
    } else if (scoringContext == "multiBase") {
        sMethods <- c("default", "simpleMean", 
                      "proportionWeightedMean")
    } else if (scoringContext == "both") {
        sMethods <- c("default", "regionMean", "simpleMean", 
                     "regionMedian",
                     # "meanDiff", "rankSum", 
                     "proportionWeightedMean")
    } else if (scoringContext == "metaRegionProfile") {
        sMethods <- c("default", "regionMean", "simpleMean", 
                     "proportionWeightedMean")
    } else {
        stop("Invalid scoringContext specified")
    }
    
    return(sMethods)
}


# Converts a list of data.tables into GRanges.
# @param dtList A list of data.tables, 
# Each should have "chr", "start", "methylCount", and "coverage" columns.
# Error results if missing "chr", "start" but if methylCount and coverage are
# missing, it will still work, just not have that info in the output.
# @return a list of GRanges objects, strand has been set to "*", 
# "start" and "end" have both been set to "start" of the DT.
# methylCount and coverage info is not included in GRanges object.
BSdtToGRanges <- function(dtList) {
    
    # if input are data.frame but not data.table then convert to data.table
    convertInd = sapply(X = dtList, FUN = function(x) is(x, "data.frame") & !is(x, "data.table"))
    dtList[convertInd] = lapply(X = dtList[convertInd], FUN = as.data.table)  
    
    gList <- list();
    for (i in seq_along(dtList)) {
        # dt <- dtList[[i]];
        if ("end" %in% colnames(dtList[[i]])) {
            #message("end is present")  # DEBUG
            setkey(dtList[[i]], chr, start, end)
            # convert the data into granges object
            gList[[i]] <- GRanges(seqnames = dtList[[i]]$chr, 
                                  ranges = IRanges(start = dtList[[i]]$start, 
                                                   end = dtList[[i]]$end), 
                                  strand = rep("*", nrow(dtList[[i]])))
        } else {
            setkey(dtList[[i]], chr, start)
            # convert the data into granges object
            gList[[i]] <- GRanges(seqnames = dtList[[i]]$chr, 
                                  ranges = IRanges(start = dtList[[i]]$start, 
                                                   end = dtList[[i]]$start), 
                                  strand = rep("*", nrow(dtList[[i]])))
            
            # I used to use end = start + 1, but this targets CG instead of just 
            # a C, and it's causing edge-effects problems when I assign Cs to 
            # tiled windows using (within). Aug 2014 I'm changing to start/end at 
            # the same coordinate.
        }
    }
    return(gList);
}




# Convert a GRanges object into a data.table
#
# Also can convert GPos objects to a data.table.
# 
# @param GR A GRanges object
# @param includeStrand "logical" object, whether to include strand from GR in output DT
# @return A data.table object with columns:
# "chr", "start", and "end" (possibly strand)
grToDt <- function(GR, includeStrand = FALSE) {
    DF <- as.data.frame(elementMetadata(GR))
    if ( ncol(DF) > 0) {
        if (includeStrand) {
            DT <- data.table(chr = as.vector(seqnames(GR)), 
                             start = start(GR), 
                             end = end(GR), 
                             strand = as.vector(strand(GR), mode = "character"), 
                             DF)    
        } else{
            DT <- data.table(chr = as.vector(seqnames(GR)), 
                             start = start(GR), 
                             end = end(GR), 
                             DF)    
        }
    } else {
        if (includeStrand) {
            DT <- data.table(chr = as.vector(seqnames(GR)), 
                             start = start(GR), 
                             end = end(GR), 
                             strand = as.vector(strand(GR), mode = "character")) 
        } else{
            DT <- data.table(chr = as.vector(seqnames(GR)), 
                             start = start(GR), 
                             end = end(GR))    
        }
    }
    return(DT)
}





# Convert a data.table to GRanges object.
# 
# @param DT a data.table with at least "chr" and "start" columns
# 
# @return gr A genomic ranges object derived from DT
dtToGr <- function(DT, chr = "chr", start = "start", 
                   end = NA, strand = NA, name = NA, 
                   splitFactor = NA, metaCols = NA) {
    
    if (is.na(splitFactor)) {
        return(dtToGrInternal(DT, chr, start, end, strand, name, metaCols));
    }
    if ( length(splitFactor) == 1 ) { 
        if (splitFactor %in% colnames(DT)) {
            splitFactor <- DT[, get(splitFactor)];
        }
    }
    lapply(split(seq_len(nrow(DT)), splitFactor), 
           function(x) { 
               dtToGrInternal(DT[x, ], chr, start, end, strand, name, metaCols)
           }
    )
}

dtToGR <- dtToGr;


# Internal part of a utility to convert data.tables into GRanges objects
# 
# @param DT A data.table with at least "chr" and "start" columns
# @return gr A genomic ranges object derived from DT
dtToGrInternal <- function(DT, chr, start, 
                           end = NA, strand = NA, name = NA, metaCols = NA) {
    
    if (is.na(end)) {
        if ("end" %in% colnames(DT)) {
            end <- "end"
        } else {
            end <- start;
        }
    }
    if (is.na(strand)) {
        if ("strand" %in% colnames(DT)) { # checking if strand info is in DT
            strand <- "strand"
        }
    }
    if (is.na(strand)) {
        gr <- GRanges(seqnames = DT[[`chr`]], 
                      ranges = IRanges(start = DT[[`start`]], end = DT[[`end`]]), 
                      strand = "*")
    } else {
        # GRanges can only handle '*' for no strand, so replace any non-accepted
        # characters with '*'
        DT[, strand := as.character(strand)]
        DT[strand == "1", strand := "+"]
        DT[strand == "-1", strand := "-"]
        DT[[`strand`]] <- gsub("[^+-]", "*", DT[[`strand`]])
        # chr, start, end and strand should be strings with the 
        # name of the corresponding columns
        gr <- GRanges(seqnames = DT[[`chr`]], 
                      ranges = IRanges(start = DT[[`start`]], end = DT[[`end`]]), 
                      strand = DT[[`strand`]])
    }
    if (! is.na(name)) {
        names(gr) <- DT[[`name`]];
    } else {
        names(gr) <- seq_along(gr);
    }
    if (! is.na(metaCols)) {
        for(x in metaCols) {
            elementMetadata(gr)[[`x`]] <- DT[[`x`]]
        }
    }
    gr;
}


# cleanws takes multi-line, code formatted strings and just formats them
# as simple strings
# @param string string to clean
# @return A string with all consecutive whitespace characters, including
# tabs and newlines, merged into a single space.
cleanws <- function(string) {
    return(gsub('\\s+'," ", string))
}


# Function to run lapply or mclapply, depending on the option set in
# getOption("mc.cores"), which can be set with setLapplyAlias().
#
# @param ... Arguments passed lapply() or mclapply()
# @param mc.preschedule Argument passed to mclapply
# @return Result from lapply or parallel::mclapply
lapplyAlias <- function(..., mc.preschedule = TRUE) {
    if (is.null(getOption("mc.cores"))) { setLapplyAlias(1) }
    if (getOption("mc.cores") > 1) {
        return(parallel::mclapply(..., mc.preschedule = mc.preschedule))
    } else {
        return(lapply(...))
    }
}



# To make parallel processing a possibility but not required, 
# I use an lapply alias which can point at either the base lapply
# (for no multicore), or it can point to mclapply, 
# and set the options for the number of cores (what mclapply uses).
# With no argument given, returns intead the number of cpus currently selected.
#
# @param cores Number of cpus
# @return None
setLapplyAlias <- function(cores = 0) {
    if (cores < 1) {
        return(getOption("mc.cores"))
    }
    if (cores > 1) { # use multicore?
        if (requireNamespace("parallel", quietly = TRUE)) {
            options(mc.cores = cores)
        } else {
            warning(cleanws("You don't have package parallel installed. 
                            Setting cores to 1."))
            options(mc.cores = 1) # reset cores option.
        }
    } else {
        options(mc.cores = 1) # reset cores option.
    }
}


# helper function
# given a vector of columns, and the equally-sized vector of functions
# to apply to those columns, constructs a j-expression for use in
# a data.table 
# (functions applied to columns in corresponding spot in "cols" string).
# One function may be given to be applied to multiple columns.
# use it in a DT[, eval(parse(text = buildJ(cols, funcs)))]
# @param cols A string/vector of strings containing columns 
# on which to use functions.
# @param funcs Functions to use on columns.
# @return A jcommand string. After performing function on column, column 
# is reassigned the same name.
buildJ <- function(cols, funcs, newColNames=NULL) {
    if (is.null(newColNames)) {
        # previously the only option (changed 10/16/17)
        r <- paste("list(", paste(paste0(cols, "=", funcs, "(", cols, ")"), collapse = ","), ")")
    } else {
        r <- paste("list(", paste(paste0(newColNames, "=", funcs, "(", cols, ")"), collapse = ","), ")")
    }
    return(r);
}