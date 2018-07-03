
library(data.table)
library(GenomicRanges)

# generate synthetic data
# coordinateDT: positions of cytosines
chr1 = seq(from=1, to = 4000, by = 200)
chr2 = seq(from=4001, to = 6400, by = 200)
start = c(chr1, chr2)
coordinateDT = data.table(chr=c(rep("chr1", length(chr1)),
                             rep("chr2", length(chr2))),
                         start = start)
# loading values
loadingMat = matrix(data = rep(1, (length(start) * 2)), ncol = 2)
colnames(loadingMat) <- c("PC1", "PC3")

dataDT = data.table(coordinateDT, loadingMat)
# arbitrarily make some != 1
dataDT$PC1[3] = 0
dataDT$PC3[5] = 2

# region sets that will be tested
regionSet1 = data.table(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"), 
                       start = c(1, 500, 3100, 5100, 6000),
                       end = c(400, 700, 3400, 5150, 6450))
regionSet1 = MIRA:::dtToGr(regionSet1)
regionSet2 = data.table(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"), 
                                     start = c(1000, 1500, 3700, 4000, 4300),
                                     end = c(1400, 1700, 3800, 4100, 4700))
regionSet2 = MIRA:::dtToGr(regionSet2)

# running the tests


# aggregateLoadings(loadingMat = loadingMat, 
#                   coordinateDT = coordinates, 
#                   regionSet = regionSet)

# cpgOLMetrics
dataDT = cbind(coordinateDT, as.data.frame(loadingMat))
cpgOLMetrics(dataDT=dataDT, regionGR = regionSet1, metrics = c("mean", "sd"), alsoNonOLMet = TRUE)

