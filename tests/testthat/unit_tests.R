
# generate synthetic data
# loadingMat = 
coordinates = data.table(chr=rep("chr1", nrow(loadingMat)), start = 1:nrow(loadingMat))
regionSet = data.table(chr = c("chr1", "chr1", "chr1", "chr1"), 
                       start = c(1, 500, 800, 1000),
                       end = c(400, 700, 950, 1500))
regionSet = MIRA:::dtToGr(regionSet)

aggregateLoadings(loadingMat = loadingMat, 
                  coordinateDT = coordinates, 
                  regionSet = regionSet)



