# Unit tests

library(PCRSA)
library(data.table)
library(GenomicRanges)

context("Unit tests for PCRSA.")


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




# package data to use for tests
data("brcaCoord1")
data("brcaLoadings1")
data("esr1_chr1")
data("gata3_chr1")

# running the tests


# aggregateLoadings(loadingMat = loadingMat, 
#                   coordinateDT = coordinates, 
#                   regionSet = regionSet)

# cpgOLMetrics
dataDT = cbind(coordinateDT, as.data.frame(loadingMat))
# cpgOLMetrics(dataDT=dataDT, regionGR = regionSet1, metrics = c("mean", "sd"), alsoNonOLMet = TRUE)


# making test data for Wilcoxon rank sum test
chr3 = seq(from=100, to = 700, by = 100)
coordinateDTW = data.table(chr=rep("chr3", length(chr3)),
                           start = chr3)
regionSetW = data.table(chr = c("chr3", "chr3"), 
                        start = c(50, 250),
                        end = c(150, 450))
regionSetW = PCRSA:::dtToGr(regionSetW)
loadingMatW = matrix(data = c(-2:4, seq(from=8, to=2, by=-1)), ncol = 2)
colnames(loadingMatW) <- c("PC2", "PC3")
loadingMatW[7, "PC3"] = 10
dataDTW = cbind(coordinateDTW, as.data.frame(loadingMatW))

test_that("aggregateLoadings and scoring metrics", {
    

    # test wilcoxon rank sum scoring metric
    rsWResults = PCRSA:::rsWilcox(dataDT = dataDTW, regionGR = regionSetW)
    PC2W = wilcox.test(x = c(-2, 0, 1), y=c(-1, 2:4))$p.value
    PC3W = wilcox.test(x = c(8, 6, 5), y=c(7, 4, 3, 10))$p.value
    expect_equal(c(PC2W, PC3W, 3, 2, 2, mean(width(regionSetW))), 
                 c(rsWResults$PC2, rsWResults$PC3, rsWResults$cytosine_coverage, 
                   rsWResults$region_coverage, 
                   rsWResults$total_region_number, 
                   rsWResults$mean_region_size))
    
    # same test for Wilcoxon but with aggregateLoadings (absolute value
    # of loadings will be taken)
    rsWResults = PCRSA:::aggregateLoadings(loadingMat = loadingMatW, 
                      mCoord = coordinateDTW, 
                      regionSet = regionSetW, 
                      PCsToAnnotate = c("PC2", "PC3"), 
                      metric = "rankSum")    
    PC2W = wilcox.test(x = c(2, 0, 1), y=c(1, 2:4), alternative = "greater")$p.value
    PC3W = wilcox.test(x = c(8, 6, 5), y=c(7, 4, 3, 10), alternative = "greater")$p.value
    expect_equal(c(PC2W, PC3W, 3, 2, 2, mean(width(regionSetW))), 
                 c(rsWResults$PC2, rsWResults$PC3, rsWResults$cytosine_coverage, 
                   rsWResults$region_coverage, 
                   rsWResults$total_region_number, 
                   rsWResults$mean_region_size))
    
    
    
    
    # test raw scoring method
    # aggregateLoadings(loadingMat = , 
    #                   mCoord = , 
    #                   regionSet = , 
    #                   PCsToAnnotate = , 
    #                   metric = , 
    #                   pcLoadAv = )
    # expect_equal()
    # 
    
    # test mean difference scoring method
})

# test_that("", {
#     
#     expect_equal()
#     
# })
# 
# test_that("", {
#     
#     expect_equal()
#     
# })
# 
# test_that("", {
#     
#     expect_equal()
#     
# })

######################### Visualization functions #############################


