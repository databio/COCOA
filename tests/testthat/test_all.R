# Unit tests

library(COCOA)
library(data.table)
library(GenomicRanges)

context("Unit tests for COCOA.")


# generate synthetic data
# coordinateDT: positions of cytosines
chr1 <- seq(from=1, to = 4000, by = 200)
chr2 <- seq(from=4001, to = 6400, by = 200)
start <- c(chr1, chr2)
coordinateDT <- data.table(chr=c(rep("chr1", length(chr1)),
                             rep("chr2", length(chr2))),
                         start = start)
# loading values
loadingMat <- matrix(data = rep(1, (length(start) * 2)), ncol = 2)
colnames(loadingMat) <- c("PC1", "PC3")

dataDT <- data.table(coordinateDT, loadingMat)
# arbitrarily make some != 1
dataDT$PC1[3] <- 0
dataDT$PC3[5] <- 2

# region sets that will be tested
regionSet1 <- data.table(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"), 
                       start = c(1, 500, 3100, 5100, 6000),
                       end = c(400, 700, 3400, 5150, 6450))
regionSet1 <- MIRA:::dtToGr(regionSet1)
regionSet2 <- data.table(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"), 
                                     start = c(1000, 1500, 3700, 4000, 4300),
                                     end = c(1400, 1700, 3800, 4100, 4700))
regionSet2 <- MIRA:::dtToGr(regionSet2)




# package data to use for tests
# don't load unnecessary data in order to save time
data("brcaCoord1")
data("brcaLoadings1")
data("esr1_chr1")


# running the tests


# aggregateLoadings(loadingMat = loadingMat, 
#                   coordinateDT = coordinates, 
#                   regionSet = regionSet)

# cpgOLMetrics
dataDT <- cbind(coordinateDT, as.data.frame(loadingMat))


# making test data for Wilcoxon rank sum test
chr3 <- seq(from=100, to = 700, by = 100)
coordinateDTW <- data.table(chr=rep("chr3", length(chr3)),
                           start = chr3, end = chr3)
regionSetW <- data.table(chr = c("chr3", "chr3"), 
                        start = c(50, 250),
                        end = c(150, 450))
regionSetW <- COCOA:::dtToGr(regionSetW)
loadingMatW <- matrix(data = c(-2:4, seq(from=8, to=2, by=-1)), ncol = 2)
colnames(loadingMatW) <- c("PC2", "PC3")
loadingMatW[7, "PC3"] <- 10
dataDTW <- cbind(coordinateDTW, as.data.frame(loadingMatW))

test_that("aggregateLoadings, scoring metrics, and pcRegionSetEnrichment", {
    

    # test wilcoxon rank sum scoring metric
    rsWResults <- COCOA:::rsWilcox(dataDT = dataDTW, regionGR = regionSetW)
    PC2W <- wilcox.test(x = c(-2, 0, 1), y=c(-1, 2:4))$p.value
    PC3W <- wilcox.test(x = c(8, 6, 5), y=c(7, 4, 3, 10))$p.value
    expect_equal(c(PC2W, PC3W, 3, 2, 2, mean(width(regionSetW))), 
                 c(rsWResults$PC2, rsWResults$PC3, rsWResults$cytosine_coverage, 
                   rsWResults$region_coverage, 
                   rsWResults$total_region_number, 
                   rsWResults$mean_region_size))
    
    # same test for Wilcoxon but with aggregateLoadings (absolute value
    # of loadings will be taken), "greater" alternate hypothesis is used in
    # aggregateLoadings
    rsWResults <- COCOA:::aggregateLoadings(loadingMat = loadingMatW, 
                      mCoord = coordinateDTW, 
                      regionSet = regionSetW, 
                      PCsToAnnotate = c("PC2", "PC3"), 
                      metric = "rankSum")    
    PC2W <- wilcox.test(x = c(2, 0, 1), y=c(1, 2:4), alternative = "greater")$p.value
    PC3W <- wilcox.test(x = c(8, 6, 5), y=c(7, 4, 3, 10), alternative = "greater")$p.value
    expect_equal(c(PC2W, PC3W, 3, 2, 2, mean(width(regionSetW))), 
                 c(rsWResults$PC2, rsWResults$PC3, rsWResults$cytosine_coverage, 
                   rsWResults$region_coverage, 
                   rsWResults$total_region_number, 
                   rsWResults$mean_region_size))
    
    
    
    
    # test rsMean scoring method, average first within region, then
    # between regions
    rsMeanRes <- COCOA:::aggregateLoadings(loadingMat = loadingMatW, 
                              mCoord = coordinateDTW, 
                              regionSet = regionSetW, 
                              PCsToAnnotate = c("PC2", "PC3"), 
                              metric = "rsMean")
    PC2R <- mean(c(2, mean(c(0, 1))))
    PC3R <- mean(c(8, mean(c(6, 5))))
    expect_equal(c(PC2R, PC3R, 3, 2, 2, mean(width(regionSetW))), 
                 c(rsMeanRes$PC2, rsMeanRes$PC3, rsMeanRes$cytosine_coverage, 
                   rsMeanRes$region_coverage, 
                   rsMeanRes$total_region_number, 
                   rsMeanRes$mean_region_size))
    
    ########### test cpgMean scoring method ################
    # this is a mean of CpG loading values just like "raw" but instead
    # of averaging within regions then averaging regions together, this 
    # method does the simple average of all CpGs within the region set
    cpgMeanRes <- COCOA:::aggregateLoadings(loadingMat = loadingMatW, 
                                       mCoord = coordinateDTW, 
                                       regionSet = regionSetW, 
                                       PCsToAnnotate = c("PC2", "PC3"), 
                                       metric = "cpgMean")
    PC2RC <- mean(c(2, 0, 1))
    PC3RC <- mean(c(8, 6, 5))
    expect_equal(c(PC2RC, PC3RC, 3, 2, 2, mean(width(regionSetW))), 
                 c(cpgMeanRes$PC2, cpgMeanRes$PC3, cpgMeanRes$cytosine_coverage, 
                   cpgMeanRes$region_coverage, 
                   cpgMeanRes$total_region_number, 
                   cpgMeanRes$mean_region_size))
    

    # test mean difference scoring method
    PC2Num <- mean(c(2, 0, 1)) - mean(c(1, 2, 3, 4))
    PC3Num <- mean(c(8, 6, 5)) - mean(c(7, 4, 3, 10))
    # see doi: 10.1186/s13040-015-0059-z
    # pooled variance times normalization factor for size of each set
    PC2Denom <- sqrt((sd(c(2, 0, 1))^2 + sd(c(1, 2, 3, 4))^2) / 2) * sqrt(1/3 - 1/4)
    PC3Denom <- sqrt((sd(c(8, 6, 5))^2 + sd(c(7, 4, 3, 10))^2) / 2) * sqrt(1/3 - 1/4)
    
    PC2MD <- PC2Num / PC2Denom
    PC3MD <- PC3Num / PC3Denom
    mdRes <- COCOA:::aggregateLoadings(loadingMat = loadingMatW, 
                                       mCoord = coordinateDTW, 
                                       regionSet = regionSetW, 
                                       PCsToAnnotate = c("PC2", "PC3"), 
                                       metric = "meanDiff")
    expect_equal(c(PC2MD, PC3MD, 3, 2, 2, mean(width(regionSetW))), 
                 c(mdRes$PC2, mdRes$PC3, mdRes$cytosine_coverage, 
                   mdRes$region_coverage, 
                   mdRes$total_region_number, 
                   mdRes$mean_region_size))
    
    ################## test pcRegionSetEnrichment with meanDiff test data
    
    coordinateDTW2 <- copy(coordinateDTW)
    coordinateDTW2$end <- coordinateDTW$end  + 1 # so there will be an actual range
    # make sure it works if there is metadata in GRanges object
    extraCol <- rep(0, nrow(coordinateDTW2))
    mCoordW2 <- COCOA:::dtToGr(coordinateDTW2)
    mcols(mCoordW2) <- data.frame(extraCol)
    twoResults <- pcRegionSetEnrichment(loadingMat = loadingMatW, 
                          mCoord = mCoordW2, 
                          GRList = GRangesList(regionSetW, regionSetW), 
                          PCsToAnnotate = c("PC2", "PC3"), 
                          scoringMetric = "meanDiff")
    expect_equal(rep(c(PC2MD, PC3MD, 3, 2, 2, mean(width(regionSetW))), each=2), 
                 c(twoResults$PC2, twoResults$PC3, twoResults$cytosine_coverage, 
                   twoResults$region_coverage, 
                   twoResults$total_region_number, 
                   twoResults$mean_region_size))
    
    ########## testing cpgOLMetrics, used for meanDiff scoring method #######
    olMetrics <- COCOA:::cpgOLMetrics(dataDT=dataDTW, 
                 regionGR = regionSetW, 
                 metrics = c("mean", "sd"), 
                 alsoNonOLMet = TRUE)
    
    PC2Met <- data.table(t(c(mean(c(-2, 0, 1)), sd(c(-2, 0, 1)), 
                   mean(c(-1, 2:4)), sd(c(-1, 2:4)), 3, 2, 2, 151)))
    PC2Met <- cbind("PC2", PC2Met)
    colnames(PC2Met) <- colnames(olMetrics)
    expect_equal(olMetrics[1, ], PC2Met)
    
    PC3Met <- data.table(t(c(mean(c(8, 6, 5)), sd(c(8, 6, 5)), 
                            mean(c(7, 4, 3, 10)), sd(c(7, 4, 3, 10)), 3, 2, 2, 151)))
    PC3Met <- cbind("PC3", PC3Met)
    colnames(PC3Met) <- colnames(olMetrics)
    expect_equal(olMetrics[2, ], PC3Met)
    
})

test_that("aggregateLoadings and scoring metrics", {
    
    
    
    # test wilcoxon rank sum scoring metric
    rsWResults <- COCOA:::rsWilcox(dataDT = dataDTW, regionGR = regionSetW)
    PC2W <- wilcox.test(x = c(-2, 0, 1), y=c(-1, 2:4))$p.value
    PC3W <- wilcox.test(x = c(8, 6, 5), y=c(7, 4, 3, 10))$p.value
    expect_equal(c(PC2W, PC3W, 3, 2, 2, mean(width(regionSetW))), 
                 c(rsWResults$PC2, rsWResults$PC3, rsWResults$cytosine_coverage, 
                   rsWResults$region_coverage, 
                   rsWResults$total_region_number, 
                   rsWResults$mean_region_size))
    
    # same test for Wilcoxon but with aggregateLoadings (absolute value
    # of loadings will be taken), "greater" alternate hypothesis is used in
    # aggregateLoadings
    rsWResults <- COCOA:::aggregateLoadings(loadingMat = loadingMatW, 
                                           mCoord = coordinateDTW, 
                                           regionSet = regionSetW, 
                                           PCsToAnnotate = c("PC2", "PC3"), 
                                           metric = "rankSum")    
    PC2W <- wilcox.test(x = c(2, 0, 1), y=c(1, 2:4), alternative = "greater")$p.value
    PC3W <- wilcox.test(x = c(8, 6, 5), y=c(7, 4, 3, 10), alternative = "greater")$p.value
    expect_equal(c(PC2W, PC3W, 3, 2, 2, mean(width(regionSetW))), 
                 c(rsWResults$PC2, rsWResults$PC3, rsWResults$cytosine_coverage, 
                   rsWResults$region_coverage, 
                   rsWResults$total_region_number, 
                   rsWResults$mean_region_size))

    
    
    })

test_that("averageByRegion", {
    loadingMatABR <- loadingMat
    loadingMatABR[1, "PC1"] <- 3
    loadingMatABR[32, "PC1"] <- 2
    abr <- COCOA:::averageByRegion(loadingMat = loadingMatABR, mCoord = COCOA:::dtToGr(coordinateDT), 
                            regionSet = regionSet1, PCsToAnnotate = "PC1", 
                            returnQuantile = FALSE)
    expect_equal(abr$PC1, c(2, 1, 1, 1.5))
    
    # test quantile
    abrq <- COCOA:::averageByRegion(loadingMat = loadingMatABR, mCoord = COCOA:::dtToGr(coordinateDT), 
                                   regionSet = regionSet1, PCsToAnnotate = "PC1", 
                                   returnQuantile = TRUE)
    # the mean values, abr$PC1, converted to quantiles
    # converting properly to quantiles?
    expect_equal(abrq$PC1, ecdf(loadingMatABR[, "PC1"])(abr$PC1))
    
})

test_that("averageByRegion", {
    loadingMatABR <- loadingMat
    loadingMatABR[1, "PC1"] <- 3
    loadingMatABR[32, "PC1"] <- 2
    abr <- COCOA:::averageByRegion(loadingMat = loadingMatABR, mCoord = COCOA:::dtToGr(coordinateDT), 
                                   regionSet = regionSet1, PCsToAnnotate = "PC1", 
                                   returnQuantile = FALSE)
    expect_equal(abr$PC1, c(2, 1, 1, 1.5))
    
    # test quantile
    abrq <- COCOA:::averageByRegion(loadingMat = loadingMatABR, mCoord = COCOA:::dtToGr(coordinateDT), 
                                    regionSet = regionSet1, PCsToAnnotate = "PC1", 
                                    returnQuantile = TRUE)
    # the mean values, abr$PC1, converted to quantiles
    # converting properly to quantiles?
    expect_equal(abrq$PC1, ecdf(loadingMatABR[, "PC1"])(abr$PC1))
    
})

test_that("pcEnrichmentProfile", {
    loadingMatP <- matrix(data = c(1:8, 10, 10, rep(1, 20)), nrow = 30)
    colnames(loadingMatP) <- "PC1"
    # two CpGs per bin of regionSetP (with 5 bins)
    .start <- c(seq(from=50, to=950, by = 100),
              seq(from=2050, to=2950, by = 100),
              seq(from=4050, to=4950, by = 100))
    coordinateDTP <- data.frame(chr = rep("chr3", nrow(loadingMatP)), 
                               start = .start,
                               end = .start + 1,
                               extraCol = rep(1, length(.start)))
    regionSetP <- data.table(chr = rep("chr3", 3), 
                            start = c(1, 2001, 4001),
                            end = c(1010, 3010, 5010))
    regionSetP <- COCOA:::dtToGr(regionSetP)
     
    
    binnedP <- pcEnrichmentProfile(loadingMat = loadingMatP, mCoord = COCOA:::dtToGr(coordinateDTP), 
                        GRList = regionSetP, PCsToAnnotate = "PC1", 
                        binNum = 5)    
    meanPerBin <- (c(seq(from=1.5, to=7.5, by=2), 10) + rep(1, 5) + rep(1, 5)) / 3
    # BSBinAggregate averages the profile around the center
    symmetricalBin <- (meanPerBin + rev(meanPerBin)) / 2
    expect_equal(binnedP[[1]]$PC1, symmetricalBin)
    
})


################### test function inputs
# make sure various input formats work

#### take out these high level tests in order to save time
# test_that("Input formats", {
#     
#     normalOut <- pcRegionSetEnrichment(loadingMat = brcaLoadings1, 
#                           mCoord = brcaCoord1, 
#                           GRList = GRangesList(esr1_chr1), 
#                           PCsToAnnotate = "PC1", 
#                           scoringMetric = "rsMean")
#     alterOut <- pcRegionSetEnrichment(loadingMat = data.frame(brcaLoadings1), 
#                           mCoord = brcaCoord1, 
#                           GRList = esr1_chr1, 
#                           PCsToAnnotate = "PC1", 
#                           scoringMetric = "rsMean")
#     expect_equal(normalOut, alterOut)
#     
#     alterOut <- pcRegionSetEnrichment(loadingMat = data.frame(brcaLoadings1), 
#                                       mCoord = dtToGr(brcaCoord1), 
#                                       GRList = esr1_chr1, 
#                                       PCsToAnnotate = "PC1", 
#                                       scoringMetric = "rsMean")
#     expect_equal(normalOut, alterOut)
#     
#     #########
#     
#     normalOut <- pcEnrichmentProfile(loadingMat = brcaLoadings1, 
#                                        mCoord = brcaCoord1, 
#                                        GRList = GRangesList(esr1_chr1), 
#                                        PCsToAnnotate = "PC1", 
#                                        binNum = 5)
#     alterOut <- pcEnrichmentProfile(loadingMat = brcaLoadings1, 
#                                     mCoord = COCOA:::dtToGr(brcaCoord1), 
#                                     GRList = esr1_chr1, 
#                                     PCsToAnnotate = "PC1", 
#                                     binNum = 5)
#     expect_equal(normalOut, alterOut)
# 
#     
#     
# })
