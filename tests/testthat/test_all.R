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

# fake region data (for testing ATAC-seq)
regionCoordDT <- data.table(chr=c(rep("chr1",4), rep("chr2", 2)),
                            start=c(300, 650, 3000, 10000, 5075, 6150),
                            end=c(600, 1000, 4000, 11000, 5125, 6300))
regionDataDT <- data.table(PC1=c(1:6), PC2=-1:4)


# region sets that will be tested
regionSet1 <- data.table(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"), 
                       start = c(1, 500, 3100, 5100, 6000),
                       end = c(400, 700, 3400, 5150, 6450))
regionSet1 <- MIRA:::dtToGr(regionSet1)
regionSet2 <- data.table(chr = c("chr1", "chr1", "chr1", "chr2", "chr2"), 
                                     start = c(1000, 1500, 3700, 4000, 4300),
                                     end = c(1400, 1700, 3800, 4100, 4700))
regionSet2 <- MIRA:::dtToGr(regionSet2)




# don't use built in package data for tests because it will take longer

# running the tests


# aggregateSignal(signal = loadingMat, 
#                   coordinateDT = coordinates, 
#                   regionSet = regionSet)

# signalOLMetrics
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

test_that("aggregateSignal, scoring metrics, and runCOCOA", {
    

    # test wilcoxon rank sum scoring metric
    rsWResults <- COCOA:::rsWilcox(dataDT = dataDTW, regionSet = regionSetW)
    PC2W <- wilcox.test(x = c(-2, 0, 1), y=c(-1, 2:4))$p.value
    PC3W <- wilcox.test(x = c(8, 6, 5), y=c(7, 4, 3, 10))$p.value
    expect_equal(c(PC2W, PC3W, 3, 2, 2, mean(width(regionSetW))), 
                 c(rsWResults$PC2, rsWResults$PC3, rsWResults$signalCoverage, 
                   rsWResults$regionSetCoverage, 
                   rsWResults$totalRegionNumber, 
                   rsWResults$meanRegionSize))
    
    # same test for Wilcoxon but with aggregateSignal (absolute value
    # of loadings will be taken), "greater" alternate hypothesis is used in
    # aggregateSignal
    rsWResults <- COCOA:::aggregateSignal(signal = loadingMatW, 
                      signalCoord = coordinateDTW, 
                      regionSet = regionSetW, 
                      signalCol = c("PC2", "PC3"), 
                      scoringMetric = "rankSum")    
    PC2W <- wilcox.test(x = c(2, 0, 1), y=c(1, 2:4), alternative = "greater")$p.value
    PC3W <- wilcox.test(x = c(8, 6, 5), y=c(7, 4, 3, 10), alternative = "greater")$p.value
    expect_equal(c(PC2W, PC3W, 3, 2, 2, mean(width(regionSetW))), 
                 c(rsWResults$PC2, rsWResults$PC3, rsWResults$signalCoverage, 
                   rsWResults$regionSetCoverage, 
                   rsWResults$totalRegionNumber, 
                   rsWResults$meanRegionSize))
    
    
    
    
    # test regionMean scoring method, average first within region, then
    # between regions
    regionMeanRes <- COCOA:::aggregateSignal(signal = loadingMatW, 
                              signalCoord = coordinateDTW, 
                              regionSet = regionSetW, 
                              signalCol = c("PC2", "PC3"), 
                              scoringMetric = "regionMean")
    PC2R <- mean(c(2, mean(c(0, 1))))
    PC3R <- mean(c(8, mean(c(6, 5))))
    expect_equal(c(PC2R, PC3R, 3, 2, 2, mean(width(regionSetW))), 
                 c(regionMeanRes$PC2, regionMeanRes$PC3, regionMeanRes$signalCoverage, 
                   regionMeanRes$regionSetCoverage, 
                   regionMeanRes$totalRegionNumber, 
                   regionMeanRes$meanRegionSize))
    
    ########### test simpleMean scoring method ################
    # this is a mean of CpG loading values just like "raw" but instead
    # of averaging within regions then averaging regions together, this 
    # method does the simple average of all CpGs within the region set
    simpleMeanRes <- COCOA:::aggregateSignal(signal = loadingMatW, 
                                       signalCoord = coordinateDTW, 
                                       regionSet = regionSetW, 
                                       signalCol = c("PC2", "PC3"), 
                                       scoringMetric = "simpleMean")
    PC2RC <- mean(c(2, 0, 1))
    PC3RC <- mean(c(8, 6, 5))
    expect_equal(c(PC2RC, PC3RC, 3, 2, 2, mean(width(regionSetW))), 
                 c(simpleMeanRes$PC2, simpleMeanRes$PC3, simpleMeanRes$signalCoverage, 
                   simpleMeanRes$regionSetCoverage, 
                   simpleMeanRes$totalRegionNumber, 
                   simpleMeanRes$meanRegionSize))
    

    # test mean difference scoring method
    PC2Num <- mean(c(2, 0, 1)) - mean(c(1, 2, 3, 4))
    PC3Num <- mean(c(8, 6, 5)) - mean(c(7, 4, 3, 10))
    # see doi: 10.1186/s13040-015-0059-z
    # pooled variance times normalization factor for size of each set
    PC2Denom <- sqrt((sd(c(2, 0, 1))^2 + sd(c(1, 2, 3, 4))^2) / 2) * sqrt(1/3 - 1/4)
    PC3Denom <- sqrt((sd(c(8, 6, 5))^2 + sd(c(7, 4, 3, 10))^2) / 2) * sqrt(1/3 - 1/4)
    
    PC2MD <- PC2Num / PC2Denom
    PC3MD <- PC3Num / PC3Denom
    mdRes <- COCOA:::aggregateSignal(signal = loadingMatW, 
                                       signalCoord = coordinateDTW, 
                                       regionSet = regionSetW, 
                                       signalCol = c("PC2", "PC3"), 
                                       scoringMetric = "meanDiff")
    expect_equal(c(PC2MD, PC3MD, 3, 2, 2, mean(width(regionSetW))), 
                 c(mdRes$PC2, mdRes$PC3, mdRes$signalCoverage, 
                   mdRes$regionSetCoverage, 
                   mdRes$totalRegionNumber, 
                   mdRes$meanRegionSize))
    
    ################## test runCOCOA with meanDiff test data
    
    coordinateDTW2 <- copy(coordinateDTW)
    coordinateDTW2$end <- coordinateDTW$end
    # make sure it works if there is metadata in GRanges object
    extraCol <- rep(0, nrow(coordinateDTW2))
    signalCoordW2 <- COCOA:::dtToGr(coordinateDTW2)
    mcols(signalCoordW2) <- data.frame(extraCol)
    twoResults <- runCOCOA(signal = loadingMatW, 
                          signalCoord = signalCoordW2, 
                          GRList = GRangesList(regionSetW, regionSetW), 
                          signalCol = c("PC2", "PC3"), 
                          scoringMetric = "meanDiff")
    expect_equal(rep(c(PC2MD, PC3MD, 3, 2, 2, mean(width(regionSetW))), each=2), 
                 c(twoResults$PC2, twoResults$PC3, twoResults$signalCoverage, 
                   twoResults$regionSetCoverage, 
                   twoResults$totalRegionNumber, 
                   twoResults$meanRegionSize))
    
    ########## testing signalOLMetrics, used for meanDiff scoring method #######
    olMetrics <- COCOA:::signalOLMetrics(dataDT=dataDTW, 
                 regionSet = regionSetW, 
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

test_that("ATAC-seq scoring methods", {
    
    # test "regionOLWeightedMean"
    weightedAve <- regionOLWeightedMean(signalDT = regionDataDT, signalGR = COCOA:::dtToGr(regionCoordDT), 
                         regionSet = regionSet1, calcCols = c("PC1", "PC2"))
    # proportion overlap is first then PC
    correctAve <- data.frame(PC1=((101/400)*1+(101/201)*1+(51/201)*2+1*3+(26/51)*5+(151/451)*6) / 
        ((101/400)+(101/201)+(51/201)+1+(26/51)+(151/451)),
        PC2=((101/400)*-1+(101/201)*-1+(51/201)*0+1*1+(26/51)*3+(151/451)*4) / 
            ((101/400)+(101/201)+(51/201)+1+(26/51)+(151/451)), 
        signalCoverage=5, regionSetCoverage=5,
        sumProportionOverlap=(101/400)+(101/201)+(51/201)+1+(26/51)+(151/451))
    expect_equal(weightedAve, correctAve)
    
    # test "regionOLMean"
    signalAve <- regionOLMean(signalDT = regionDataDT,
                              signalGR = COCOA:::dtToGr(regionCoordDT), 
                              regionSet = regionSet1,
                              calcCols = c("PC1", "PC2"))
    # proportion overlap is first then PC
    correctAve <- data.frame(PC1=(1+1+2+3+5+6)/6,
                             PC2=(-1+-1+0+1+3+4)/6, 
                             signalCoverage=5,
                             regionSetCoverage=5)
    expect_equal(signalAve, correctAve)
    
    # test "weightedAvePerRegion"
    avePerRegion <- weightedAvePerRegion(signalDT= regionDataDT,
                                    signalCoord=COCOA:::dtToGr(regionCoordDT),
                                    regionSet=regionSet1,
                                    calcCols = c("PC1", "PC2"))
    correctAve <- data.table(PC1 = c(1*1, (101/201*1 + 51/201*2)/ (101/201 + 51/201), 1*3, 1*5, 1*6), 
                             PC2 = c(1*-1, (101/201*-1 + 51/201*0)/ (101/201 + 51/201), 1*1, 1*3, 1*4))
    expect_equal(avePerRegion$PC1, correctAve$PC1)
    expect_equal(avePerRegion$PC2, correctAve$PC2)
    
})


test_that("averagePerRegion", {
    loadingMatABR <- loadingMat
    loadingMatABR[1, "PC1"] <- 3
    loadingMatABR[32, "PC1"] <- 2
    abr <- COCOA:::averagePerRegion(signal = loadingMatABR, signalCoord = COCOA:::dtToGr(coordinateDT), 
                            regionSet = regionSet1, signalCol = "PC1", 
                            returnQuantile = FALSE)
    expect_equal(abr$PC1, c(2, 1, 1, 1.5))
    
    # test quantile
    abrq <- COCOA:::averagePerRegion(signal = loadingMatABR, signalCoord = COCOA:::dtToGr(coordinateDT), 
                                   regionSet = regionSet1, signalCol = "PC1", 
                                   returnQuantile = TRUE)
    # the mean values, abr$PC1, converted to quantiles
    # converting properly to quantiles?
    expect_equal(abrq$PC1, ecdf(loadingMatABR[, "PC1"])(abr$PC1))
    
})

test_that("getMetaRegionProfile", {
    loadingMatP <- matrix(data = c(1:8, 10, 10, rep(1, 20)), nrow = 30)
    colnames(loadingMatP) <- "PC1"
    # two CpGs per bin of regionSetP (with 5 bins)
    .start <- c(seq(from=50, to=950, by = 100),
              seq(from=2050, to=2950, by = 100),
              seq(from=4050, to=4950, by = 100))
    coordinateDTP <- data.frame(chr = rep("chr3", nrow(loadingMatP)), 
                               start = .start,
                               end = .start,
                               extraCol = rep(1, length(.start)))
    regionSetP <- data.table(chr = rep("chr3", 3), 
                            start = c(1, 2001, 4001),
                            end = c(1010, 3010, 5010))
    regionSetP <- COCOA:::dtToGr(regionSetP)
     
    # COCOA:::dtToGr(coordinateDTP)
    binnedP <- lapply(GRangesList(regionSetP), 
                     function(x) getMetaRegionProfile(signal = loadingMatP, 
                                                   signalCoord = coordinateDTP, 
                                                   regionSet = x, 
                                                   signalCol = "PC1", 
                                                   binNum = 5))    
    meanPerBin <- (c(seq(from=1.5, to=7.5, by=2), 10) + rep(1, 5) + rep(1, 5)) / 3
    # BSBinAggregate averages the profile around the center
    symmetricalBin <- (meanPerBin + rev(meanPerBin)) / 2
    expect_equal(binnedP[[1]]$PC1, symmetricalBin)
    
    # making sure that different input formats still work
    #########
    alterOut <- lapply(GRangesList(regionSetP), 
                       function(x) getMetaRegionProfile(signal = data.frame(loadingMatP), 
                                                     signalCoord = COCOA:::dtToGr(coordinateDTP), 
                                                     regionSet = x, 
                                                     signalCol = "PC1", 
                                                     binNum = 5))  
    
    expect_equal(binnedP, alterOut)    
    
})


################### test function inputs
# make sure various input formats work



#### take out these high level tests in order to save time
test_that("Input formats", {

    normalOut <- runCOCOA(signal = loadingMatW, 
                                       signalCoord = coordinateDTW, 
                                       GRList = GRangesList(regionSetW), 
                                       signalCol = "PC2", 
                                       scoringMetric = "regionMean")
    alterOut <- runCOCOA(signal = data.frame(loadingMatW), 
                                       signalCoord = COCOA:::dtToGr(coordinateDTW), 
                                       GRList = regionSetW, 
                                       signalCol = "PC2", 
                                       scoringMetric = "regionMean")
    expect_equal(normalOut, alterOut)

})

