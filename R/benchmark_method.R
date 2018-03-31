
# testing whether method works with a small(er) number of samples

source(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/00-init.R"))

# 
setwd(paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/plots/"))
patientMetadata = brcaMetadata
set.seed(1234)


# DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))
simpleCache("combinedBRCAMethyl_noXY")
brcaMList = combinedBRCAMethyl_noXY

# reading in the metadata, will be used to split data 
# into training and test set with balanced ER and PGR status
#restrict patients included in this analysis
patientMetadata = patientMetadata[patientMetadata$subject_ID %in% 
                                      colnames(brcaMList[["methylProp"]]), ]

###########################################################

# load LOLA database
lolaPath = paste0(Sys.getenv("REGIONS"), "LOLACore/hg38/")
#lolaPath = system.file("extdata", "hg19", package="LOLA")
regionSetDB = loadRegionDB(lolaPath)
loRegionAnno = regionSetDB$regionAnno
# a549Ind = grep("a549", loRegionAnno$cellType, ignore.case = TRUE)
sheff_dnaseInd = grep("sheffield_dnase", loRegionAnno$collection, ignore.case = TRUE)
# mcf7Ind = grep("mcf-7", loRegionAnno$cellType, ignore.case = TRUE)
# k562Ind = grep("k562", loRegionAnno$cellType,  ignore.case = TRUE)
# GRList = GRangesList(regionSetDB$regionGRL[c(a549Ind, mcf7Ind)])
GRList = GRangesList(regionSetDB$regionGRL[-sheff_dnaseInd])
# adding ER Chipseq dataset
erSet = fread(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/inst/extdata/",
                     "GSM2305313_MCF7_E2_peaks_hg38.bed"))
setnames(erSet, c("V1", "V2", "V3"), c("chr", "start", "end"))
GRList = c(GRangesList(MIRA:::dtToGr(erSet)), GRList)
names(GRList) = c("GSM2305313_MCF7_E2_peaks_hg38.bed", loRegionAnno$filename[-sheff_dnaseInd])


##############################################################################
# 657 patients total (that have designated + or - ER and PGR status)
############# ~10% (68 patients) #####################
# keep same proportion of patients with ER_status and PGR_status in training 
# and test sets
set.seed(1234)
dataSplit = createDataPartition(y=factor(paste0(patientMetadata$ER_status,"_", patientMetadata$PGR_status)),
                                p = .1, list=FALSE)
trainingIDs = patientMetadata[dataSplit, subject_ID]
testIDs = patientMetadata[-dataSplit, subject_ID]
trainingMData = brcaMList[["methylProp"]][, 
                                          colnames(brcaMList[["methylProp"]]) %in% trainingIDs] 

# do the PCA
simpleCache("allMPCA_10", {
    prcomp(t(trainingMData), center = TRUE)
})
# plot(allMPCA$x[,c("PC1", "PC3")])
allMPCAWeights10 = as.data.table(allMPCA_10$rotation)
# mIQR = apply(trainingMData, 1, IQR)
# simpleCache("top10MPCA", {
#     prcomp(t(trainingMData[mIQR >= quantile(mIQR, 0.9), ]), center = TRUE)
# })
coordinates = brcaMList[["coordinates"]]
# top10Coord = coordinates[mIQR >= quantile(mIQR, 0.9), ]
# top10PCWeights = as.data.table(top10MPCA$rotation)
# top10TSNE = Rtsne(X = top10MPCA$x[, 1:50], pca = FALSE, max_iter=5000,
#                   perplexity = 30)
# plot(top10TSNE$Y)
# run PC region set enrichment analysis
simpleCache("rsEnrichment_10", {
    rsEnrichment_10 = pcRegionSetEnrichment(loadingMat=allMPCAWeights10, coordinateDT = coordinates, 
                                         GRList, 
                                         PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"), permute=FALSE)
    # rsNames = c("Estrogen_Receptor", loRegionAnno$filename[c(a549Ind, mcf7Ind)])
    rsNames = c("Estrogen_Receptor", loRegionAnno$filename[-sheff_dnaseInd])
    rsEnrichment_10[, rsNames:= rsNames]
    rsEnrichment_10
    
})
View(rsEnrichment_10[order(PC1,decreasing = TRUE)])

View(rsEnrichment_10[order(PC2,decreasing = TRUE)])
# PC2


###############################################
############# ~5% (35 patients) ###
set.seed(1234)
dataSplit = createDataPartition(y=factor(paste0(patientMetadata$ER_status,"_", patientMetadata$PGR_status)),
                                p = .05, list=FALSE)
trainingIDs = patientMetadata[dataSplit, subject_ID]
testIDs = patientMetadata[-dataSplit, subject_ID]
trainingMData = brcaMList[["methylProp"]][, 
                                          colnames(brcaMList[["methylProp"]]) %in% trainingIDs] 

# do the PCA
simpleCache("allMPCA_05", {
    prcomp(t(trainingMData), center = TRUE)
})
# plot(allMPCA$x[,c("PC1", "PC3")])
allMPCAWeights05 = as.data.table(allMPCA_05$rotation)
# mIQR = apply(trainingMData, 1, IQR)
# simpleCache("top10MPCA", {
#     prcomp(t(trainingMData[mIQR >= quantile(mIQR, 0.9), ]), center = TRUE)
# })
coordinates = brcaMList[["coordinates"]]
# top10Coord = coordinates[mIQR >= quantile(mIQR, 0.9), ]
# top10PCWeights = as.data.table(top10MPCA$rotation)
# top10TSNE = Rtsne(X = top10MPCA$x[, 1:50], pca = FALSE, max_iter=5000,
#                   perplexity = 30)
# plot(top10TSNE$Y)
# run PC region set enrichment analysis
simpleCache("rsEnrichment_05", {
    rsEnrichment_05 = pcRegionSetEnrichment(loadingMat=allMPCAWeights05, coordinateDT = coordinates, 
                                            GRList, 
                                            PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"), permute=FALSE)
    # rsNames = c("Estrogen_Receptor", loRegionAnno$filename[c(a549Ind, mcf7Ind)])
    rsNames = c("Estrogen_Receptor", loRegionAnno$filename[-sheff_dnaseInd])
    rsEnrichment_05[, rsNames:= rsNames]
    rsEnrichment_05
    
})
View(rsEnrichment_05[order(PC1,decreasing = TRUE)])

View(rsEnrichment_05[order(PC2,decreasing = TRUE)])

####### 20 patients




###############################################################################
# Summary of benchmarking
# PC1, top region sets (not counting different cell types as unique, just the factor/hist. mod.)
# 50%: Gata3, H3R17me2, ERa, Foxa1, AR, Ezh2, Znf217, Suz12, Tcf7l2, ...
# 10%: Gata3, H3R17me2, ERa, Foxa1, AR, Znf217, Tcf7l2, Jund, Cebpa, GR, ...
# 5%: 
# 20 s: 

## PC2
# 50%: 
# 10%: 
# 5%: 
# 20 s:

## PC3
# 50%: 
# 10%: 
# 5%: 
# 20 s:







# testing what happens if you include different subsets of patients
# each different PCA might give different region sets
# subset by clinical metadata