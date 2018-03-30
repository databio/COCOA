
library(LOLA)
library(simpleCache)
library(data.table)
library(GenomicRanges)
library(caret)
library(RGenomeUtils)
library(gridExtra) #marrangeGrob for colorClusterPlots()
# some of the environmental variables from aml/.../00-init.R will need to be reset
source(paste0(Sys.getenv("CODE"), "aml_e3999/src/00-init.R" ))
source(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/PRA.R"))


# 
setwd(paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/plots/"))
patientMetadata = fread(paste0(Sys.getenv("CODE"), 
                               "PCARegionAnalysis/metadata/brca_metadata.csv"))
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
# keep same proportion of patients with ER_status and PGR_status in training 
# and test sets
dataSplit = createDataPartition(y=factor(paste0(patientMetadata$ER_status,"_", patientMetadata$PGR_status)),
                                p = .5, list=FALSE)
trainingIDs = patientMetadata[dataSplit, subject_ID]
testIDs = patientMetadata[-dataSplit, subject_ID]
trainingMData = brcaMList[["methylProp"]][, 
                colnames(brcaMList[["methylProp"]]) %in% trainingIDs] 
testMData = brcaMList[["methylProp"]][, 
                colnames(brcaMList[["methylProp"]]) %in% testIDs]

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

########################################################3333
# do the PCA
simpleCache("allMPCA", {
    prcomp(t(trainingMData), center = TRUE)
})
# plot(allMPCA$x[,c("PC1", "PC3")])
allMPCAWeights = as.data.table(allMPCA$rotation)
mIQR = apply(trainingMData, 1, IQR)
simpleCache("top10MPCA", {
    prcomp(t(trainingMData[mIQR >= quantile(mIQR, 0.9), ]), center = TRUE)
})
coordinates = brcaMList[["coordinates"]]
top10Coord = coordinates[mIQR >= quantile(mIQR, 0.9), ]
top10PCWeights = as.data.table(top10MPCA$rotation)
# top10TSNE = Rtsne(X = top10MPCA$x[, 1:50], pca = FALSE, max_iter=5000,
#                   perplexity = 30)
# plot(top10TSNE$Y)
# run PC region set enrichment analysis
rsEnrichment = pcRegionSetEnrichment(loadingMat=allMPCAWeights, coordinateDT = coordinates, 
                      GRList, 
                      PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"), permute=FALSE)
# rsNames = c("Estrogen_Receptor", loRegionAnno$filename[c(a549Ind, mcf7Ind)])
rsNames = c("Estrogen_Receptor", loRegionAnno$filename[-sheff_dnaseInd])
rsEnrichment[, rsNames:= rsNames]
View(rsEnrichment[order(PC1,decreasing = TRUE)])

# check whether is enrichment is specific to this region set by
# seeing if loading values have a spike in the center of these region sets
# compared to surrounding genome 
GRList = lapply(GRList, resize, width = 10000, fix="center")

pcProf = pcEnrichmentProfile(loadingMat = allMPCAWeights, coordinateDT = coordinates,
                             GRList=GRList, PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                             binNum = 21)
pcP = pcProf
# plot(pcP$PC1, type="l")
# plot(pcP$PC2, type="l")
# plot(pcP$PC3, type="l")
# plot(pcP$PC4, type="l")
# plot(pcP$PC5, type="l")

rsNames = c("Estrogen_Receptor", loRegionAnno$filename[c(a549Ind, mcf7Ind)])

grDevices::pdf(paste0(Sys.getenv("PLOTS"), "allMPCProfilesDNase300.pdf"))
for (i in 1:length(pcProf)) {
    plot(pcP[[i]]$PC1, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC2, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC3, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC4, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC5, type="l") + title(rsNames[i])
}
dev.off()


# top 10% mofst variable
pcProf10 = pcEnrichmentProfile(loadingMat = top10PCWeights, coordinateDT = top10Coord,
                    GRList=GRList, PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                    binNum = 21)
pcP = pcProf10
grDevices::pdf(paste0(Sys.getenv("PLOTS"), "top10MPCProfiles.pdf"))
for (i in 1:length(pcProf)) {
    plot(pcP[[i]]$PC1, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC2, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC3, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC4, type="l") + title(rsNames[i])
    plot(pcP[[i]]$PC5, type="l") + title(rsNames[i])
}
dev.off()

# initResults = list(rsEnrichment, rsEnrichment2, pcProf, ret[c(46:72, 623:634)])
# save(initResults, file=paste0(Sys.getenv("PROCESSED"), "brca_PCA/", "initialPRAresults.RData"))







# confirm that results make sense
# PCs that are enriched for ATACseq from a certain cell type should be 
# able to separate that cell type from the others
# visualize


######################################################
# visualizing PCA

# shorten metadata for visualization
patientMetadata$menopause_status = sub(pattern = "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)", 
            "Pre", x = patientMetadata$menopause_status, fixed = TRUE)
patientMetadata$menopause_status = sub(pattern = "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)", 
                                       "Post", x = patientMetadata$menopause_status, fixed = TRUE)
patientMetadata$menopause_status = sub(pattern = "Indeterminate (neither Pre or Postmenopausal)", 
                                       "Indeterminate", x = patientMetadata$menopause_status, fixed = TRUE)
patientMetadata$menopause_status = sub(pattern = "Peri (6-12 months since last menstrual period)", 
                                       "Peri", x = patientMetadata$menopause_status, fixed = TRUE)
patientMetadata$race = sub(pattern = "BLACK OR AFRICAN AMERICAN", 
                                       "BLACK", x = patientMetadata$race, fixed = TRUE)
patientMetadata$race = sub(pattern = "AMERICAN INDIAN OR ALASKA NATIVE", 
                           "NATIVE AM. OR ALASKAN", x = patientMetadata$race, fixed = TRUE)

# add annotation information
pcaWithAnno = cbind(as.data.table(allMPCA$x), patientMetadata[dataSplit, ])

colorByCols = colnames(patientMetadata)[!(colnames(patientMetadata) %in% "subject_ID")]
for (i in 2:6) {
    multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                             plotCols = c("PC1", paste0("PC", i)), 
                                             colorByCols=colorByCols)
    ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), paste0("multiColorPCAPlots1", i), 
                                    ".pdf"), plot = multiColorPCAPlots, device = "pdf",
                    limitsize=FALSE)
}
for (i in c(2, 4:6)) {
    multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                           plotCols = c("PC3", paste0("PC", i)), 
                                           colorByCols=colorByCols)
    ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), paste0("multiColorPCAPlots3", i), 
                                    ".pdf"), plot = multiColorPCAPlots, device = "pdf",
                    limitsize=FALSE)
}




######################################################################################
# analysis of whether certain subsets of region sets are the variable ones

# erSet = GRangesList(MIRA:::dtToGr(erSet))
# genomeLoadings = cbind(coordinates, as.data.table(allMPCA$rotation))
# regionAv = averageByRegion(BSDT = genomeLoadings[, .(chr, start, PC1, PC2, PC3, PC4)], regionsGRL = erSet, 
#                 jCommand = MIRA:::buildJ(c("PC1", "PC2", "PC3", "PC4"), "mean"),
#                 hasCoverage = FALSE)
simpleCache("regionAv101", {
    regionAv = lapply(GRList[1:100], function(x) averageByRegion(loadingMat = allMPCA$rotation, coordinateDT= coordinates, GRList = x, 
                                                                  PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")))
    names(regionAv) <- names(GRList)[1:100]
    regionAv
})
grep(pattern = "A549P300", names(regionAv))
# two region sets that show small peaks in the middle of dips
# these have some regions that have very high loadings which probably drive that
hist(regionAv$wgEncodeAwgTfbsHaibA549P300V0422111Etoh02UniPk.narrowPeak$PC5)
hist(regionAv$wgEncodeAwgTfbsHaibA549Tcf12V0422111Etoh02UniPk.narrowPeak$PC2)
# a region set with a strong peak, has a long tail but very robust
hist(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`$PC6)
# are the high regions in one PC the same as the high regions in another PC?
View(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[order(PC6, decreasing = TRUE), ])
plot(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[, .(PC1, PC3)])
plot(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[, .(PC1, PC4)])
plot(regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[, .(PC3, PC4)])
cor(x = regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[ PC1 > .002, PC1],
    y = regionAv$`Human_MCF-7_GATA3_No-treatment_White.bed`[ PC1 > .002, PC3])

# LOLA analysis of regions that are high in one PC compared to another for 
# Human_MCF−7_ESR1_E2−45min_Brown.bed which had good peaks for PC1, 3, 4, 5
plot(regionAv$`Human_MCF-7_ESR1_E2-45min_Brown.bed`[, .(PC1, PC3)])
esr1 = regionAv$`Human_MCF-7_ESR1_E2-45min_Brown.bed`
pc1Reg = esr1[PC1 > 0.001, .(chr, start, end)]
pc2Reg = esr1[PC2 > 0.001, .(chr, start, end)]
pc3Reg = esr1[PC3 > 0.001, .(chr, start, end)]
RGenomeUtils::writeBed(esr1[, .(chr, start, end)], filename = "esr1_all.bed")
RGenomeUtils::writeBed(pc1Reg, filename = "esr1_PC1.bed")
RGenomeUtils::writeBed(pc2Reg, filename = "esr1_PC2.bed")
RGenomeUtils::writeBed(pc3Reg, filename = "esr1_PC3.bed")
fos = regionAv$`Human_MCF-7_c-Fos_E2-45min-3hr_Liu.bed`
pc1Reg = fos[PC1 > 0.001, .(chr, start, end)]
pc2Reg = fos[PC2 > 0.001, .(chr, start, end)]
pc3Reg = fos[PC3 > 0.001, .(chr, start, end)]
RGenomeUtils::writeBed(fos[, .(chr, start, end)], filename = "fos_all.bed")
RGenomeUtils::writeBed(pc1Reg, filename = "fos_PC1.bed")
RGenomeUtils::writeBed(pc2Reg, filename = "fos_PC2.bed")
RGenomeUtils::writeBed(pc3Reg, filename = "fos_PC3.bed")
# findOverlaps(dtToGr(pc1RegEsr), dtToGr(pc1RegFos))
    
##################################################################################
# could also do the analysis for hormone receptors in breast cancer:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4754852/
load("allBRCAexpression.RData")
# figuring out which patients are ER+ (inexact estimate based on gene expression data)
ER1Ind = grep(pattern = "ENSG00000091831", x = myExprDT[, "Gene"],
     ignore.case = TRUE) # ESR1
hist(as.numeric(myExprDT[ER1Ind, 2:ncol(myExprDT)]),
     breaks = seq(0,300,1))
sum(as.numeric(myExprDT[ER1Ind, 2:ncol(myExprDT)]) < 1) / (ncol(myExprDT)-1)
# progesterone receptoer positive
PGRInd = grep(pattern = "ENSG00000082175", x = myExprDT[, "Gene"],
              ignore.case = TRUE) # PGR
hist(as.numeric(myExprDT[PGRInd, 2:ncol(myExprDT)]),
     breaks = seq(0,400,1))
sum(as.numeric(myExprDT[PGRInd, 2:ncol(myExprDT)]) < 2) / (ncol(myExprDT)-1)


# References
# https://www.ncbi.nlm.nih.gov/pubmed/17616709/: ER transcriptional network
