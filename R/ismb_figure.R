source(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/00-init.R"))


scriptID = "ismb"
# 
setwd(paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/plots/"))

set.seed(1234)

patientMetadata = brcaMetadata

# DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))
simpleCache("combinedBRCAMethyl_noXY")
brcaMList = combinedBRCAMethyl_noXY

# reading in the metadata, will be used to split data 
# into training and test set with balanced ER and PGR status
#restrict patients included in this analysis
patientMetadata = patientMetadata[patientMetadata$subject_ID %in% 
                                      colnames(brcaMList[["methylProp"]]), ]
# include all samples
trainingIDs = patientMetadata[, subject_ID]
trainingMData = brcaMList[["methylProp"]][, 
                                          colnames(brcaMList[["methylProp"]]) %in% trainingIDs] 

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
simpleCache("allMPCA_657", {
    prcomp(t(trainingMData), center = TRUE)
})
allMPCA = allMPCA_657
# plot(allMPCA$x[,c("PC1", "PC3")])
allMPCAWeights = as.data.table(allMPCA$rotation)
mIQR = apply(trainingMData, 1, IQR)
simpleCache("top10MPCA_657", {
    prcomp(t(trainingMData[mIQR >= quantile(mIQR, 0.9), ]), center = TRUE)
})
top10MPCA = top10MPCA_657
coordinates = brcaMList[["coordinates"]]
top10Coord = coordinates[mIQR >= quantile(mIQR, 0.9), ]
top10PCWeights = as.data.table(top10MPCA$rotation)

########################################################3333
# run PC region set enrichment analysis
simpleCache("rsEnrichment_657", {
    rsEnrichment = pcRegionSetEnrichment(loadingMat=allMPCAWeights, coordinateDT = coordinates, 
                                         GRList, 
                                         PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"), permute=FALSE)
    # rsNames = c("Estrogen_Receptor", loRegionAnno$filename[c(a549Ind, mcf7Ind)])
    rsNames = c("Estrogen_Receptor", loRegionAnno$filename[-sheff_dnaseInd])
    rsEnrichment[, rsNames:= rsNames]
    rsEnrichment
    
})
rsEnrichment = rsEnrichment_657
View(rsEnrichment[order(PC1,decreasing = TRUE)])

# check whether is enrichment is specific to this region set by
# seeing if loading values have a spike in the center of these region sets
# compared to surrounding genome 
GRList = lapply(GRList, resize, width = 14000, fix="center")

simpleCache("pcProf14k",{
    pcProf = pcEnrichmentProfile(loadingMat = allMPCAWeights, coordinateDT = coordinates,
                                 GRList=GRList, PCsToAnnotate = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                                 binNum = 21)
    # set names by reference
    setattr(pcProf, "names", names(GRList))
    pcProf
})
pcP = pcProf14k
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


# top 10% most variable
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
pcaWithAnno = cbind(as.data.table(allMPCA$x), patientMetadata)

colorByCols = colnames(patientMetadata)[!(colnames(patientMetadata) %in% "subject_ID")]
for (i in 2:6) {
    multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                           plotCols = c("PC1", paste0("PC", i)), 
                                           colorByCols=colorByCols)
    ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), paste0("multiColorPCAPlots1", i, "_657_", scriptID), 
                                    ".pdf"), plot = multiColorPCAPlots, device = "pdf",
                    limitsize=FALSE)
}
for (i in c(2, 4:6)) {
    multiColorPCAPlots = colorClusterPlots(pcaWithAnno, 
                                           plotCols = c("PC3", paste0("PC", i)), 
                                           colorByCols=colorByCols)
    ggplot2::ggsave(filename=paste0(Sys.getenv("PLOTS"), paste0("multiColorPCAPlots3", i, "_657_", scriptID), 
                                    ".pdf"), plot = multiColorPCAPlots, device = "pdf",
                    limitsize=FALSE)
}





######################################################################################
