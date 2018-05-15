source(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/R/00-init.R"))


scriptID = "ismb"
# 
setwd(paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/"))
Sys.setenv("PLOTS"=paste0(Sys.getenv("PROCESSED"), "brca_PCA/analysis/plots/"))

set.seed(1234)

patientMetadata = brcaMetadata

# DNA methylation data
setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))
simpleCache("combinedBRCAMethyl_noXY", assignToVariable = "brcaMList")

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
# reading in the region sets
# load LOLA database
lolaPath1 = paste0(Sys.getenv("REGIONS"), "/LOLACore/hg38/")
regionSetDB = loadRegionDB(lolaPath1)
loRegionAnno = regionSetDB$regionAnno
lolaCoreRegionAnno = loRegionAnno
# a549Ind = grep("a549", loRegionAnno$cellType, ignore.case = TRUE)
sheff_dnaseInd = grep("sheffield_dnase", loRegionAnno$collection, ignore.case = TRUE)
# mcf7Ind = grep("mcf-7", loRegionAnno$cellType, ignore.case = TRUE)
# k562Ind = grep("k562", loRegionAnno$cellType,  ignore.case = TRUE)
# GRList = GRangesList(regionSetDB$regionGRL[c(a549Ind, mcf7Ind)])
fInd1 = filterFetal(lolaCoreRegionAnno)
lolaCoreRegionAnno = lolaCoreRegionAnno[-sort(unique(c(sheff_dnaseInd, fInd1)))]
GRList1 = GRangesList(regionSetDB$regionGRL[-sort(unique(c(sheff_dnaseInd, fInd1)))])

# ROADMAP Epigenome project and Jaspar motifs
lolaPath2 = paste0(Sys.getenv("REGIONS"), "/LOLAExt/hg38/")
regionSetDB2 = loadRegionDB(lolaPath2, useCache = TRUE)
loRegionAnno2 = regionSetDB2$regionAnno
roadmapRegionAnno = loRegionAnno2[loRegionAnno2$collection == "roadmap_epigenomics", ]
GRList2 = GRangesList(regionSetDB2$regionGRL[loRegionAnno2$collection == "roadmap_epigenomics"])
fInd2 = filterFetal(roadmapRegionAnno)
roadmapRegionAnno = roadmapRegionAnno[-fInd2]
GRList2 = GRList2[-fInd2]

# processing Jaspar motif regions
# resizing regions and filtering to only open chromatin in K562
# size was ~999
motifRegionAnno = loRegionAnno2[loRegionAnno2$collection == "jaspar_motifs", ]
GRList3 = GRangesList(regionSetDB2$regionGRL[loRegionAnno2$collection == "jaspar_motifs"])
# filtering based on chromatin accessibility data from blood/AML
load(paste0(Sys.getenv("PROCESSED"), "/aml_e3999/prjResources/","bloodAccessibleRegions.RData"))
GRList3 = getOLRegions(GRList = GRList3, intGR=bloodAccessibleRegions, removeOL = FALSE)
GRList3 = GRangesList(GRList3)
#I'm not sure that center is where motif is: GRList3 = resize(GRList3, width = 200, fix="center") 

# combine into one GRList
GRList = c(GRList1, GRList2, GRList3)

#adding ER Chipseq dataset
erSet = fread(paste0(Sys.getenv("CODE"), "PCARegionAnalysis/inst/extdata/",
                     "GSM2305313_MCF7_E2_peaks_hg38.bed"))
setnames(erSet, c("V1", "V2", "V3"), c("chr", "start", "end"))
GRList = c(GRangesList(MIRA:::dtToGr(erSet)), GRList)
# names(GRList) = c("GSM2305313_MCF7_E2_peaks_hg38.bed", names(GRList))

#################################################################
# cleaning up since there were many large objects
rm(list = c("GRList1", "GRList2", "GRList3", "regionSetDB", "regionSetDB2"))

#################################################################

allMPCAString = "pca_top100C_02cVM_pQC"
top10MPCAString = "pca_top10C_02cVM_pQC"

# gives output of rsEnrichment from PCA of all shared cytosines
# and rsEnrichmentTop10 from PCA of 10% most variable shared cytosines
source(paste0(Sys.getenv("CODE"),"/aml_e3999/src/PCRSA_pipeline.R"))

rsEnrichment = 
rsEnrichmentTop10 = 

write.csv(x = rsEnrichment, 
          file = dirData("analysis/sheets/PC_Enrichment_All_Shared_Cs.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(x = rsEnrichmentTop10, 
          file = dirData("analysis/sheets/PC_Enrichment_Top_10%_Variable_Cs.csv"),
          quote = FALSE, row.names = FALSE)

########################################################
# do the PCA
simpleCache("allMPCA_657", {
    prcomp(t(trainingMData), center = TRUE)
}, assignToVariable = "allMPCA")
# plot(allMPCA$x[,c("PC1", "PC3")])
allMPCAWeights = as.data.table(allMPCA$rotation)
mIQR = apply(trainingMData, 1, IQR)
simpleCache("top10MPCA_657", {
    prcomp(t(trainingMData[mIQR >= quantile(mIQR, 0.9), ]), center = TRUE)
}, assignToVariable = "top10MPCA")
coordinates = brcaMList[["coordinates"]]
top10Coord = coordinates[mIQR >= quantile(mIQR, 0.9), ]
top10PCWeights = as.data.table(top10MPCA$rotation)

########################################################

rsNames = c("GSM2305313_MCF7_E2_peaks_hg38.bed", 
            lolaCoreRegionAnno$filename, roadmapRegionAnno$filename, 
            motifRegionAnno$filename)
rsDescription = c("GSM2305313_MCF7_E2_peaks_hg38.bed",
                  lolaCoreRegionAnno$description, roadmapRegionAnno$description, 
                  motifRegionAnno$description)

# run PC region set enrichment analysis
simpleCache("rsEnrichment_657", {
    rsEnrichment = pcRegionSetEnrichment(loadingMat=allMPCAWeights, coordinateDT = coordinates, 
                                         GRList, 
                                         PCsToAnnotate = paste0("PC", 1:10), permute=FALSE)
    # rsNames = c("Estrogen_Receptor", loRegionAnno$filename[c(a549Ind, mcf7Ind)])
    rsEnrichment[, rsNames:= rsNames]
    rsEnrichment[, rsDescription := rsDescription]
    
}, recreate = TRUE)
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
