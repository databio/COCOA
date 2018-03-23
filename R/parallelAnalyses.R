
# parallel/confirmatory analyses
library(data.table)

# run cytosines with extreme loading values through LOLA
setCacheDir(paste0(Sys.getenv("PROCESSED"), "brca_PCA/RCache/"))
simpleCache("combinedBRCAMethyl_noXY")
coordinates = combinedBRCAMethyl_noXY[["coordinates"]]

simpleCache("allMPCA")

mLoadings = allMPCA$rotation[, 1:6]

hist(mLoadings[,"PC1"])
pc1 = abs(mLoadings[, "PC1"])
# highest 25% of cytosines
pc1KeepCoords = coordinates[pc1 >= quantile(pc1, probs = 0.75), ]
pc1KeepCoords[, end := start]
write.table(x=pc1KeepCoords, file = "PC1.bed", quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)

#PC2
pc2KeepCoords = coordinates[pc1 >= quantile(pc1, probs = 0.75), ]
pc2KeepCoords[, end := start]
write.table(x=pc1KeepCoords, file = "PC1.bed", quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)



# correlate MIRA scores with principal component values for relevant components
# if there is a peak for estrogen receptor for PC1, correlate PC1 values
# with MIRA values for that same region set



# correlate expression of transcription factors with principal component values
# eg. if estrogent r