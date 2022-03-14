require(corrplot)
require(RColorBrewer)
require(kmed)

setwd("Z:/NGS_related/BulkRNA/ISB017/Analysis/")
mat = as.matrix(read.csv("DE/ISB017_merged_for_partek.csv", row = 1, as.is = T, check.names = F))
mat = mat[, !grepl(pattern = "^2D.*", x = colnames(mat))]
corrMat = cor(na.omit(mat))

ColorPalette = colorRampPalette(
  colors = c(
    "#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", 
    "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"
  )
)

# CEPT Corrplot
mat_subset = mat[, grepl(pattern = "^CEPT.*", x = colnames(mat))]
corrMat <- cor(na.omit(mat_subset))

groups = gsub(pattern = "_[0-9]+$", replacement = "", x = colnames(mat_subset))
sampCols = c("black", "blue", "red")
names(sampCols) = unique(groups)
textCols = sampCols[groups]

pdf("Correlation/ISB017_CorrPlot_3D_CEPT_CI95_03062020.pdf", width = 12, height = 12)
corrplot(
  corr = corrMat^2, 
  type = "upper",
  bg = "lightgrey", 
  col = rev(ColorPalette(60)), 
  insig = "blank",
  cl.lim = c(0.5, 1), 
  title = bquote("ISB017 - CEPT | " ~ italic(R)^2),
  mar = c(0, 0, 4, 0),
  tl.col = textCols
)
dev.off()

# Y27 Corrplot
mat_subset = mat[, grepl(pattern = "^Y.*", x = colnames(mat))]
corrMat <- cor(na.omit(mat_subset))

groups = gsub(pattern = "_[0-9]+$", replacement = "", x = colnames(mat_subset))
sampCols = c("black", "blue", "red")
names(sampCols) = unique(groups)
textCols = sampCols[groups]

pdf("ISB017_CorrPlot_3D_Y_CI95_03062020.pdf", width = 12, height = 12)
corrplot(
  corr = corrMat^2, 
  type = "upper",
  bg = "lightgrey", 
  col = rev(ColorPalette(60)), 
  insig = "blank",
  cl.lim = c(0.5,1), 
  title = bquote("ISB017 - Y27 | " ~ italic(R)^2),
  mar = c(0,0,4,0),
  tl.col = textCols
)
dev.off()

# Silhouette plot
corrMat = cor(na.omit(mat))
groups = gsub(pattern = "_[0-9]+$", replacement = "", x = colnames(mat))
distMat = 1-corrMat

# Return medoid indices from distance matrix
Meds = function(distMat, groups) {
  # distMat[lower.tri(distMat, diag = T)] = NA 
  diag(distMat) = NA
  uniqueGroups = unique(groups)
  names(groups) = colnames(distMat)
  meds = vector("integer", length = length(uniqueGroups))
  sampleIndices = 1:ncol(distMat)
  names(sampleIndices) = colnames(distMat)
  for(i in 1:length(uniqueGroups)) {
    distTemp = distMat[names(groups)[groups == uniqueGroups[i]], 
                       names(groups)[groups == uniqueGroups[i]]]
    rowSumsTemp = rowSums(distTemp, na.rm = T)
    meds[i] = sampleIndices[rownames(distTemp)[rowSumsTemp == min(rowSumsTemp)][1]]
  }
  return(meds)
}

idMed = Meds(distMat = distMat, groups = groups)
uniqueGroups = 1:length(unique(groups))
names(uniqueGroups) = unique(groups)
idClust = uniqueGroups[groups]
silResult = sil(distMat, idmedoid = idMed, idcluster = idClust,
                title = "ISB017 - Silhouette Widths")
silResult$plot$data$cluster = groups
pdf("Correlation/ISB017_Sil_03062020.pdf", width = 10, height = 8)
silResult$plot
dev.off()
