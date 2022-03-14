library(data.table)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(stringr)
library(dendextend)
library(ComplexHeatmap)
library(gtools)
library(DESeq2)
library(scales)

load("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/VJ4001_D0-30_astrocyte_DDS.RData")
dds
as.data.frame(dds@colData)

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

genes <- c('BCL11B', 'TBR1', 'NEUROG2', 'EOMES', 'SATB2', 'CUX1', 'CALB2', 'RELN', 'FOXG1', 'NEUROD6', 'DCX', 'MAP2')

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genes]
length(found)
length(genes)
genes[!genes %in% found] #0

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genes)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)

names(mat)<- gsub('_',' ',coldata$condition)

mat <- mat[-c(1:7)]

names(mat)<- gsub('_',' ',coldata$condition[8:length(coldata$condition)])
names(mat) <- gsub('Lonza', 'LiPSC_GR1.1', names(mat))

cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red

scaled_mat <- t(scale(t(mat)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2)) # this will produce NaN's if there is no variance in a gene. must replace NaNs with 0 after to avoid an error next, or remove the gene.

h1 <- Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = TRUE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

h1

h2 <- Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

h2
