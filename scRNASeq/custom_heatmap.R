library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
library(RColorBrewer)
library(data.table)

load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData')

cluster.markers.0.2_merged <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Res0.2_merged_cluster_DE.csv')

clusters_up <- setorder(setDT(cluster.markers.0.2_merged), -avg_logFC)[, head(.SD, 10), keyby = cluster]

h1.seurat <- DoHeatmap(is019, features = c(clusters_up$gene), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)
#ggplot_build(h1.seurat)

#is019.scaledata <- ggplot_build(h1.seurat)$plot$data




#is019.scaledata <- as.data.frame(as.matrix(is019@assays$RNA@scale.data))
#is019.scaledata <- subset(is019.scaledata, row.names(is019.scaledata) %in% clusters_up$gene)
#is019.scaledata[,1:10]

is019.scaledata.gg <- ggplot_build(h1.seurat)

is019.scaledata.plotdata <- is019.scaledata.gg$plot$data

is019.scaledata.plotdata

data_wide <- dcast(is019.scaledata.plotdata, Cell +Identity ~ Feature, value.var="Expression")
head(data_wide)

data_wide <- na.omit(data_wide)

data_wide_t <- as.data.table(na.omit(data_wide))
data_wide_t <- as.data.table(t(data_wide_t[,-c(1:2)]))
data_wide_t[1:10,1:10]
names(data_wide_t) <- c(as.character(data_wide$Cell))
data_wide_t <- data.matrix(data_wide_t)
row.names(data_wide_t) <- names(data_wide)[3:length(data_wide)]
dim(data_wide_t)


data_wide_noNA <- na.omit(data_wide)

names(data_wide_t) <- data_wide_noNA$Cell

metadata <- as.data.frame(data_wide_noNA$Identity)
row.names(metadata) <- data_wide_noNA$Cell
names(metadata)[1] <- 'Cluster'

metadata$Cluster_new <- metadata$Cluster
metadata$Cluster_new <- gsub('CEPT0_2','CEPT1', metadata$Cluster_new)
metadata$Cluster_new <- gsub('CEPT3','CEPT2', metadata$Cluster_new)
metadata$Cluster_new <- gsub('CEPT_4-7','CEPT3', metadata$Cluster_new)
metadata$Cluster_new <- gsub('Y1','Y4', metadata$Cluster_new)
metadata$Cluster_new <- gsub('Y5','Y6', metadata$Cluster_new)
metadata$Cluster_new <- gsub('Y0_3','Y5', metadata$Cluster_new)
unique(metadata$Cluster_new)
metadata <- subset(metadata, select=c('Cluster_new'))
names(metadata)[1] <- 'Cluster'

metadata.copy <- metadat


ha <- HeatmapAnnotation(df = metadata, which='column', col=list(Cluster=c('CEPT1'='red', 'CEPT2'='orange', 'CEPT3'='lightgreen','Y4'='darkgreen','Y5'='skyblue', 'Y6'='purple', 'Y7'='#4b0082')))

draw(ha)

cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100)

h1 <- Heatmap(as.matrix(data_wide_t), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = F,
              cluster_rows = TRUE, cluster_columns = TRUE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Scaled expression'))

h1

save(h1, file='Custom_heatmap_h1_obj.RData')

#draw(h1 + ha, column_title = "IS019 CEPT, Y clusters and top upregulated genes")
draw(h1, annotation_legend_list = list(ha))


#column clustering
Heatmap(as.matrix(data_wide_t), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = F,
        cluster_rows = TRUE, cluster_columns = TRUE, top_annotation = ha,
        row_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Scaled expression'))

# no column clustering
Heatmap(as.matrix(data_wide_t), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = F,
        cluster_rows = TRUE, cluster_columns = FALSE, top_annotation = ha,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Scaled expression'))

# no column clustering, and order by cluster number.-----

unique(metadata$Cluster)

metadata.reorder <- metadata
metadata.reorder$cell <- row.names(metadata.reorder)
metadata.reorder <- as.data.table(metadata.reorder)
metadata.reorder[,original_order := 1:8950]
metadata.reorder <- metadata.reorder[order(Cluster)]
metadata.reorder
unique(metadata$Cluster)
metadata.reorder[,new_order := 1:8950]
#metadata.reorder <- metadata.reorder[,c('new_order')]
metadata.reorder <- as.data.frame(metadata.reorder)
row.names(metadata.reorder) <- metadata.reorder$cell
metadata.reorder <- subset(metadata.reorder, select=c('Cluster'))
metadata.reorder

data_wide_t_reorder <- subset(data_wide_t, select=c(metadata.reorder$cell))
all(names(data_wide_t_reorder) == metadata.reorder$cell)

ha_reorder <- HeatmapAnnotation(df = metadata.reorder, which='column', col=list(Cluster=c('CEPT1'='red', 'CEPT2'='orange', 'CEPT3'='lightgreen','Y4'='darkgreen','Y5'='skyblue', 'Y6'='purple', 'Y7'='#4b0082')))


Heatmap(as.matrix(data_wide_t_reorder), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = F,
        cluster_rows = TRUE, cluster_columns = FALSE, top_annotation = ha_reorder,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Scaled expression'))

