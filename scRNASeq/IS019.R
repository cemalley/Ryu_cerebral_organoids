library(Seurat)
library(data.table)
library(biomaRt)
library(stringr)
library(ggplot2)

setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Countfiles_ENSG')

files <- Sys.glob('*csv')
for (file in files){
  dt <- as.data.frame(fread(file))
  names(dt)[1] <- 'ENSG'
  dt <- as.data.table(dt)
  ensg.genes <- data.table("ENSG" = dt$ENSG)
  genes <- ensg.genes$ENSG
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="uswest.ensembl.org", ensemblRedirect = FALSE)
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
  G_list <- as.data.table(G_list)
  ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG", by.y="ensembl_gene_id")
  ensg.genes <- na.omit(ensg.genes)
  dt <- subset(dt, dt$ENSG %in% ensg.genes$ENSG)
  dt <- merge(dt, ensg.genes, by="ENSG")
  dt <- dt[,ENSG:=NULL]
  dt <- dt[!duplicated(external_gene_name),]
  
  sample_id <- str_split_fixed(file, "-dense-expression-matrix.csv",2)[1]
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/Chromium/IS019/Countfiles_gene_symbol/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
  
}

setwd('../Countfiles_gene_symbol/')
files <- Sys.glob('*gene_symbol.csv')

for (file in files){
  dt <- fread(file, header=T)
  cols <- names(dt)[c(1: (length(names(dt))-1) )]
  dt<- subset(dt, select=c("external_gene_name",cols))
  sample_id <- str_split_fixed(file, "_gene_symbol.csv",2)[1]
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/Chromium/IS019/Countfiles_gene_symbol/", sample_id, "_gene_symbol.csv"), col.names = T, row.names = F, quote=F, sep=",")
}


reformat_for_seurat <- function(x, samplename){
  x <- x[!duplicated(external_gene_name),]
  x <- x[external_gene_name != "NA",]
  x <- subset(x, select=c("external_gene_name", unlist(names(x))[2:(length(x)-1)] ))
  x <- as.data.frame(x)
  row.names(x) <- x$external_gene_name
  x <- x[-1]
  barcodes <- names(x)
  barcodes <- paste0(barcodes, paste0(".", samplename))
  names(x) <- barcodes
  x <- CreateSeuratObject(counts = x, project = "IS019")
  x@meta.data$sample <- samplename
  # x <- NormalizeData(x)
  # x <- FindVariableFeatures(x)
  # x <- ScaleData(x)
  # x <- RunPCA(x)
  # x <- FindNeighbors(x)
  # x <- FindClusters(x)
  # x <- RunTSNE(x)
  return(x)
}

is019.CEPT <- reformat_for_seurat(fread(files[1]), "CEPT")
is019.CE <- reformat_for_seurat(fread(files[2]), "CE")
is019.Y <- reformat_for_seurat(fread(files[3]), "Y")


is019 <- merge(is019.CEPT, is019.CE, normalize=F)
is019 <- merge(is019, is019.Y, normalize=F)

is019 <- NormalizeData(is019)
is019 <- FindVariableFeatures(is019)
is019 <- ScaleData(is019, features=c(row.names(as.data.frame(as.matrix(is019@assays$RNA@data)))))
is019 <- RunPCA(is019)
is019 <- FindNeighbors(is019)
is019 <- FindClusters(is019)
is019 <- RunTSNE(is019)

save(is019, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019.Seurat.RData')

Idents(is019) <- 'sample'

DoHeatmap(is019, features = c('NEUROD6', 'DCX', 'NEUROD2', 'FOXG1', 'MAP2', 'FABP7', 'GPM6A', 'BASP1', 'PAFAH1B3', 'FDFT1', 'CRMP1', 'FEZ1', 'NCAM1', 'BLOC1S1', 'STMN4', 'STMN2', 'MAPT', 'GPM6B', 'SOX11', 'TUBB2B'),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

CEPT.markers <- FindMarkers(is019, ident.1=c('CEPT'), ident.2=c('Y'))
write.table(CEPT.markers, '../Analysis/CEPT.markers.csv', row.names = T, col.names = T)


DoHeatmap(is019, features = c(row.names(CEPT.markers)[1:20]),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)



# heatmap july 10-----

genelist <- c('NEUROD6', 'DCX', 'NEUROD2', 'FOXG1', 'MAP2', 'FABP7', 'GPM6A', 'BASP1', 'PAFAH1B3', 'TUBB', 'CRMP1', 'FEZ1', 'NCAM1', 'BLOC1S1', 'STMN4', 'STMN2', 'MAPT', 'GPM6B', 'SOX11', 'TUBB2B', 'TUBB2A', 'NFIA', 'NSG1', 'SYT1', 'MAP1B')
load("/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019.Seurat.RData")

library(Seurat)
library(RColorBrewer)
library(ggplot2)

Idents(is019) <- 'sample'

DoHeatmap(is019, features = c(genelist),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

# clustering, sept 12-----
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/')
setwd('../Countfiles_gene_symbol/')
files <- Sys.glob('*gene_symbol.csv')

is019.CEPT <- reformat_for_seurat(fread(files[1]), "CEPT")
is019.Y <- reformat_for_seurat(fread(files[3]), "Y")

is019 <- merge(is019.CEPT, is019.Y)
is019 <- NormalizeData(is019)
is019 <- FindVariableFeatures(is019)
is019 <- ScaleData(is019, features=c(row.names(as.data.frame(as.matrix(is019@assays$RNA@data)))))
is019 <- RunPCA(is019)
is019 <- FindNeighbors(is019)
is019 <- FindClusters(is019)
is019 <- RunTSNE(is019)


unique(is019@meta.data$sample)

TSNEPlot(is019)
Idents(is019) <- 'sample'
TSNEPlot(is019)

is019 <- FindClusters(is019, resolution = 0.1)
TSNEPlot(is019)

all.markers <- FindAllMarkers(is019)
all.markers <- as.data.table(all.markers)
all.markers <- all.markers[abs(avg_logFC) > 1 ,]
all.markers <- all.markers[,c(7,1:6)]
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/')
fwrite(all.markers,'IS019_clustering_markers.csv', sep=',', col.names = T, row.names = F, quote=T)

plot.genes <- all.markers[abs(avg_logFC) >= 2,'gene']
plot.genes <- plot.genes$gene

DoHeatmap(is019, features = c(plot.genes),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

Idents(is019) <- 'sample'
DoHeatmap(is019, features = c(plot.genes),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

is019@meta.data$label <- paste0(is019@meta.data$sample,'_', is019@meta.data$RNA_snn_res.0.1)

Idents(is019) <- 'label'
DoHeatmap(is019, features = c(plot.genes),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

FeaturePlot(is019, features=c(plot.genes))

#new cluster numbers----
is019 <- FindClusters(is019, resolution = 0.209)
TSNEPlot(is019)

save(is019, file='../IS019_CEPT_Y_only.Seurat.RData')

load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only.Seurat.RData')

is019 <- FindClusters(is019, resolution = 0.2)
TSNEPlot(is019)

Idents(is019) <- 'RNA_snn_res.0.8'
TSNEPlot(is019)
Idents(is019) <- 'RNA_snn_res.0.2'
TSNEPlot(is019)

Idents(is019) <- 'sample'
TSNEPlot(is019)


genes <- c('TTR', 'NEUROD6', 'TPH1', 'APOE', 'IGFBP7', 'GNB3', 'SERPINF1', 'RBP1', 'LGALS1', 'NEUROD2', 'SOX11', 'ID3', 'NFIB', 'GSG1', 'PTGDS', 'NDUFA4L2', 'FOXG1', 'IFITM3', 'COL3A1', 'MGP', 'SOX4', 'DCX', 'ID1', 'NEUROD1', 'LINC01551', 'COL1A1', 'PCAT4', 'COL1A2', 'CRABP2', 'BCL11A', 'NEAT1', 'SLA', 'CLU', 'DCN', 'ISOC1', 'C1orf61', 'TMEM161B-AS1', 'DCT', 'SPARCL1', 'FDFT1', 'RCVRN', 'TCF4', 'S100B', 'MEF2C', 'FABP7', 'SST', 'CSRP2', 'TIMP1', 'BASP1', 'ANXA2', 'TUBB', 'OTX2', 'LUM', 'SFRP1', 'GPM6A', 'NEFL', 'CD63', 'MAP1B')

Idents(is019) <- 'RNA_snn_res.0.2'
DoHeatmap(is019, features = c(genes),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

cluster.markers.0.2 <- FindAllMarkers(is019)

cluster.markers.0.2
cluster.markers.0.2 <- as.data.table(cluster.markers.0.2)
fwrite(cluster.markers.0.2, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Res0.2_cluster_DE.csv')

clusters_up <- setorder(setDT(cluster.markers.0.2), -avg_logFC)[, head(.SD, 25), keyby = cluster]
clusters_down <- setorder(setDT(cluster.markers.0.2), avg_logFC)[, head(.SD, 25), keyby = cluster]

is019@meta.data$res0.2_labels <- paste0(is019@meta.data$sample,is019@meta.data$RNA_snn_res.0.2)

Idents(is019) <- 'res0.2_labels'

DoHeatmap(is019, features = c(clusters_up$gene), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

DoHeatmap(is019, features = c(clusters_down$gene), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

is019@meta.data$res0.2_labels_merge <- is019@meta.data$res0.2_labels
is019@meta.data$res0.2_labels_merge <- gsub('CEPT2', 'CEPT_new', is019@meta.data$res0.2_labels_merge)
is019@meta.data$res0.2_labels_merge <- gsub('CEPT0', 'CEPT_new', is019@meta.data$res0.2_labels_merge)
is019@meta.data$res0.2_labels_merge <- gsub('CEPT_new', 'CEPT0_2', is019@meta.data$res0.2_labels_merge)

Idents(is019) <- 'res0.2_labels_merge'

TSNEPlot(is019)

CellSelector(TSNEPlot(is019))

is019@meta.data$res0.2_labels_merge <- gsub('CEPT4', 'CEPT_new', is019@meta.data$res0.2_labels_merge)
is019@meta.data$res0.2_labels_merge <- gsub('CEPT5', 'CEPT_new', is019@meta.data$res0.2_labels_merge)
is019@meta.data$res0.2_labels_merge <- gsub('CEPT6', 'CEPT_new', is019@meta.data$res0.2_labels_merge)
is019@meta.data$res0.2_labels_merge <- gsub('CEPT7', 'CEPT_new', is019@meta.data$res0.2_labels_merge)
is019@meta.data$res0.2_labels_merge <- gsub('CEPT_new', 'CEPT_other', is019@meta.data$res0.2_labels_merge)
is019@meta.data$res0.2_labels_merge <- gsub('CEPT_other', 'CEPT_4-7', is019@meta.data$res0.2_labels_merge)


is019@meta.data$res0.2_labels_merge <- gsub('Y0', 'Y0_3', is019@meta.data$res0.2_labels_merge)
is019@meta.data$res0.2_labels_merge <- gsub('Y3', 'Y0_3', is019@meta.data$res0.2_labels_merge)

unique(is019@meta.data$res0.2_labels_merge)

Idents(is019)<-'res0.2_labels_merge'

save(is019, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData')

TSNEPlot(is019)

cluster.markers.0.2_merged <- FindAllMarkers(is019)

fwrite(cluster.markers.0.2_merged, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Res0.2_merged_cluster_DE.csv')

clusters_up <- setorder(setDT(cluster.markers.0.2_merged), -avg_logFC)[, head(.SD, 10), keyby = cluster]
clusters_down <- setorder(setDT(cluster.markers.0.2_merged), avg_logFC)[, head(.SD, 10), keyby = cluster]

clusters_bigFC <- cluster.markers.0.2_merged[order(-abs(avg_logFC)),head(.SD,60)]
clusters_bigFC <- clusters_bigFC[order(-avg_logFC),]


VlnPlot(is019, features = c('TTR'))


h1.seurat <- DoHeatmap(is019, features = c(clusters_up$gene), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)
ggplot_build(h1.seurat)

is019.scaledata <- ggplot_build(h1.seurat)$plot$data

DoHeatmap(is019, features = c(clusters_down$gene), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

DoHeatmap(is019, features = c(clusters_bigFC$gene), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

# making the same heatmap but with complexheatmap----
is019.scaledata <- as.data.frame(as.matrix(is019@assays$RNA@scale.data))
is019.scaledata <- subset(is019.scaledata, row.names(is019.scaledata) %in% clusters_up$gene)
is019.scaledata[,1:10]

is019.scaledata.gg <- ggplot_build(h1.seurat)

is019.scaledata.plotdata <- is019.scaledata.gg$plot$data

is019.scaledata.plotdata

data_wide <- dcast(is019.scaledata.plotdata, Cell +Identity ~ Expression, value.var="Feature")
head(data_wide)

metadata <- as.data.frame(is019@meta.data)
metadata <- subset(metadata, select=c('res0.2_labels_merge'))
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


ha <- HeatmapAnnotation(df = metadata, which='row', width=unit(1, 'cm'), col=list(Cluster=c('CEPT1'='red', 'CEPT2'='orange', 'CEPT3'='lightgreen','Y4'='darkgreen','Y5'='skyblue', 'Y6'='purple', 'Y7'='#4b0082')))

draw(ha)

cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100)

h1 <- Heatmap(as.matrix(is019.scaledata), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = F,
              cluster_rows = TRUE, cluster_columns = TRUE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Scaled expression'))

h1



# subpopulation analysis for CEPT only------

is019.CEPT <- reformat_for_seurat(fread(files[1]), "CEPT")

finish_seurat <- function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x <- ScaleData(x, features = row.names(as.data.frame(x@assays$RNA@data)))
  Idents(x) <- 'sample'
  x <- RunPCA(x)
  x <- FindNeighbors(x)
  x <- FindClusters(x)
  x <- RunTSNE(x) # error for Noci sample.
  return(x)
}

is019.CEPT <- finish_seurat(is019.CEPT)


clusters.info <- as.data.frame(is019@meta.data)
clusters.info$barcode <- row.names(clusters.info)
clusters.info <- as.data.table(clusters.info)
clusters.info <- clusters.info[sample=='CEPT',]
clusters.info <- clusters.info[,c('barcode','res0.2_labels_merge')]
clusters.info

all(row.names(is019.CEPT@meta.data) %in% clusters.info$barcode)

is019.CEPT@meta.data$res0.2_labels_merge <- clusters.info$res0.2_labels_merge

Idents(is019.CEPT) <- 'res0.2_labels_merge'

TSNEPlot(is019.CEPT)

save(is019.CEPT, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019.CEPT.Seurat.RData')


load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019.CEPT.Seurat.RData')

is019.CEPT <- RunUMAP(is019.CEPT, reduction = "pca", dims=1:5)

save(is019.CEPT, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019.CEPT.Seurat.RData')

load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019.CEPT.Seurat.RData')

UMAPPlot(is019.CEPT)


cept3.vs.0_2 <- FindMarkers(is019.CEPT, ident.1 = 'CEPT3', ident.2='CEPT0_2')
fwrite(cept3.vs.0_2, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/CEPT3.vs.0_2.csv')

cept3.vs.0_2$gene <- row.names(cept3.vs.0_2)
cept3.vs.0_2 <- as.data.table(cept3.vs.0_2)

cept3.vs.0_2.up <- cept3.vs.0_2[order(-avg_logFC),head(.SD,25)]
cept3.vs.0_2.down <- cept3.vs.0_2[order(avg_logFC),head(.SD,25)]


DoHeatmap(is019.CEPT, features = c(cept3.vs.0_2.up$gene, cept3.vs.0_2.down$gene), raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


#looking for neuronal/cortical genes in IS019 CEPT only-----
neural.genes <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/marker_sets/glia neuron gene list.xlsx'))

neural.genes <- neural.genes[,c(1:2)]
neural.genes

cluster.averages <- AverageExpression(is019.CEPT, return.seurat = FALSE)
cluster.averages <- as.data.frame(cluster.averages)
names(cluster.averages) <- gsub('RNA.', 'Mean.', names(cluster.averages))
cluster.averages$GeneId <- row.names(cluster.averages)

cortical.avgs <- AverageExpression(is019.CEPT, features = unlist(neural.genes[category=='Cortical neurons',gene]) )
cortical.avgs.sums <- colSums(cortical.avgs$RNA)
cortical.avgs.sums


#all.genes.present <- row.names(as.data.frame(as.matrix(is019.CEPT@assays$RNA@counts)))
all.genes.present <- row.names(as.data.frame(as.matrix(is019@assays$RNA@counts)))
progenitor.genes <- neural.genes[category=='Progenitors',gene]

progenitor.genes <- progenitor.genes[progenitor.genes %in% all.genes.present]


progenitor.avgs <- AverageExpression(is019.CEPT, features = progenitor.genes)
progenitor.avgs.sums <- colSums(progenitor.avgs$RNA)
progenitor.avgs.sums


glia.genes <- neural.genes[category=='Glia',gene]
glia.genes <- glia.genes[glia.genes %in% all.genes.present]
glia.avgs <- AverageExpression(is019.CEPT, features = glia.genes)
glia.avgs.sums <- colSums(glia.avgs$RNA)
glia.avgs.sums

lower.cortex.genes <- neural.genes[category=='Lower cortex',gene]
lower.cortex.genes <- lower.cortex.genes[lower.cortex.genes %in% all.genes.present]
lower.cortex.avgs <- AverageExpression(is019.CEPT, features = lower.cortex.genes)
lower.cortex.avgs.sums <- colSums(lower.cortex.avgs$RNA)
lower.cortex.avgs.sums

upper.cortex.genes <- neural.genes[category=='Upper cortex',gene]
upper.cortex.genes <- upper.cortex.genes[upper.cortex.genes %in% all.genes.present]
upper.cortex.avgs <- AverageExpression(is019.CEPT, features = upper.cortex.genes)
upper.cortex.avgs.sums <- colSums(upper.cortex.avgs$RNA)
upper.cortex.avgs.sums

other.genes <- neural.genes[category=='Other',gene]
other.genes <- other.genes[other.genes %in% all.genes.present]
other.avgs <- AverageExpression(is019.CEPT, features = other.genes)
other.avgs.sums <- colSums(other.avgs$RNA)
other.avgs.sums

neural.crest.genes <- neural.genes[category=='Neural crest',gene]
neural.crest.genes <- neural.crest.genes[neural.crest.genes %in% all.genes.present]
neural.crest.avgs <- AverageExpression(is019.CEPT, features = neural.crest.genes)
neural.crest.avgs.sums <- colSums(neural.crest.avgs$RNA)
neural.crest.avgs.sums

cortical.genes <- neural.genes[category=='Cortical neurons',gene]
cortical.genes <- cortical.genes[cortical.genes %in% all.genes.present]
cortical.avgs <- AverageExpression(is019.CEPT, features = cortical.genes)
cortical.avgs.sums <- colSums(cortical.avgs$RNA)
cortical.avgs.sums

# multiple tSNE of cell gene means for only CEPT-----

tsne.blend <- as.data.frame(is019.CEPT@reductions$tsne@cell.embeddings)
tsne.blend$barcode <- row.names(tsne.blend)
tsne.blend <- as.data.table(tsne.blend)

tsne.blend

is019.CEPT.data <- as.data.frame(is019.CEPT@assays$RNA@data)
is019.CEPT.data$GeneId <- row.names(is019.CEPT.data)
is019.CEPT.data <- as.data.table(is019.CEPT.data)
glia.avgs <- colMeans(is019.CEPT.data[GeneId %in% glia.genes,-c('GeneId')])


tsne.blend[,glia:=colMeans(is019.CEPT.data[GeneId %in% glia.genes,-c('GeneId')])]
tsne.blend[,lower_cortical:=colMeans(is019.CEPT.data[GeneId %in% lower.cortex.genes,-c('GeneId')])]
tsne.blend[,upper_cortical:=colMeans(is019.CEPT.data[GeneId %in% upper.cortex.genes,-c('GeneId')])]
tsne.blend[,cortical:=colMeans(is019.CEPT.data[GeneId %in% cortical.genes,-c('GeneId')])]
tsne.blend[,progenitor:=colMeans(is019.CEPT.data[GeneId %in% progenitor.genes,-c('GeneId')])]
tsne.blend[,other:=colMeans(is019.CEPT.data[GeneId %in% other.genes,-c('GeneId')])]


library(cowplot)

glia.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=glia)) + geom_point() + theme_bw() + scale_color_gradient(low="gray", high="red") +theme(legend.title = element_blank())

low.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=lower_cortical)) + geom_point() + theme_bw()+ scale_color_gradient(low="gray", high="red")+theme(legend.title = element_blank())

upper.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=upper_cortical)) + geom_point() + theme_bw()+ scale_color_gradient(low="gray", high="red")+theme(legend.title = element_blank())

cortical.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=cortical)) + geom_point() + theme_bw()+ scale_color_gradient(low="gray", high="red")+theme(legend.title = element_blank())

progenitor.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=progenitor)) + geom_point() + theme_bw()+ scale_color_gradient(low="gray", high="orange")+theme(legend.title = element_blank())

other.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=other)) + geom_point() + theme_bw()+ scale_color_gradient(low="gray", high="red")+theme(legend.title = element_blank())

cowplot::plot_grid(progenitor.plot, glia.plot, cortical.plot, low.plot, upper.plot, other.plot)


# multiple tSNE of cell gene means for all cells -----

load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData')
tsne.blend <- as.data.frame(is019@reductions$tsne@cell.embeddings)
tsne.blend$barcode <- row.names(tsne.blend)
tsne.blend <- as.data.table(tsne.blend)
tsne.blend
is019.data <- as.data.frame(is019@assays$RNA@data)
is019.data$GeneId <- row.names(is019.data)
is019.data <- as.data.table(is019.data)
glia.avgs <- colMeans(is019.data[GeneId %in% glia.genes,-c('GeneId')])
tsne.blend[,glia:=colMeans(is019.data[GeneId %in% glia.genes,-c('GeneId')])]
tsne.blend[,lower_cortical:=colMeans(is019.data[GeneId %in% lower.cortex.genes,-c('GeneId')])]
tsne.blend[,upper_cortical:=colMeans(is019.data[GeneId %in% upper.cortex.genes,-c('GeneId')])]
tsne.blend[,cortical:=colMeans(is019.data[GeneId %in% cortical.genes,-c('GeneId')])]
tsne.blend[,progenitor:=colMeans(is019.data[GeneId %in% progenitor.genes,-c('GeneId')])]
tsne.blend[,other:=colMeans(is019.data[GeneId %in% other.genes,-c('GeneId')])]
library(cowplot)



# re theming
new.theme <- theme(legend.title = element_blank(), axis.text = element_text(size=12, color='black'), text = element_text(size=14, color='black'), axis.line=element_line(size=0.5, color='black'), panel.background = element_blank(), axis.ticks=element_line(size=0.5), axis.ticks.length = unit(3.5, "pt"))

# glia.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=glia)) + geom_point() + scale_color_gradient(low="gray", high="red") + new.theme
# low.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=lower_cortical)) + geom_point() + scale_color_gradient(low="gray", high="red") + new.theme
# upper.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=upper_cortical)) + geom_point() + scale_color_gradient(low="gray", high="red")+ new.theme
# cortical.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=cortical)) + geom_point() + scale_color_gradient(low="gray", high="red")+ new.theme
# progenitor.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=progenitor)) + geom_point() + scale_color_gradient(low="gray", high="red")+ new.theme
# other.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=other)) + geom_point() + scale_color_gradient(low="gray", high="red")+ new.theme
# cowplot::plot_grid(progenitor.plot, glia.plot, cortical.plot, low.plot, upper.plot, other.plot)

glia.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=glia)) + geom_point(size=1) + scale_color_gradient(low="gray", high="#c200ff", limits = c(0,3), breaks = c(0, 1, 2, 3)) +new.theme + guides(color=FALSE)

low.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=lower_cortical)) + geom_point(size=1) + scale_color_gradient(low="gray", high="#b81e24", limits = c(0,3), breaks = c(0, 1, 2, 3))+new.theme + guides(color=FALSE)

upper.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=upper_cortical)) + geom_point(size=1) + scale_color_gradient(low="gray", high="#086500", limits = c(0,3), breaks = c(0, 1, 2, 3))+new.theme + guides(color=FALSE)

cortical.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=cortical)) + geom_point(size=1) + scale_color_gradient(low="gray", high="#505da9", limits = c(0,3), breaks = c(0, 1, 2, 3))+new.theme + guides(color=FALSE) # it includes SOX2 neuronal epithelium

progenitor.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=progenitor)) + geom_point(size=1) + scale_color_gradient(low="gray", high="#d95f0e", limits = c(0,3), breaks = c(0, 1, 2, 3))+new.theme + guides(color=FALSE)

# other.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=other)) + geom_point() + scale_color_gradient(low="gray", high="red")+new.theme

tsne.plot.is019 <- TSNEPlot(is019, pt.size=1, label=TRUE) + guides(color=FALSE)



#progenitor.plot <- ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=progenitor)) + geom_point() + scale_color_gradient(low="gray", high="red")+ new.theme
#progenitor.plot
#COWPLOT object-----
cowplot.obj <- cowplot::plot_grid(tsne.plot.is019, progenitor.plot, glia.plot, cortical.plot, low.plot, upper.plot)
cowplot.obj
save(cowplot.obj, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/tSNE/Cowplot_obj.RData')

# look for specific brain regions expression in IS019-----
load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData')
genelists <- as.data.table(readxl::read_xlsx("/Volumes/ncatssctl/NGS_related/marker_sets/D'Alessio_tissue_TFs_brain_only.xlsx"))
genelists <- genelists[1:25,]

genelists <- as.data.table(stack(genelists))
names(genelists) <- c('GeneId', 'tissue')

tsne.blend <- as.data.frame(is019@reductions$tsne@cell.embeddings)
tsne.blend$barcode <- row.names(tsne.blend)
tsne.blend <- as.data.table(tsne.blend)
all(tsne.blend$barcode %in% row.names(is019@meta.data))
tsne.blend[,cluster := is019@meta.data$res0.2_labels_merge]
is019.data <- as.data.frame(is019@assays$RNA@data)
is019.data$GeneId <- row.names(is019.data)
is019.data <- as.data.table(is019.data)
#tsne.blend[,glia:=colMeans(is019.data[GeneId %in% glia.genes,-c('GeneId')])]

#test.table <- data.table()

for (tissue in unique(genelists$tissue)){
  
  column_name <- quote(tissue)
  
  list_current <- unique(unlist(genelists[tissue==tissue,GeneId], use.names=F))
  
  data.relevant <- is019.data[GeneId %in% list_current,]
  
  tsne.blend[,eval(column_name):= colMeans(data.relevant[,-c('GeneId')])]

}
# check table

melt(tsne.blend, id.vars = c('barcode'))

ggplot(tsne.blend)


# one plot
ggplot(data=tsne.blend, aes(x=tSNE_1, y=tSNE_2, color=`Schwann cell`)) + geom_point() + theme_bw() + scale_color_gradient(low="gray", high="red")

# convert to heatmappable table

heat.tbl <- tsne.blend[,c(5:46)]
heat.tbl <- as.data.table(t(heat.tbl))
dim(heat.tbl)

dim(tsne.blend[,c(5:46)])

names(heat.tbl) <- tsne.blend$barcode
heat.tbl[1:5,1:5]


#

DoHeatmap(is019, features = c(unique(genelists$GeneId)),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


# too many genes..

# MEGENA: multiscale clustering of geometrical network----

install.packages('MEGENA')
library(MEGENA)
install.packages('BRETIGEA')
library(BRETIGEA)

str(aba_marker_expression,list.len = 5)
str(aba_pheno_data, list.len = 5)

ct_res = brainCells(aba_marker_expression, nMarker = 50)
table(head(ct_res))

aba_marker_expression


is019.data.df <- as.data.frame(is019.data)
row.names(is019.data.df) <- is019.data.df$GeneId
is019.data.df <- subset(is019.data.df, select=c(1:8950))

is019.data.df.brain <- subset(is019.data.df, row.names(is019.data.df) %in% row.names(aba_marker_expression))

str(is019.data.df.brain, list.len=10)

names(is019.data.df.brain) <- gsub('-',"_", names(is019.data.df.brain))

is019.data.df.brain <- as.data.frame(is019.data.df.brain)

ct_res = brainCells(is019.data.df.brain, nMarker = 50, scale=F, species="human", method="PCA")
colnames(ct_res)

#astrocytes, endothelial cells, microglia, neurons, oligodendrocytes, and OPCs

ct_res_dt <- as.data.frame(ct_res)
ct_res_dt$barcode <- row.names(ct_res_dt)
ct_res_dt <- as.data.table(ct_res_dt)

names(ct_res_dt)[1:6] <- c("astrocyte", "endothelial", "microglia", "neuron", "oligodendrocyte", "oligodendrocyte_precursor")

ct_res_dt

#cor_mic = cor.test(ct_res[, "mic"], as.numeric(aba_pheno_data$ihc_iba1_ffpe),
#                   method = "spearman")
#print(cor_mic)

ct_res_dt[,tSNE_1 := tsne.blend$tSNE_1]
ct_res_dt[,tSNE_2 := tsne.blend$tSNE_2]

astro.plot <- ggplot(data=ct_res_dt, aes(x=tSNE_1, y=tSNE_2, color=astrocyte)) + geom_point() + theme_bw() + scale_color_gradient(low="gray", high="red")

endo.plot <- ggplot(data=ct_res_dt, aes(x=tSNE_1, y=tSNE_2, color=endothelial)) + geom_point() + theme_bw() + scale_color_gradient(low="gray", high="red")

microglia.plot <- ggplot(data=ct_res_dt, aes(x=tSNE_1, y=tSNE_2, color=microglia)) + geom_point() + theme_bw() + scale_color_gradient(low="gray", high="red")

neuron.plot <- ggplot(data=ct_res_dt, aes(x=tSNE_1, y=tSNE_2, color=neuron)) + geom_point() + theme_bw() + scale_color_gradient(low="gray", high="red")

oligodendrocyte.plot <- ggplot(data=ct_res_dt, aes(x=tSNE_1, y=tSNE_2, color=oligodendrocyte)) + geom_point() + theme_bw() + scale_color_gradient(low="gray", high="red")

oligodendrocyte.prec.plot <- ggplot(data=ct_res_dt, aes(x=tSNE_1, y=tSNE_2, color=oligodendrocyte_precursor)) + geom_point() + theme_bw() + scale_color_gradient(low="gray", high="red")


cowplot::plot_grid(astro.plot, endo.plot, microglia.plot, neuron.plot, oligodendrocyte.plot, oligodendrocyte.prec.plot)


# try again to make heatmap----

heat.tbl <- ct_res_dt[,c(1:6)]
heat.tbl <- as.data.frame(t(heat.tbl))
dim(heat.tbl)

names(heat.tbl) <- ct_res_dt$barcode
heat.tbl[1:5,1:5]

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

cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red


all(ct_res_dt$barcode %in% row.names(is019@meta.data))

meta.data <- is019@meta.data
meta.data$barcode <- row.names(meta.data)
meta.data <- as.data.table(meta.data)
meta.data
all(ct_res_dt$barcode %in% meta.data$barcode)
meta.data <- meta.data[,c('barcode','res0.2_labels_merge')]
meta.data$barcode <- gsub('-','_',meta.data$barcode)
all(ct_res_dt$barcode %in% meta.data$barcode)

ct_res_dt$cluster <- meta.data$res0.2_labels_merge

Heatmap(as.matrix(heat.tbl), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = F,
        cluster_rows = TRUE, cluster_columns = TRUE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Scaled enrichment'))



h2 <- Heatmap(as.matrix(heat.tbl), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = F,
              cluster_rows = TRUE, cluster_columns = TRUE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Cell type\nproportion'),
              top_annotation = ha)

h2

mycolanno <- as.data.frame(ct_res_dt[,c('cluster')])

ha <- HeatmapAnnotation(df = mycolanno, which='column', col=list(cluster=c('CEPT3'='#F8766D', 'CEPT0_2'='#D39200', 'CEPT_4-7'='#93AA00', 'Y6'='#00BA38', 'Y5'='#00C19F', 'Y1'='#00B9E3', 'Y4'='#619CFF', 'Y7'='#DB72FB', 'Y0_3'='#FF61C3')))

draw(h2 + ha)

# violin plot format?----

heat.tbl[1:5,1:5]
ct_res_dt[1:5,]

library(cowplot)

ct_res_dt[,cluster_refactored := factor(cluster, levels = c('CEPT3', 'CEPT0_2','CEPT_4-7','Y6','Y5','Y1','Y4','Y7','Y0_3'))]


ggplot(data=ct_res_dt) +geom_jitter(aes(x=cluster, y=astrocyte), width=0.1)+ geom_violin(aes(x=cluster, y=astrocyte))
ggplot(data=ct_res_dt) +geom_violin(aes(x=cluster_refactored, y=neuron, fill=cluster_refactored))+geom_jitter(aes(x=cluster_refactored, y=neuron), width=0.05) + theme_cowplot(12) + guides(fill=FALSE) + labs(x='Cluster', y='Proportion',title='Neuron relative cell type proportion in IS019 Y/CEPT')


ct_res_dt[,sample := tstrsplit(barcode, '\\.')[2]]

ggplot(data=ct_res_dt) +geom_violin(aes(x=sample, y=neuron, fill=sample))+geom_jitter(aes(x=sample, y=neuron), width=0.05) + theme_cowplot(12) + guides(fill=FALSE) + labs(x='Sample', y='Proportion',title='Neuron relative cell type proportion in IS019 Y/CEPT')


# plot all cell scores in violins----

ct_res_dt_melt <- melt(ct_res_dt, id.vars = c('sample', 'barcode'), measure.vars = c('astrocyte','endothelial','microglia','neuron','oligodendrocyte','oligodendrocyte_precursor'))

ct_res_dt_melt[,sample_reorder:=factor(sample, levels=c('Y','CEPT'))]

ggplot(data=ct_res_dt_melt) +geom_boxplot(aes(x=sample_reorder, y=value, fill=variable))+ theme_cowplot(12) + guides(fill=FALSE) + labs(x='Sample', y='Proportion',title='Brain relative cell type proportion in IS019 Y/CEPT') + facet_wrap(nrow=3, facets = vars(variable))

fwrite(ct_res_dt_melt, 'BRETIGEA_cell_type_scoring_data.csv')

# +geom_jitter(aes(x=sample, y=value), width=0.01, alpha=0.1) 

#


# plot select genes for cortical layers -------
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/')
load('../IS019_CEPT_Y_only_res0.2merged.Seurat.RData')
vlnplot <- VlnPlot(is019, features = c('CALB2','RELN','CUX1','BCL11B','SATB2','TBR1','EOMES','PAX6','NEUROG2', 'SOX2')) + scale_fill_manual()

str(vlnplot)

#ventricular zone (sox2), subventricular zone(TBR2), lower layer neuron (CTIP2, TBR1), upper layer neuron (SATB2)

#


# RE-DO seurat object without the ribo/mt/pseudo genes.-----
unique(is019@meta.data$sample)
unique(is019@meta.data$res0.2_labels_merge)

is019.raw <- as.data.frame(as.matrix(is019@assays$RNA@counts))
is019.metadata <- as.data.frame(is019@meta.data)

is019.raw[1:5,1:5]
nrow(is019.raw)

genelist.filter <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/28309_gene_IDs_to_exclude_full_list.txt', header=F)

is019.raw <- subset(is019.raw, row.names(is019.raw) %nin% genelist.filter$V1)
nrow(is019.raw) #21240

is019.raw[grep('^MT-', row.names(is019.raw))]

is019.filtered <- CreateSeuratObject(counts = is019.raw, project = "IS019")
is019.filtered <- NormalizeData(is019.filtered)
is019.filtered <- FindVariableFeatures(is019.filtered)
is019.filtered <- ScaleData(is019.filtered, features=c(row.names(as.data.frame(as.matrix(is019.filtered@assays$RNA@data)))))
is019.filtered <- RunPCA(is019.filtered)
is019.filtered <- FindNeighbors(is019.filtered)
is019.filtered <- FindClusters(is019.filtered)
is019.filtered <- RunTSNE(is019.filtered)

Idents(is019.filtered) <- 'sample'
TSNEPlot(is019.filtered)
is019.filtered <- FindClusters(is019.filtered, resolution = 0.2)
TSNEPlot(is019.filtered)
is019.filtered@meta.data$cluster <- is019.metadata$res0.2_labels_merge
Idents(is019.filtered) <- 'cluster'
TSNEPlot(is019.filtered, do.label=TRUE)

Markers <- FindAllMarkers(is019.filtered, return.thresh = 0.001)

Markers <- as.data.table(Markers)
Markers <- Markers[p_val_adj < 0.001,]


clusters_up <- setorder(setDT(Markers), -avg_logFC)[, head(.SD, 10), keyby = cluster]

DoHeatmap(is019.filtered, features = c(clusters_up$gene),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

clusters_up <- setorder(setDT(Markers), -avg_logFC)[, head(.SD, 100), keyby = cluster]
fwrite(clusters_up, 'IS019_genefiltered_markers.csv')

TSNEPlot(is019.filtered)

# trying with res=0.1------

is019.filtered <- FindClusters(is019.filtered, resolution = 0.1)
TSNEPlot(is019.filtered)

all(row.names(is019.filtered@meta.data) %in% row.names(is019))
is019.filtered@meta.data$sample <- is019@meta.data$sample

is019.filtered@meta.data$cluster_0.1 <- paste(is019.filtered@meta.data$sample, is019.filtered@meta.data$RNA_snn_res.0.1, sep='_')
is019.filtered@meta.data
Idents(is019.filtered) <- 'cluster_0.1'
TSNEPlot(is019.filtered)


Markers <- FindAllMarkers(is019.filtered, return.thresh = 0.001)

Markers <- as.data.table(Markers)

fwrite(Markers, 'IS019_filtered_res0.1_markers.csv')

Markers.sig <- Markers[p_val_adj < 0.001,]


clusters_up <- setorder(setDT(Markers), -avg_logFC)[, head(.SD, 10), keyby = cluster]

DoHeatmap(is019.filtered, features = c(clusters_up$gene),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


Markers[GeneId=='TTR',]

genes.filter[grep('TTR', genes.filter)]


VlnPlot(is019.filtered, features='TTR')

meta.new <- is019.filtered@meta.data
meta.new$cluster_0.1_merged <- meta.new$cluster_0.1
meta.new$cluster_0.1_merged <- gsub('CEPT_1','CEPT_new',meta.new$cluster_0.1_merged)
meta.new$cluster_0.1_merged <- gsub('CEPT_2','CEPT_new',meta.new$cluster_0.1_merged)
meta.new$cluster_0.1_merged <- gsub('CEPT_4','CEPT_new',meta.new$cluster_0.1_merged)
meta.new$cluster_0.1_merged <- gsub('CEPT_new','CEPT_1-4',meta.new$cluster_0.1_merged)

Idents(is019.filtered)
is019.filtered@meta.data$cluster_0.1_merged <- meta.new$cluster_0.1_merged
Idents(is019.filtered) <- 'cluster_0.1_merged'


save(is019.filtered, file='IS019.filtered.RData')

# combine enrichr output from each cluster-----
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Enrichr/10-17 Enrichr runs/Gene_filtered_markers/')

table <- data.table('Term'=numeric(0), 'Overlap'=numeric(0), 'P-value'=numeric(0), 'Adjusted P-value'=numeric(0), 'Old P-value'=numeric(0), 'Old Adjusted P-value'=numeric(0), 'Odds Ratio'=numeric(0), 'Combined Score'=numeric(0), 'Genes'=numeric(0), 'cluster'=numeric(0), 'library'=numeric(0))

for (i in c('c0','c1-2-4','c3','y1','y2','y3','y4')){
  files <- Sys.glob('*txt')
  current_set <- files[grep(paste0('^',i), files)]
  
  for (file in current_set){
    table.temp <- fread(file)
    table.temp[,cluster := i]
    table.temp[,library := str_split_fixed(str_split_fixed(file,'_table.txt', Inf)[1], paste0(i, '_'),Inf)[2]]
    table <- rbind(table, table.temp)
    table <- unique(table)
  }
}

table[,c('P-value','Old P-value', 'Old Adjusted P-value','Odds Ratio'):=NULL]
table <- table[`Adjusted P-value` <= 0.001,]
#fwrite(table, 'Combined_enrichr_IS019.csv')

table.subset <- setorder(setDT(table), -`Combined Score`)[, head(.SD, 3), keyby = c('cluster','library')]


# plot average of genes in terms -----
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/')
load('IS019.filtered.RData')

is019 <- is019.filtered
rm(is019.filtered)
TSNEPlot(is019)

ClusterAvgs <- AverageExpression(is019) # average per cluster for every gene.

# heatmap of average per cell for genes per term.------

#set.merged <- data.table('barcode'=numeric(0), 'avg.exp'=numeric(0), 'library'=numeric(0), 'Term'=numeric(0), 'cluster'=numeric(0))

set.merged <- data.table()

for (row in 1:nrow(table.subset)){
  current.genes <- as.vector(table.subset[row,str_split_fixed(Genes,';',Inf)])
  set <- as.data.frame(colMeans(x = as.data.frame(as.matrix(is019@assays$RNA@scale.data[current.genes, ])), na.rm = TRUE))
  set <- data.table('barcode'=row.names(set), 'avg.exp'=unlist(set, use.names=F))
  #set.merged <- rbind(set.merged, set)
  
  set <- as.data.table(t(set[,c('barcode', 'avg.exp')]))
  names(set) <- unlist(set[1,], use.names=F)
  set <- set[-1,]
  set[,Term := table.subset[row,Term]]
  set[,cluster := table.subset[row,cluster]]
  set[,library:= table.subset[row,library]]
  
  set <- set[,c(8951:8953, 1:8950)]
  
  if(nrow(set.merged) >0){
    set.merged <- rbind(set.merged, set)
  }
  if(nrow(set.merged) ==0){set.merged <- set}
  
  
}

fwrite(set.merged, 'Cell_avgs_top_terms.csv')

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

cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(50) # reversed RdBu, creates Blue-white-red



mat <- set.merged
mat[1:5,1:5]
mat[,replicate:= cluster]
mat[, id := seq_len(.N), by = replicate]
mat[, id := rowid(replicate)]

mat[,replicate:= paste0(replicate,'_',id)]
mat[,id:=NULL]
mat[,c('replicate')]

mat <- as.data.frame(mat)
row.names(mat) <- mat$replicate
mat <- mat[-c(1:3,8954)]
mat[1:5,1:5]

#scaled_mat <- t(scale(t(as.numeric(as.matrix(mat)))))

#rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

col_fun <- circlize::colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

h1 <- Heatmap(as.matrix(mat), row_names_side = "right",
              column_names_side = "top", col=col_fun, show_column_names = F,
              cluster_rows = FALSE, cluster_columns = TRUE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

h1 # all looks the same color.



# try with cluster averages per term.-----

cluster.set.merged <- data.table()

for (row in 1:nrow(table.subset)){
  current.genes <- as.vector(table.subset[row,str_split_fixed(Genes,';',Inf)])
  set <- as.data.frame(ClusterAvgs <- AverageExpression(is019, features=c(current.genes)))
  names(set) <- gsub('RNA.','',names(set))
  set <- data.table(set)
  set[,gene:=current.genes]
  set[,Term := table.subset[row,Term]]
  set[,cluster := table.subset[row,cluster]]
  set[,library:= table.subset[row,library]]
  
  if(nrow(cluster.set.merged) >0){
    cluster.set.merged <- rbind(cluster.set.merged, set)
  }
  if(nrow(cluster.set.merged) ==0){cluster.set.merged <- set}
  
  
}

cluster.set.merged[,term_cluster:= paste(cluster, gsub(' ','_',Term), sep='_')]

colmeans <- cluster.set.merged[, lapply(.SD, mean),.SDcols=names(cluster.set.merged)[1:7], by=term_cluster]

fwrite(colmeans, 'Cluster_avgs_top_terms.csv')
colmeans <- fread('Cluster_avgs_top_terms.csv')

mat2 <- as.matrix(as.data.frame(colmeans)[-1])
row.names(mat2) <- colmeans$term_cluster
mat2 <- as.matrix(mat2)

library(scales)

scaled_mat <- t(scale(t(as.numeric(mat2))))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

h2 <- Heatmap(as.matrix(mat2), row_names_side = "right",
              column_names_side = "top", col=col_fun, show_column_names = T,
              cluster_rows = T, cluster_columns = T,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Avg. expr.'))

h2

cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(1000)
scaled_mat <- t(scale(t(mat2)))
rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top", col=cols.use, show_column_names = T,
        cluster_rows = T, cluster_columns = T,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-score\navg. expr.'))



#look at just one library at a time-----
table.subset.ARCHS4 <- table[library=='ARCHS4_Tissues',]
table.subset.ARCHS4 <- setorder(setDT(table.subset.ARCHS4), -`Combined Score`)[, head(.SD, 3), keyby = c('cluster','library')]
View(table.subset.ARCHS4)

cluster.set.merged.ARCHS4 <- data.table()

for (row in 1:nrow(table.subset.ARCHS4)){
  current.genes <- as.vector(table.subset.ARCHS4[row,str_split_fixed(Genes,';',Inf)])
  set <- as.data.frame(ClusterAvgs <- AverageExpression(is019, features=c(current.genes)))
  names(set) <- gsub('RNA.','',names(set))
  set <- data.table(set)
  set[,gene:=current.genes]
  set[,Term := table.subset.ARCHS4[row,Term]]
  set[,cluster := table.subset.ARCHS4[row,cluster]]
  set[,library:= table.subset.ARCHS4[row,library]]
  
  if(nrow(cluster.set.merged.ARCHS4) >0){
    cluster.set.merged.ARCHS4 <- rbind(cluster.set.merged.ARCHS4, set)
  }
  if(nrow(cluster.set.merged.ARCHS4) ==0){cluster.set.merged.ARCHS4 <- set}
}

cluster.set.merged.ARCHS4[,term_cluster:= paste(cluster, gsub(' ','_',Term), sep='_')]

colmeans.ARCHS4 <- cluster.set.merged.ARCHS4[, lapply(.SD, mean),.SDcols=names(cluster.set.merged.ARCHS4)[1:7], by=term_cluster]

fwrite(colmeans.ARCHS4, 'Cluster_avgs_top_terms_ARCHS4.csv')
colmeans.ARCHS4 <- fread('Cluster_avgs_top_terms_ARCHS4.csv')

mat3 <- as.matrix(as.data.frame(colmeans.ARCHS4)[-1])
row.names(mat3) <- colmeans.ARCHS4$term_cluster
mat3 <- as.matrix(mat3)

library(scales)

h3 <- Heatmap(as.matrix(mat3), row_names_side = "right",
              column_names_side = "top", col=col_fun, show_column_names = T,
              cluster_rows = T, cluster_columns = T,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Avg. expr.'))

h3

# just GO library------

table.subset.GO <- table[library=='GO_Biological_Process_2018',]
table.subset.GO <- setorder(setDT(table.subset.GO), -`Combined Score`)[, head(.SD, 3), keyby = c('cluster','library')]
View(table.subset.GO)

cluster.set.merged.GO <- data.table()

for (row in 1:nrow(table.subset.GO)){
  current.genes <- as.vector(table.subset.GO[row,str_split_fixed(Genes,';',Inf)])
  set <- as.data.frame(ClusterAvgs <- AverageExpression(is019, features=c(current.genes)))
  names(set) <- gsub('RNA.','',names(set))
  set <- data.table(set)
  set[,gene:=current.genes]
  set[,Term := table.subset.GO[row,Term]]
  set[,cluster := table.subset.GO[row,cluster]]
  set[,library:= table.subset.GO[row,library]]
  
  if(nrow(cluster.set.merged.GO) >0){
    cluster.set.merged.GO <- rbind(cluster.set.merged.GO, set)
  }
  if(nrow(cluster.set.merged.GO) ==0){cluster.set.merged.GO <- set}
}

cluster.set.merged.GO[,term_cluster:= paste(cluster, gsub(' ','_',Term), sep='_')]

colmeans.GO <- cluster.set.merged.GO[, lapply(.SD, mean),.SDcols=names(cluster.set.merged.GO)[1:7], by=term_cluster]

fwrite(colmeans.GO, 'Cluster_avgs_top_terms_GO.csv')
colmeans.GO <- fread('Cluster_avgs_top_terms_GO.csv')

mat4 <- as.matrix(as.data.frame(colmeans.GO)[-1])
row.names(mat4) <- colmeans.GO$term_cluster
mat4 <- as.matrix(mat4)

library(scales)

h4 <- Heatmap(as.matrix(mat4), row_names_side = "right",
              column_names_side = "top", col=col_fun, show_column_names = T,
              cluster_rows = T, cluster_columns = T,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Avg. expr.'))

h4

#WikiPathways_2019_Human------

table.subset.Pathways <- table[library=='WikiPathways_2019_Human',]
table.subset.Pathways <- setorder(setDT(table.subset.Pathways), -`Combined Score`)[, head(.SD, 3), keyby = c('cluster','library')]
View(table.subset.Pathways)

cluster.set.merged.Pathways <- data.table()

for (row in 1:nrow(table.subset.Pathways)){
  current.genes <- as.vector(table.subset.Pathways[row,str_split_fixed(Genes,';',Inf)])
  set <- as.data.frame(ClusterAvgs <- AverageExpression(is019, features=c(current.genes)))
  names(set) <- gsub('RNA.','',names(set))
  set <- data.table(set)
  set[,gene:=current.genes]
  set[,Term := table.subset.Pathways[row,Term]]
  set[,cluster := table.subset.Pathways[row,cluster]]
  set[,library:= table.subset.Pathways[row,library]]
  
  if(nrow(cluster.set.merged.Pathways) >0){
    cluster.set.merged.Pathways <- rbind(cluster.set.merged.Pathways, set)
  }
  if(nrow(cluster.set.merged.Pathways) ==0){cluster.set.merged.Pathways <- set}
}

cluster.set.merged.Pathways[,term_cluster:= paste(cluster, gsub(' ','_',Term), sep='_')]

colmeans.Pathways <- cluster.set.merged.Pathways[, lapply(.SD, mean),.SDcols=names(cluster.set.merged.Pathways)[1:7], by=term_cluster]

fwrite(colmeans.Pathways, 'Cluster_avgs_top_terms_WikiPathways.csv')
colmeans.Pathways <- fread('Cluster_avgs_top_terms_WikiPathways.csv')

mat5 <- as.matrix(as.data.frame(colmeans.Pathways)[-1])
row.names(mat5) <- colmeans.Pathways$term_cluster
mat5 <- as.matrix(mat5)

library(scales)

h5 <- Heatmap(as.matrix(mat5), row_names_side = "right",
              column_names_side = "top", col=col_fun, show_column_names = T,
              cluster_rows = T, cluster_columns = T,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Avg. expr.'))

h5

cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(1000)

scaled_mat <- t(scale(t(mat5)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

h5.recolor <- Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top", col=cols.use, show_column_names = T,
        cluster_rows = T, cluster_columns = T,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Avg. expr.'))
h5.recolor
#
# going back to IS019 res0.2 Y/CEPT original clustering with no gene filtering------
# make scaled GO term plot
load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData')

setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Enrichr/9-25 Enrichr runs/')

table <- data.table('Term'=numeric(0), 'Overlap'=numeric(0), 'P-value'=numeric(0), 'Adjusted P-value'=numeric(0), 'Old P-value'=numeric(0), 'Old Adjusted P-value'=numeric(0), 'Odds Ratio'=numeric(0), 'Combined Score'=numeric(0), 'Genes'=numeric(0), 'cluster'=numeric(0), 'library'=numeric(0))

for (i in c('CEPT3','CEPT0_2','CEPT4-7','Y0_3','Y1','Y4','Y5','Y6','Y7')){
  files <- Sys.glob('*txt')
  current_set <- files[grep(paste0('^',i), files)]
  
  for (file in current_set){
    table.temp <- fread(file)
    table.temp[,cluster := i]
    table.temp[,library := str_split_fixed(str_split_fixed(file,'_table.txt', Inf)[1], paste0(i, '.'),Inf)[2]]
    table <- rbind(table, table.temp)
    table <- unique(table)
  }
}

table[,c('P-value','Old P-value', 'Old Adjusted P-value','Odds Ratio'):=NULL]
table <- table[`Adjusted P-value` <= 0.001,]
fwrite(table, 'Combined_enrichr_IS019_oldclusters.csv')
table <- fread('Combined_enrichr_IS019_oldclusters.csv')
table.subset <- setorder(setDT(table), -`Combined Score`)[, head(.SD, 3), keyby = c('cluster','library')]

cluster.set.merged <- data.table()

for (row in 1:nrow(table.subset)){
  current.genes <- as.vector(table.subset[row,str_split_fixed(Genes,';',Inf)])
  set <- as.data.frame(ClusterAvgs <- AverageExpression(is019, features=c(current.genes)))
  names(set) <- gsub('RNA.','',names(set))
  set <- data.table(set)
  set[,gene:=current.genes]
  set[,Term := table.subset[row,Term]]
  set[,cluster := table.subset[row,cluster]]
  set[,library:= table.subset[row,library]]
  
  if(nrow(cluster.set.merged) >0){
    cluster.set.merged <- rbind(cluster.set.merged, set)
  }
  if(nrow(cluster.set.merged) ==0){cluster.set.merged <- set}
  
  
}

cluster.set.merged[,term_cluster:= paste(cluster, gsub(' ','_',Term), sep='_')]

colmeans <- cluster.set.merged[, lapply(.SD, mean),.SDcols=names(cluster.set.merged)[1:9], by=term_cluster]

fwrite(colmeans, 'Cluster_avgs_top_terms_oldclusters.csv')
#colmeans <- fread('Cluster_avgs_top_terms_oldclusters.csv')

mat2 <- as.data.frame(colmeans)
mat2 <- mat2[-1]
mat2 <- data.matrix(mat2)

library(scales)

rownames <- colmeans$term_cluster

scaled_mat <- t(scale(t(mat2)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
row.names(rescaled_mat) <- rownames
cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(1000)

h7 <- Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
              column_names_side = "top", col=cols.use, show_column_names = T,
              cluster_rows = T, cluster_columns = T,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-score\nnorm. expr.'))
h7

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top", col=cols.use, show_column_names = T,
        cluster_rows = T, cluster_columns = F,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-score\nnorm. expr.'))



# cluster only Y----
colmeans <- as.data.frame(read_xlsx('IS019_Y_enrichr_terms.xlsx'))
mat2 <- colmeans[-1]

mat2 <- data.matrix(mat2)

library(scales)

rownames <- colmeans$term_cluster

scaled_mat <- t(scale(t(mat2)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
row.names(rescaled_mat) <- rownames
cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(1000)

h7 <- Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
              column_names_side = "top", col=cols.use, show_column_names = T,
              cluster_rows = T, cluster_columns = T,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-score\nnorm. expr.'))
h7


# tSNE of just Y----
is019.raw <- as.data.frame(as.matrix(is019@assays$RNA@counts))
is019.metadata <- as.data.frame(is019@meta.data)

Y.metadata <- subset(is019.metadata, is019.metadata$sample =='Y')

Y.raw <- subset(is019.raw, select=c(names(is019.raw) %in% row.names(Y.metadata)))

is019.Y <- CreateSeuratObject(counts = Y.raw, project='IS019')
is019.Y <- NormalizeData(is019.Y)
is019.Y <- FindVariableFeatures(is019.Y)
is019.Y <- ScaleData(is019.Y, features=c(row.names(as.data.frame(as.matrix(is019.Y@assays$RNA@data)))))
is019.Y <- RunPCA(is019.Y)
is019.Y <- FindNeighbors(is019.Y)
is019.Y <- FindClusters(is019.Y)
is019.Y <- RunTSNE(is019.Y)

TSNEPlot(is019.Y)
save(is019.Y, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/IS019.Y.Seurat.RData')

is019.metadata$barcode <- row.names(is019.metadata)
Y.metadata$barcode <- row.names(Y.metadata)

Y.metadata <- merge(Y.metadata, is019.metadata, by='barcode', all.x=T, all.Y=F)

all(Y.metadata$barcode == row.names(is019.Y@meta.data))

is019.Y@meta.data$new_clusters <- Y.metadata$new_clusters
save(is019.Y, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/IS019.Y.Seurat.RData')

Idents(is019.Y) <- 'new_clusters'
TSNEPlot(is019.Y, label=T, label.size=10, pt.size=1)

FeaturePlot(is019.Y, features=c('FOXG1','OTX2','HOXA2', 'NEUROG1'), label=T, reduction = 'tsne', pt.size = 1)
VlnPlot(is019.Y, features=c('FOXG1','OTX2','HOXA2', 'NEUROG1'))
#

# tSNE of just Y clusters TBD
Y4.vs.Y5 <- FindMarkers(is019.Y, ident.1='Y4', ident.2='Y5')
Y4.vs.Y6 <- FindMarkers(is019.Y, ident.1='Y4', ident.2='Y6')
Y4.vs.Y7 <- FindMarkers(is019.Y, ident.1='Y4', ident.2='Y7')

Y.allmarkers <- FindAllMarkers(is019.Y)
Y.allmarkers <- as.data.table(Y.allmarkers)
Y.allmarkers <- Y.allmarkers[p_val_adj <= 0.001 & abs(avg_logFC) >= 1,]
Y.allmarkers
fwrite(Y.allmarkers, 'IS019_Y_only_clustering_markers.csv')

Y.allmarkers.up <- setorder(setDT(Y.allmarkers), -avg_logFC)[avg_logFC >0 & p_val_adj <= 0.001, head(.SD, 25), keyby = cluster]
Y.allmarkers.up

DoHeatmap(is019.Y, features = c(Y.allmarkers.up$gene),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


# new tSNE with decided cell clustering----
is019@meta.data$new_clusters <- is019@meta.data$res0.2_labels_merge
is019@meta.data$new_clusters <- gsub('CEPT0_2', 'CEPT1', is019@meta.data$new_clusters)
is019@meta.data$new_clusters <- gsub('CEPT3', 'CEPT2', is019@meta.data$new_clusters)
is019@meta.data$new_clusters <- gsub('CEPT_4-7', 'CEPT3', is019@meta.data$new_clusters)
is019@meta.data$new_clusters <- gsub('Y1', 'Yfour_new', is019@meta.data$new_clusters)
is019@meta.data$new_clusters <- gsub('Y4', 'Yfour_new', is019@meta.data$new_clusters)
is019@meta.data$new_clusters <- gsub('Yfour_new', 'Y4', is019@meta.data$new_clusters)
is019@meta.data$new_clusters <- gsub('Y5', 'Ysix_new', is019@meta.data$new_clusters)
is019@meta.data$new_clusters <- gsub('Y6', 'Ysix_new', is019@meta.data$new_clusters)
is019@meta.data$new_clusters <- gsub('Ysix_new', 'Y6', is019@meta.data$new_clusters)
is019@meta.data$new_clusters <- gsub('Y0_3', 'Y5', is019@meta.data$new_clusters)

Idents(is019) <- 'new_clusters'
TSNEPlot(is019, label=T, label.size = 8,pt.size=1)
save(is019, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/IS019.Seurat.newclusters.RData')


# doing DE tests again to be sure after new cluster assignment------

all.markers <- FindAllMarkers(is019)
all.markers <- as.data.table(all.markers)
all.markers <- all.markers[abs(avg_logFC) > 1 ,]
all.markers <- all.markers[,c(7,1:6)]
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/')
fwrite(all.markers,'IS019_clustering_markers_102319.csv', sep=',', col.names = T, row.names = F, quote=T)
all.markers

clusters_up <- setorder(setDT(all.markers), -avg_logFC)[avg_logFC >0 & p_val_adj <= 0.001, head(.SD, 25), keyby = cluster]
clusters_up


DoHeatmap(is019, features = c(clusters_up$gene),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

is019@meta.data$new_clusters_refactor <- is019@meta.data$new_clusters
is019@meta.data$new_clusters_refactor <- gsub('CEPT2', '2_CEPT2',is019@meta.data$new_clusters_refactor)
is019@meta.data$new_clusters_refactor <- gsub('CEPT1', '1_CEPT1',is019@meta.data$new_clusters_refactor)
is019@meta.data$new_clusters_refactor <- gsub('CEPT3', '3_CEPT3',is019@meta.data$new_clusters_refactor)
is019@meta.data$new_clusters_refactor <- gsub('Y6', '6_Y6',is019@meta.data$new_clusters_refactor)
is019@meta.data$new_clusters_refactor <- gsub('Y4', '4_Y4',is019@meta.data$new_clusters_refactor)
is019@meta.data$new_clusters_refactor <- gsub('Y5', '5_Y5',is019@meta.data$new_clusters_refactor)
is019@meta.data$new_clusters_refactor <- gsub('Y7', '7_Y7',is019@meta.data$new_clusters_refactor)

Idents(is019) <- 'new_clusters_refactor'
is019@meta.data$new_clusters_refactor <- factor(is019@meta.data$new_clusters_refactor, levels=c('1_CEPT1','2_CEPT2','3_CEPT3','4_Y4','5_Y5','6_Y6','7_Y7'))


DoHeatmap(is019, features = c(clusters_up$gene),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


is019.scaledata <- as.data.frame(as.matrix(is019@assays$RNA@scale.data))

mat <- subset(is019.scaledata, row.names(is019.scaledata) %in% clusters_up$gene)
mat[1:5,1:5]




Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top", col=cols.use, show_column_names = T,
        cluster_rows = T, cluster_columns = T,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-score\nnorm. expr.'))


# feature plots quickly for specific genes----

FeaturePlot(is019, features=c('FOXG1','OTX2','HOXA2', 'NEUROG1'), label=T, reduction = 'tsne', pt.size = 1)
VlnPlot(is019, features=c('FOXG1','OTX2','HOXA2', 'NEUROG1'))

# tSNE of just CEPT0_2 (new CEPT1) and CEPT3 (new CEPT2)----
is019.raw <- as.data.frame(as.matrix(is019@assays$RNA@counts))
is019.metadata <- as.data.frame(is019@meta.data)

CEPT.metadata <- subset(is019.metadata, is019.metadata$new_clusters =='CEPT1' | is019.metadata$new_clusters =='CEPT2')

CEPT.raw <- subset(is019.raw, select=c(names(is019.raw) %in% row.names(CEPT.metadata)))

is019.CEPT <- CreateSeuratObject(counts = CEPT.raw, project='IS019')
is019.CEPT <- NormalizeData(is019.CEPT)
is019.CEPT <- FindVariableFeatures(is019.CEPT)
is019.CEPT <- ScaleData(is019.CEPT, features=c(row.names(as.data.frame(as.matrix(is019.CEPT@assays$RNA@data)))))
is019.CEPT <- RunPCA(is019.CEPT)
is019.CEPT <- FindNeighbors(is019.CEPT)
is019.CEPT <- FindClusters(is019.CEPT)
is019.CEPT <- RunTSNE(is019.CEPT)
is019.CEPT <- RunUMAP(is019.CEPT, verbose=F)

TSNEPlot(is019.CEPT)

all(row.names(is019.CEPT@meta.data)== row.names(CEPT.metadata))
is019.CEPT@meta.data$new_clusters <- CEPT.metadata$new_clusters
Idents(is019.CEPT) <- 'new_clusters'
TSNEPlot(is019.CEPT, label=TRUE, label.size=10)

save(is019.CEPT, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/IS019.CEPT.Seurat.RData')

FeaturePlot(is019.CEPT, features=c('FOXG1','OTX2','HOXA2', 'NEUROG1'), label=T, reduction = 'tsne', pt.size = 1)
VlnPlot(is019.CEPT, features=c('FOXG1','OTX2','HOXA2', 'NEUROG1'))

# heatmap of specific genes for Seungmi-----
genes.Seungmi <- c('NEUROD6', 'NEUROD2', 'SOX11', 'SOX4', 'DCX', 'NFIB', 'BCL11A', 'STMN2', 'NELL2', 'CRMP1', 'FOXG1', 'GRIA2', 'SFRP1', 'SOX2', 'C1orf61', 'FABP7', 'HMGB2', 'FABP5', 'SOX3', 'EMX2', 'HES6', 'VIM', 'DOK5', 'DLK1', 'HILPDA', 'GDF15', 'HSPA5', 'SLC2A1', 'ADM', 'TIMP1', 'SERPINH1', 'DDIT3', 'HSPB1', 'P4HB', 'TPH1', 'NEUROD1', 'RCVRN', 'SST', 'NEFL', 'TTR', 'NRL', 'CPLX3', 'HMGB2', 'NUSAP1', 'PCLAF', 'CENPF', 'CCNB1', 'CDK1', 'HMGN2', 'MAD2L1', 'CENPW', 'CCNB2', 'CKAP2', 'APOE', 'CLU', 'ID3', 'FRZB', 'PMEL', 'BST2')

#genes.Seungmi <- all.markers[gene %in% genes.Seungmi,]
#genes.Seungmi <- setorder(setDT(genes.Seungmi), -avg_logFC)[, head(.SD, 25), keyby = cluster]

DoHeatmap(is019, features = c(genes.Seungmi),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)

# integrated with Kanton data----

IS019_Kanton <- readRDS("/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_Kanton_integration/IS019_Kanton-60D-and-less_Integrated_Processed_11062019.rds")

IS019_Kanton
UMAPPlot(IS019_Kanton)

genelist <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Latest_figures/IS019_heatmap_102319_shortlist_genes.txt', header=F)

DoHeatmap(IS019_Kanton, features = c(genelist$V1),
          raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)


# specific genes expression for Seungmi -----

DoHeatmap(is019, features = c('APOA1', 'APOA2', 'APOA4', 'APOB', 'APOC1', 'APOC2', 'APOC3', 'APOD', 'APOE', 'APOH', 'APOL1','APOL2', 'APOL3','APOL4','APOL5','APOL6','L1CAM', 'NCAM1', 'NRCAM', 'ICAM1', 'CDH2', 'CDH1', 'TJP1', 'NCAN', 'ITGB1', 'ITGA6', 'PXN', 'TTR', 'CST3', 'KRT7', 'RAB7A', 'TIMP1')
, raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)



# add updated cluster name column for seungmi----
metadata <- as.data.frame(is019@meta.data)
metadata$Cluster_new <- metadata$res0.2_labels_merge
metadata$Cluster_new <- gsub('CEPT0_2','CEPT1', metadata$Cluster_new)
metadata$Cluster_new <- gsub('CEPT3','CEPT2', metadata$Cluster_new)
metadata$Cluster_new <- gsub('CEPT_4-7','CEPT3', metadata$Cluster_new)
metadata$Cluster_new <- gsub('Y1','Y4', metadata$Cluster_new)
metadata$Cluster_new <- gsub('Y5','Y6', metadata$Cluster_new)
metadata$Cluster_new <- gsub('Y0_3','Y5', metadata$Cluster_new)

is019@meta.data$clusters <- factor(metadata$Cluster_new, levels=c('CEPT1', 'CEPT2', 'CEPT3', 'Y4','Y5','Y6','Y7'))
save(is019, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData')


# pseudotime of IS019 with Kanton data included-------
is019_kanton <- readRDS('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_Kanton_integration/IS019_Kanton_H9_0-128d_Integrated_11262019.rds')
Idents(is019_kanton)

is019_kanton <- FindVariableFeatures(is019_kanton)

library(dyno)
library(dynutils)
library(tidyverse)
library(data.table)
library(ggplot2)


mycolanno <- as.data.frame(is019_kanton@meta.data)
mycolanno <- subset(mycolanno, select=c('cellLabels'))
head(mycolanno)
names(mycolanno) <- 'cell_info'
mycolanno$cell_ids <- row.names(mycolanno)
mycolanno <- subset(mycolanno, select=c('cell_ids', 'cell_info'))

grouping <- mycolanno
head(grouping)
grouping <- as.vector(grouping)
grouping <- setNames(as.character(grouping$cell_info), row.names(grouping))
head(grouping)

subsample.genes <- is019_kanton@assays$RNA@var.features
is019_kanton.subset.raw <- as.matrix(is019_kanton@assays$RNA@counts)[subsample.genes,]
is019_kanton.subset.data <- as.matrix(is019_kanton@assays$RNA@data)[subsample.genes,]


task <- wrap_expression(
  counts = t(is019_kanton.subset.raw),
  expression= t(is019_kanton.subset.data),
  cell_info = mycolanno
)

model.embeddr <- infer_trajectory(task, 'embeddr') 
 # not enough memory

# library(scater)
# library(embeddr)
# pd <- new('AnnotatedDataFrame', data=mycolanno)
# sce <- SummarizedExperiment(assays=SimpleList(exprs=is019_kanton.subset.raw), colData = mycolanno)
# sce <- embeddr(sce)

#try with phenopath-----
set.seed(123L)
library(Seurat)
library(phenopath)
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
is019_kanton <- readRDS('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_Kanton_integration/IS019_Kanton_H9_0-128d_Integrated_11262019.rds')
is019_kanton <- FindVariableFeatures(is019_kanton)
subsample.genes <- is019_kanton@assays$RNA@var.features
is019_kanton.subset.data <- as.matrix(is019_kanton@assays$RNA@data)[subsample.genes,]
is019_kanton.subset.data.offset <- is019_kanton.subset.data + 0.000001

mycolanno <- as.data.frame(is019_kanton@meta.data)
mycolanno <- subset(mycolanno, select=c('cellLabels'))
names(mycolanno) <- 'cell_info'
head(mycolanno)

sce <- SummarizedExperiment(assays=list(exprs=is019_kanton.subset.data.offset), colData = mycolanno)
sce

fit <- phenopath(sce, ~ cell_info) 
save(fit, file='/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_Kanton_integration/pseudotime/Phenopath_fit.RData')
#https://bioconductor.org/packages/release/bioc/vignettes/phenopath/inst/doc/introduction_to_phenopath.html#using-an-summarizedexperiment-as-input
#https://kieranrcampbell.github.io/phenopath/phenopath_shalek_vignette.html
#
plot_elbo(fit)

zdf <- data_frame(z = trajectory(fit), cluster = sce$cell_info)


ggplot(zdf, aes(x = cluster, y = log(z+1), fill = cluster)) + 
  geom_violin(alpha = 0.8) +
  theme(legend.position = "none") + 
  scale_fill_brewer(palette = "Set1") +
  xlab("Cluster") +
  ylab("Pathway score\n(pseudotime)") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)) 


ggplot(zdf, aes(x = cluster, y = log(z+1), fill = cluster)) + 
  geom_boxplot(alpha = 0.8) +
  theme(legend.position = "none") +
  xlab("Cluster") +
  ylab("Pathway score\n(pseudotime)") +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)) 

qplot(zdf$z, trajectory(fit)) +
  xlab("True z") + ylab("Phenopath z")

interaction_effects(fit)
gene_names <- paste0("gene", seq_len(ncol(fit$m_beta)))

df_beta <- data_frame(beta = interaction_effects(fit),
                      beta_sd = interaction_sds(fit),
                      is_sig = significant_interactions(fit),
                      gene = gene_names)

#ints <- interactions(fit)
#ints$feature <- fData(sce_hvg)$mgi_symbol


# try plotting UMAP/TSNE with pseudotime values----
is019_kanton@meta.data$phenopath <- trajectory(fit)

umap_data <- as.data.frame(is019_kanton@reductions$umap@cell.embeddings)

umap_data$barcode <- row.names(is019_kanton@meta.data)

all(row.names(mycolanno) == umap_data$barcode)

umap_data$cluster <- mycolanno$cell_info
umap_data$phenopath <- as.numeric(trajectory(fit))

umap_data <- subset(umap_data, umap_data$cluster != "Y7")

ggplot(data=umap_data, aes(x=UMAP_1, y=UMAP_2, color=phenopath))+geom_point()

ggplot(data=umap_data, aes(x=cluster, y=phenopath))+geom_jitter() + coord_flip()


# plot specific markers from other researcher---
load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData')
is019

markers <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/markers_to_check.txt', header=F)

VlnPlot(is019, features=c(markers$V1))

# updated violin plots, may 20, 2020-----
is019@meta.data$sample <- factor(is019@meta.data$sample, levels=c('Y', 'CEPT'))
Idents(is019)<- 'sample'
is019@active.ident

is019.data <- as.data.frame(is019@assays$RNA@data)
is019.data <- subset(is019.data, row.names(is019.data) %in% c('SOX2', 'PAX6', 'EOMES', 'TBR1', 'BCL11B', 'CUX1', 'SATB2', 'CALB2', 'RELN'))
is019.data[1:2,1:2]
is019.data$GeneId <- row.names(is019.data)
is019.data <- as.data.table(is019.data)
is019.data <- is019.data[,c(8951, 1:8950)]
is019.data[1:5,1:5]

melt(is019.data)[1:2,1:2]
is019.data <- melt(is019.data, id.vars='GeneId')
is019.data[1:2,1:2]
is019.data[,Treatment:= tstrsplit(variable, '\\.')[2]]
names(is019.data) <- c('GeneId','cell_id', 'expression', 'treatment')
is019.data[,gene_color:= GeneId]
is019.data$gene_color <- gsub('SOX2', '#505da9',is019.data$gene_color)
is019.data$gene_color <- gsub('PAX6', '#c200ff',is019.data$gene_color)
is019.data$gene_color <- gsub('EOMES', '#c200ff',is019.data$gene_color)
is019.data$gene_color <- gsub('TBR1', '#b81e24',is019.data$gene_color)
is019.data$gene_color <- gsub('BCL11B', '#b81e24',is019.data$gene_color)
is019.data$gene_color <- gsub('CUX1', '#71ad59',is019.data$gene_color)
is019.data$gene_color <- gsub('SATB2', '#71ad59',is019.data$gene_color)
is019.data$gene_color <- gsub('CALB2', '#88b5ed',is019.data$gene_color)
is019.data$gene_color <- gsub('RELN', '#88b5ed',is019.data$gene_color)

is019.data$GeneId <- factor(is019.data$GeneId,
                            levels=c('SOX2', 'PAX6', 'EOMES', 'TBR1', 'BCL11B', 'CUX1', 'SATB2', 'CALB2', 'RELN'))

is019.data$treatment <- factor(is019.data$treatment, levels=c('Y','CEPT'))


is019.data.nonzero <- is019.data[expression>0,]
is019.data.nonzero <- na.omit(is019.data.nonzero)

attach(is019.data.nonzero)
ggplot(is019.data.nonzero, aes(x=treatment, y=expression, fill=GeneId, color=GeneId))+
  geom_jitter(alpha=0.4)+
  geom_boxplot(width=0.5, alpha=0.8, color='black', fill='white', outlier.shape=NA)+
  facet_wrap(vars(GeneId), ncol=3)+theme_bw()+
  theme(panel.grid = element_blank(), axis.text=element_text(size=15),
        strip.background = element_rect(fill='white', colour = 'white'),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size=15),
        axis.title=element_text(size=15))+
  guides(fill=FALSE, color=FALSE, alpha=FALSE)+ylim(c(0,3.75))+
  labs(x='Treatment',y='Expression level', title='')+
  scale_fill_manual(values=c('SOX2' = '#505da9', 'PAX6' = '#c200ff', 'EOMES' = '#c200ff', 'TBR1' = '#b81e24', 'BCL11B' = '#b81e24', 'CUX1' = '#71ad59', 'SATB2' = '#71ad59', 'CALB2' = '#88b5ed', 'RELN' = '#88b5ed'))+
  scale_color_manual(values=c('SOX2' = '#505da9', 'PAX6' = '#c200ff', 'EOMES' = '#c200ff', 'TBR1' = '#b81e24', 'BCL11B' = '#b81e24', 'CUX1' = '#71ad59', 'SATB2' = '#71ad59', 'CALB2' = '#88b5ed', 'RELN' = '#88b5ed'))

fwrite(is019.data, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/VlnPlot/IS019_boxplot_data.csv')


# ion channel heatmaps using tao's marker list and IS019 ANOVA test results from Seungmi-----
setwd('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Manuscript_heatmaps/')
ion.channels <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/marker_sets/Ion channels group_updated.xlsx'))

de.cortical <- na.omit(fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Partek_DE_CEPT12_vs_Y456.txt'))
de.cortical
de.cortical$`Gene symbol` [de.cortical$`Gene symbol` %in% ion.channels$`Approved symbol`]

de.treatment <- na.omit(fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/Partek_DE_all_CEPT_vs_Y.txt'))
de.treatment$`Gene symbol`[de.treatment$`Gene symbol` %in% ion.channels$`Approved symbol`]


#partek----


is019.metadata <- as.data.frame(is019@meta.data)
is019.metadata[1:2,1:2]
is019.metadata <- subset(is019.metadata, select=c('clusters'))
is019.metadata$cell_id <- row.names(is019.metadata)
is019.metadata<- as.data.table(is019.metadata)
is019.metadata <- is019.metadata[,c('cell_id','clusters')]
fwrite(is019.metadata, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/IS019_CEPT_Y_metadata_clusters_for_Partek.csv')
is019.metadata <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/IS019_CEPT_Y_metadata_clusters_for_Partek.csv')
# is019.metadata.full <- as.data.frame(is019@meta.data)
# is019.metadata.full <- subset(is019.metadata.full, sample=='CE')
# is019.metadata.full <- subset(is019.metadata.full, select=c('sample'))
# write.table(is019.metadata.full, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/IS019_CE_metadata_clusters_for_Partek.csv')

is019.metadata[,cell_type:=clusters]
is019.metadata$cell_type <- gsub('CEPT1','cortical neuron', is019.metadata$cell_type)
is019.metadata$cell_type <- gsub('CEPT2','radial glia', is019.metadata$cell_type)
is019.metadata$cell_type <- gsub('CEPT3','non-neuronal', is019.metadata$cell_type)
is019.metadata$cell_type <- gsub('Y4','cortical neuron', is019.metadata$cell_type)
is019.metadata$cell_type <- gsub('Y5','radial glia', is019.metadata$cell_type)
is019.metadata$cell_type <- gsub('Y6','neural progenitor', is019.metadata$cell_type)
is019.metadata$cell_type <- gsub('Y7','non-neuronal', is019.metadata$cell_type)

is019.metadata[,treatment := tstrsplit(cell_id, '\\.')[2]]
fwrite(is019.metadata, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/IS019_CEPT_Y_metadata_clusters_for_Partek.csv')

#
is019.raw.counts <- as.data.frame(is019@assays$RNA@counts)
is019.raw.counts$GeneId <- row.names(is019.raw.counts)
is019.raw.counts <- as.data.table(is019.raw.counts)
is019.raw.counts <- is019.raw.counts[,c(8951, 1:8950)]
is019.raw.counts[1:5,1:5]

fwrite(is019.raw.counts, '/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Clustering/IS019_CEPT_Y_counts_for_Partek.csv')


# june 8, 2020 heatmap update with VIM moved to Y6 cluster.----
load('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/IS019_CEPT_Y_only_res0.2merged.Seurat.RData')
genelist <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS019/Analysis/Manuscript_heatmaps/IS019_unbiased_up-down_heatmap_genes_v3.txt', header=F)
DoHeatmap(is019, features = c(genelist$V1), size=6, raster=F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE)+theme(axis.text=element_text(size=6))
