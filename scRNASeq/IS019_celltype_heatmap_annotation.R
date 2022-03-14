myrowanno <- as.data.frame(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/Genes for the temporal bulkRNA seq heatmap.xlsx', sheet=2))
row.names(myrowanno) <- myrowanno$GeneId
myrowanno <- subset(myrowanno, select=c('Geneset'))
myrowanno

ha <- HeatmapAnnotation(df = myrowanno, which='row', width=unit(1, 'cm'), col=list(Geneset=c('Pluripotency'='#a6cee3', 'Neuroepithelium'='#1f78b4', 'Radial Glial'='#b2df8a', 'Astrocyte Precursors'='#33a02c', 'Astrocytes'='#fb9a99', 'Neurons'='#e31a1c', 'Oligodendrocytes'='#fdbf6f', 'Endothelial'='#ff7f00', 'Microglia'='#cab2d6', 'Proliferation'='#6a3d9a')))



mycolanno <- as.data.frame(ct_res_dt[,c('cluster')])


p <- TSNEPlot(is019)
p
q <- ggplot_build(p)$data
q <- as.data.table(q)
q
unique(q$colour)
c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3")


ha <- HeatmapAnnotation(df = mycolanno, which='column', col=list(cluster=c('CEPT3'='#F8766D', 'CEPT0_2'='#D39200', 'CEPT_4-7'='#93AA00', 'Y6'='#00BA38', 'Y5'='#00C19F', 'Y1'='#00B9E3', 'Y4'='#619CFF', 'Y7'='#DB72FB', 'Y0_3'='#FF61C3')))
draw(ha)
ha
draw(ha, 1:8950)

