library(data.table)
library(readxl)
library(stringr)

indir <- '/data/ISB013/'
outdir <- '/data/ISB013/outdir'
setwd('/data/ISB013/')

# swarm -f htseq_swarm.sh -t 4 -g 4 --module htseq --time=10:00:00

# convert gene names---------
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(stringr)
library(dendextend)
library(ComplexHeatmap)
library(gtools)
library(doParallel)
library(foreach)
library(scales)


## getting the intersection of genes for each file.
setwd('/home/inmanjm/data/ISB013/all_countfiles/')
files <- Sys.glob('*txt')

new.ex <- fread(files[1], header=F)
old.ex <- fread(files[10], header=F)

genes.intersection <- intersect(new.ex$V1, old.ex$V1)
genes.intersection <- genes.intersection[-grep('A[A-Z]+\\d+\\.\\d', genes.intersection)]
genes.intersection <- genes.intersection[-grep('^RP',genes.intersection)]
genes.intersection <- genes.intersection[-grep('^MT', genes.intersection)]
genes.intersection <- genes.intersection[-grep('^LINC', genes.intersection)]

for (countfile in files){
  data <- fread(countfile, header=F)
  data <- data[V1 %in% genes.intersection,]
  #print(nrow(data))
  outfile <- gsub('_gene_symbol_counts.txt', '_intersection_gene_symbol_counts.txt', countfile)
  data <- data[order(V1)]
  fwrite(data, paste0('/home/inmanjm/data/ISB013/all_countfiles_intersection/',outfile), col.names = F, row.names = F, quote=F, sep='\t')
  message("Done with ", countfile)
  
}
##

setwd('/home/inmanjm/data/ISB013/')

sampletable <- as.data.table(readxl::read_xlsx('ISB003-13_sampletable.xlsx'))
sampletable
directory <- "/home/inmanjm/data/ISB013/all_countfiles_intersection/"
dds <- DESeqDataSetFromHTSeqCount(sampletable,
                                  directory = directory,
                                  design = ~ condition)

dds <- DESeq(dds)
#save(dds, file='ISB012.Tao.DDS.pre_correction.RData')

# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# dds <- DESeq(dds)
# save(dds, file='ISB012.DDS.mincounts.RData')

# 

# RUVSeq batch correction-----
library(RUVSeq) # http://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf
library(EDASeq)
library(RColorBrewer)

counts <- as.data.frame(counts(dds, normalized=F))
names(counts) <- sampletable$replicate
counts <- as.matrix(counts)

x <- as.factor(sampletable$condition)
set1 <- newSeqExpressionSet(counts,
                           phenoData = data.frame(x, row.names= c(sampletable$replicate) ))
set1

colors <- brewer.pal(8, "Set2")
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set1, col=colors[x], cex=1.2)


spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3

spikes <- spikes[spikes %in% row.names(counts)]

set2 <- RUVg(set1, spikes, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)


#### YOU"RE HERE, JASON

#plot before and after side by side
par(mfcol=c(1,2))
plotRLE(set1, outline=FALSE, ylim=c(-2, 2), col=colors[x])
plotRLE(set2, outline=FALSE, ylim=c(-2, 2), col=colors[x])

# set up new dds with batch correction factor included in the design
condition <- x

names(pData(set2)) <- c('condition', 'W_1')

dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)

#save(dds, file='ISB003-13.Seungmi.DDS.after_correction.RData')

View(as.data.frame(dds@colData))

# DE tests------
library(doParallel)

source('/mnt/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')

tests_table <- as.data.table(readxl::read_xlsx('/home/inmanjm/data/ISB013/teststable.xlsx', sheet=1))

testdir <- '/home/inmanjm/data/ISB013/DE'
volcanodir <- '/home/inmanjm/data/ISB013//Volcano_plots'

keep <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 3
dds <- dds[keep,]

save(dds, file='ISB013.Seungmi.DDS.after_correction.mincounts.RData')

cl <- makeCluster(10)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]
  
  testdir <- testdir
  volcanodir <- volcanodir
  
  DESeq2_pipeline(dds, condition1, condition2, testdir, volcanodir)
}

stopCluster(cl)

# new heatmap with new and old D28 samples-----
setwd('/home/inmanjm/data/ISB013/')
load('ISB013.Seungmi.DDS.after_correction.mincounts.RData')
dds

genelist <- c('POU5F1', 'NANOG', 'SOX2', 'SOX21', 'OTX2', 'LHX2', 'PAX3', 'NEUROG1', 'TFAP2A', 'TFAP2B', 'SOX10', 'POU4F1', 'ISL1', 'DRGX', 'SLC17A6', 'PRPH', 'RUNX1', 'TAC1', 'RET', 'P2RX3', 'TRPV1', 'SCN9A', 'SCN10A', 'SCN11A', 'NTRK2', 'NRG2', 'PDGFRA', 'BDNF', 'NGF', 'SOX17', 'GATA4', 'GATA6', 'TBXT', 'EOMES')

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")


dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genelist]
length(found)
length(genelist)
#genes[!genes %in% found] #

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genelist)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)

# rearranged genes.
mat <- as.data.frame(readxl::read_xlsx('mat_for_heatmap.xlsx'))
row.names(mat) <- mat$GeneId
mat <- mat[-1]
mat

#names(mat)<- gsub('_',' ',coldata$condition)

mat

cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red

scaled_mat <- t(scale(t(mat)))

h4 <- Heatmap(as.matrix(scaled_mat), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
              cluster_rows = FALSE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

h4

library("scales")

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

h5 <- Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
              cluster_rows = FALSE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))
h5

# save data for Tao to use-----
data <- as.data.frame(counts(dds, normalized=TRUE))
newnames <- row.names(dds@colData)

names(data) <- newnames
data
#fwrite(data, 'ISB012.Tao.normalized.counts.csv', sep=',', quote=F, col.names = T, row.names = T)
data <- as.data.frame(counts(dds, normalized=FALSE))
names(data) <- newnames
data
#fwrite(data, 'ISB012.Tao.raw.counts.csv', sep=',', quote=F, col.names = T, row.names = T)

# try taking out H9now_D28_new_1 because it is skewing the colors----

mat.df <- as.data.frame(mat)
mat.sub <- mat.df[,c(1:18, 20:24)]

rescaled_mat_sub <- rescale(t(base::scale(t(mat.sub))), to=c(-2,2))

Heatmap(as.matrix(rescaled_mat_sub), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))


# new heatmaps request, august 9-----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/')
load('ISB012.Tao.DDS.after_correction.RData')
genelist <- c('POU5F1', 'NANOG', 'SOX2', 'SOX21', 'OTX2', 'LHX2', 'PAX3', 'NEUROG1', 'TFAP2A', 'TFAP2B', 'SOX10', 'POU4F1', 'ISL1', 'DRGX', 'SLC17A6', 'PRPH', 'RUNX1', 'TAC1', 'RET', 'P2RX3', 'TRPV1', 'SCN9A', 'SCN10A', 'SCN11A', 'NTRK2', 'NRG2', 'PDGFRA', 'BDNF', 'NGF')

dds <- dds[ ,dds$condition %in% c('H9noc_D0', 'H9noc_D4', 'H9noc_D8', 'H9noc_D12', 'H9noc_D21', 'H9noc_D28_new') ]

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genelist]
length(found)
length(genelist)

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genelist)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)
mat$GeneId <- row.names(mat)

genelist.dt <- data.table(GeneId=genelist, order=c(1:length(genelist)))

mat <- merge(mat, genelist.dt, by='GeneId', all=T)
mat <- mat[order(mat$order),]
mat <- mat[,-c(20)]
mat <- as.data.frame(mat)
row.names(mat) <- mat$GeneId
mat <- mat[-1]
mat

scaled_mat <- t(scale(t(mat)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = TRUE, cluster_columns = TRUE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))


# before scaling

Heatmap(as.matrix(mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))



# another heatmap with all samples...-----

load('ISB012.Tao.DDS.after_correction.RData')
genelist <- c('POU5F1', 'NANOG', 'SOX2', 'SOX21', 'OTX2', 'LHX2', 'PAX3', 'NEUROG1', 'TFAP2A', 'TFAP2B', 'SOX10', 'POU4F1', 'ISL1', 'DRGX', 'SLC17A6', 'PRPH', 'RUNX1', 'TAC1', 'RET', 'P2RX3', 'TRPV1', 'SCN9A', 'SCN10A', 'SCN11A', 'NTRK2', 'NRG2', 'PDGFRA', 'BDNF', 'NGF')

dds <- dds[ ,dds$condition %in% c('H9noc_D0', 'H9noc_D4', 'H9noc_D8', 'H9noc_D12', 'H9noc_D21', 'H9noc_D28','H9noc_D28_new') ]

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genelist]
length(found)
length(genelist)

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genelist)

# 013120 making a DDS for all samples in ISB012.

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)
mat$GeneId <- row.names(mat)

genelist.dt <- data.table(GeneId=genelist, order=c(1:length(genelist)))

mat <- merge(mat, genelist.dt, by='GeneId', all=T)
mat <- mat[order(mat$order),]
mat <- mat[,-c(23)] # remove laste 'order' column
mat <- as.data.frame(mat)
row.names(mat) <- mat$GeneId
mat <- mat[-1]
mat

scaled_mat <- t(scale(t(mat)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))


# take out LHX2, PDGFRA (still no endo/meso): only new D28; 
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/')
load('ISB012.Tao.DDS.after_correction.RData')

genelist <- c('POU5F1', 'NANOG', 'SOX2', 'SOX21', 'OTX2', 'PAX3', 'NEUROG1', 'TFAP2A', 'TFAP2B', 'SOX10', 'POU4F1', 'ISL1', 'DRGX', 'SLC17A6', 'PRPH', 'RUNX1', 'TAC1', 'RET', 'P2RX3', 'TRPV1', 'SCN9A', 'SCN10A', 'SCN11A', 'NTRK2', 'NRG2', 'BDNF', 'NGF')

dds <- dds[ ,dds$condition %in% c('H9noc_D0', 'H9noc_D4', 'H9noc_D8', 'H9noc_D12', 'H9noc_D21', 'H9noc_D28_new') ]

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genelist]
length(found)
length(genelist)

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genelist)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)
mat$GeneId <- row.names(mat)

genelist.dt <- data.table(GeneId=genelist, order=c(1:length(genelist)))

mat <- merge(mat, genelist.dt, by='GeneId', all=T)
mat <- mat[order(mat$order),]
mat <- mat[,-c(20)]
mat <- as.data.frame(mat)
row.names(mat) <- mat$GeneId
mat <- mat[-1]
mat

scaled_mat <- t(scale(t(mat)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

#all samples including axol and endoderm and mesoderm but without LHX2, PDGFRA, or SOX17.

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/')
load('ISB012.Tao.DDS.after_correction.RData')

genelist <- c('POU5F1', 'NANOG', 'SOX2', 'SOX21', 'OTX2', 'PAX3', 'NEUROG1', 'TFAP2A', 'TFAP2B', 'SOX10', 'POU4F1', 'ISL1', 'DRGX', 'SLC17A6', 'PRPH', 'RUNX1', 'TAC1', 'RET', 'P2RX3', 'TRPV1', 'SCN9A', 'SCN10A', 'SCN11A', 'NTRK2', 'NRG2', 'BDNF', 'NGF', 'GATA4', 'GATA6', 'TBXT', 'EOMES')

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genelist]
length(found)
length(genelist)

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genelist)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)
mat$GeneId <- row.names(mat)

genelist.dt <- data.table(GeneId=genelist, order=c(1:length(genelist)))

mat <- merge(mat, genelist.dt, by='GeneId', all=T)
mat <- mat[order(mat$order),]
mat <- mat[,-c(26)]
mat <- as.data.frame(mat)
row.names(mat) <- mat$GeneId
mat <- mat[-1]
mat

scaled_mat <- t(scale(t(mat)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

rescaled_mat <- rescaled_mat[,c(4:24,1:3)]

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = FALSE, cluster_columns = FALSE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))

# 013120 making a DDS for all ISB012 samples----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/')

sampletable <- as.data.table(readxl::read_xlsx('ISB012_full_sampletable.xlsx'))
sampletable
directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Countfiles_gene_symbols/"
dds <- DESeqDataSetFromHTSeqCount(sampletable,
                                  directory = directory,
                                  design = ~ condition)

dds <- DESeq(dds)
save(dds, file='ISB012.DDS.RData')
