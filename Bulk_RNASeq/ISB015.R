library(data.table)
library(readxl)
library(stringr)
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/')
samples <- fread('samples.txt', header=F)

indir <- '/data/NCATS_ifx/data/mRNASeq/ISB015/'
outdir <- '/data/NCATS_ifx/data/mRNASeq/ISB015/htseq'
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/')

for (sample in samples){
  sample <- gsub('_star.genome.sorted.bam', '', sample)
  
  cat(paste0('htseq-count -f bam -r pos -s no -t exon -m union ', indir , 'bam/' , sample  , '_star.genome.sorted.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > ', outdir, '/', sample , '_htseq_counts.txt' , "\n\n"))
}

# swarm -f htseq_swarm.sh -t 4 -g 4 --module htseq --time=10:00:00

# format count files for partek------
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/Countfiles_gene_symbol/')

files <- Sys.glob('*.txt')
files



for (file in files[1]){
  data <- fread(file, header=F)
  samplename <- str_split_fixed(file, '_gene_symbol_counts.txt', n=Inf)[1]
  samplename <- str_split_fixed(samplename, '_', n=Inf)
  samplename <- paste(samplename[2:length(samplename)],collapse="_")
  names(data) <- c('GeneId', paste0(samplename))
  
  genes <- data$GeneId[-grep('[A-Z][A-Z]+\\d+\\.\\d', data$GeneId)]
  genes <- genes[-grep('^RP',genes)]
  genes <- genes[-grep('^MT', genes)]
  genes <- genes[-grep('^LINC', genes)]
  
  data <- data[GeneId %in% genes,]
}


merged_data <- data.table('GeneId'=data$GeneId)

for (file in files){
  data <- fread(file, header=F)
  samplename <- str_split_fixed(file, '_gene_symbol_counts.txt', n=Inf)[1]
  samplename <- str_split_fixed(samplename, '_', n=Inf)
  samplename <- paste(samplename[2:length(samplename)],collapse="_")
  names(data) <- c('GeneId', paste0(samplename))
  
  data <- data[GeneId %in% genes,]
  
  merged_data <- merge(merged_data, data, by='GeneId', all=T)
  
}

fwrite(merged_data, 'ISB015_merged_for_partek.csv', sep=',', col.names = T, row.names = F, quote=T)

# set up DESeq-----

merged_data <- read.csv('ISB015_merged_for_partek.csv', row.names = 1)
keep.samples <- names(merged_data)[grep('^BLBP_', names(merged_data))]
data <- subset(merged_data, select=c(keep.samples))

cts <- as.data.frame(data)
cts <- as.matrix(cts)
cts <- apply(cts, 2, as.numeric)
row.names(cts) <- row.names(data)

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/')

sampletable <- as.data.table(read_xlsx('teststable.xlsx', sheet=2))

condition <- data.table(names = sampletable$replicate)
condition[,condition := sampletable$condition]
condition <- condition$condition

coldata <- data.frame(condition = condition, row.names = names(data))
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
dds

save(dds, file='ISB015.DDS.RData')

#RUVSeq-----

library(RUVSeq) # http://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf
library(EDASeq)
library(RColorBrewer)

counts <- as.data.frame(counts(dds, normalized=F))
names(counts) <- sampletable$replicate
counts <- as.matrix(counts)

x <- as.factor(sampletable$condition)
set <- newSeqExpressionSet(counts,
                           phenoData = data.frame(x, row.names= c(sampletable$replicate) ))
set

colors <- brewer.pal(12, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3

spikes <- spikes[spikes %in% row.names(counts)]

set2 <- RUVg(set, spikes, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)

set3 <- RUVg(set, empirical, k=1)
pData(set3)


#plot before and after side by side----
par(mfcol=c(1,3))

colors <- c('#a6cee3', '#a6cee3', '#a6cee3', '#1f78b4', '#1f78b4', '#1f78b4', '#b2df8a', '#b2df8a', '#b2df8a', '#33a02c', '#33a02c', '#33a02c', '#fb9a99', '#fb9a99', '#fb9a99', '#e31a1c', '#e31a1c', '#e31a1c', '#fdbf6f', '#fdbf6f', '#fdbf6f', '#ff7f00', '#ff7f00', '#ff7f00', '#cab2d6', '#cab2d6', '#cab2d6', '#6a3d9a', '#6a3d9a', '#6a3d9a', '#ffff99', '#ffff99', '#ffff99', '#b15928', '#b15928', '#b15928')

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors)
plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors)
# set up new dds with batch correction factor included in the design
condition <- x

names(pData(set3)) <- c('condition', 'W_1')

dds <- DESeqDataSetFromMatrix(countData = counts(set3),
                              colData = pData(set3),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)

save(dds, file='ISB015.DDS.after_empirical_correction.RData')

condition <- x
names(pData(set2)) <- c('condition', 'W_1')
dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)
save(dds, file='ISB015.DDS.after_housekp_correction.RData')


# empirical genes correction DE------
load("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/ISB015.DDS.after_empirical_correction.RData")
load("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/ISB015.DDS.after_housekp_correction.RData")

library(foreach)
library(doParallel)
library(data.table)
library(readxl)

source('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')

tests_table <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/teststable.xlsx', sheet=1))

testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/Housekp'
volcanodir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/Housekp/Volcano_plots'


#optional to filter low expressed genes.
keep <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 3
dds <- dds[keep,]

cl <- makeCluster(6)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]
  
  testdir <- testdir
  volcanodir <- volcanodir
  
  DESeq2_pipeline(dds, condition1, condition2, testdir, volcanodir)
}

stopCluster(cl)


# RRHO2------
library(devtools)
install_github("RRHO2/RRHO2")
library(RRHO2)
library(lattice)


list.length <- 100
list.names <- paste('Gene',1:list.length, sep='')

data1 <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/Housekp/DE.BLBP_Lonza_d8.vs.BLBP_Lonza_d60.csv')
data1[,log2FC_manual:= log2((`BLBP_Lonza_d60`+0.0001)/(`BLBP_Lonza_d8`+0.0001) )]

data2 <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/Empirical/DE.BLBP_Lonza_d8.vs.BLBP_Lonza_d60.csv')
data2[,log2FC_manual:= log2((`BLBP_Lonza_d60`+0.0001)/(`BLBP_Lonza_d8`+0.0001) )]

#gene.list1 <- data.frame(list.names, data1$GeneId[1:100])
#gene.list2 <- data.frame(list.names, data2$GeneId[1:100])

#RRHO.example <- RRHO2(gene.list1, gene.list2, BY=TRUE, alternative='enrichment')

plotFolder <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/RRHO2'
#system(paste('mkdir -p', plotFolder))
set.seed(15213)

#gene.list1<- data.frame(list.names, sample(100)*sample(c(1,-1),100,replace=TRUE))
#gene.list2<- data.frame(list.names, sample(100)*sample(c(1,-1),100,replace=TRUE))


#data1 <- data1[abs(log2FoldChange) >=1,]
#data2 <- data1[abs(log2FoldChange) >=1,]

gene.list1<- data.frame(data1$GeneId, data1$log2FC_manual)
gene.list2<- data.frame(data2$GeneId, data2$log2FC_manual)


RRHO.example <-  RRHO2(gene.list1, gene.list2, labels=c('Housekp','Empirical'), plots=TRUE, outputdir=plotFolder, BY=TRUE, log10.ind=TRUE)
levelplot(RRHO.example$hypermat)


sig.in.housekp <- data1[abs(log2FC_manual)>=1 & padj <= 0.01,]
sig.in.emp <- data2[abs(log2FC_manual)>=1 & padj <= 0.01,]

#RRHO <-  RRHO2(data.frame(sig.in.housekp$GeneId, sig.in.housekp$log2FoldChange),
 #              data.frame(sig.in.emp$GeneId, sig.in.emp$log2FoldChange),
  #             labels=c('Housekp','Empirical'), plots=TRUE, outputdir=plotFolder, BY=TRUE, log10.ind=TRUE)


#data1[,rank:= c(1:nrow(data1) )]

#overlap <- merge(data1, data2, by=c("GeneId", "BLBP_Lonza_d8", "BLBP_Lonza_d60", "log2FC_manual" ), all=F)
#overlap

sig.in.housekp <- sig.in.housekp[order(-abs(log2FC_manual)),]
sig.in.emp <- sig.in.emp[order(-abs(log2FC_manual)),]

sig.in.housekp[,rank:= c(1:nrow(sig.in.housekp))]
sig.in.emp[,rank:= c(1:nrow(sig.in.emp))]

sig.in.housekp
sig.in.emp

fwrite(sig.in.emp, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/RRHO2/BLBP_Lonza_d8 vs d60/sig.in.emp.csv')
fwrite(sig.in.housekp, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/RRHO2/BLBP_Lonza_d8 vs d60/sig.in.housekp.csv')

gene.list1<- data.frame(sig.in.housekp$GeneId, sig.in.housekp$log2FC_manual)
gene.list2<- data.frame(sig.in.emp$GeneId, sig.in.emp$log2FC_manual)


RRHO.example <-  RRHO2(gene.list1, gene.list2, labels=c('Housekp_sig','Empirical_sig'), plots=TRUE, outputdir=plotFolder, BY=TRUE, log10.ind=TRUE)

housekp.only <- unlist(fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/RRHO2/BLBP_Lonza_d8 vs d60/housekp-only.txt', header=F), use.names = F)

housekp.only

cat(sig.in.housekp[GeneId %in% housekp.only,GeneId][1:100], sep='\n')

sig.in.housekp[GeneId %in% housekp.only,]
# set up DESeq2 for Seungmi's samples-----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/')
merged_data <- read.csv('Countfiles_gene_symbol/ISB015_merged_for_partek.csv', row.names = 1)
sampleTable <- as.data.table(readxl::read_xlsx('ISB015_seungmi_sampletable.xlsx'))
keep.samples <- names(merged_data)[names(merged_data) %in% sampleTable$replicate ]
data <- subset(merged_data, select=c(keep.samples))
data <- subset(data, row.names(data) %nin% list.new)


cts <- as.data.frame(data)
cts <- as.matrix(cts)
cts <- apply(cts, 2, as.numeric)
row.names(cts) <- row.names(data)

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/')

condition <- data.table(names = sampleTable$replicate)
condition[,condition := sampleTable$condition]
condition <- condition$condition

coldata <- data.frame(condition = condition, row.names = names(data))
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
dds

save(dds, file='ISB015.Seungmi.DDS.RData')
#
tests_table <- as.data.table(readxl::read_xlsx('ISB015_seungmi_sampletable.xlsx',
                            sheet='tests'))
tests_table

source('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')

testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/Seungmi_samples/'
volcanodir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/Seungmi_samples/Volcano_plots'

keep <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 3
dds <- dds[keep,]

cl <- makeCluster(6)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]
  
  testdir <- testdir
  volcanodir <- volcanodir
  
  DESeq2_pipeline(dds, condition1, condition2, testdir, volcanodir)
}

stopCluster(cl)

# look for ecto/endo/meso genes in ISB015-----
load('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/ISB015.Seungmi.DDS.RData')
dds

meso <- c('FGF4', 'GDF3', 'NPPB', 'NR5A2', 'PTHLH', 'T', 'ABCA4', 'ALOX15', 'BMP10', 'CDH5', 'CDX2', 'COLEC10', 'ESM1', 'FCN3', 'FOXF1', 'HAND1', 'HAND2', 'HEY1', 'HOPX', 'IL6ST', 'NKX2-5', 'ODAM', 'PDGFRA', 'PLVAP', 'RGS4', 'SNAI2', 'TBX3', 'TM4SF1', 'TBXT', 'EOMES', 'FOXC1', 'GSC', 'MIXL1', 'SNAI1', 'TBX6', 'TWIST1', 'TWIST2', 'MESD', 'NCLN', 'CFC1', 'BMP2', 'BMP4', 'BMP6', 'BMP7', 'INHBA', 'INHBB', 'INHBC', 'FABP4', 'FGF5', 'NODAL', 'TGFB1', 'TGFB2', 'TGFB3', 'WNT3A', 'WNT8A', 'CDH2', 'SST', 'KLF5', 'GATA4')

endo <- c('CABP7', 'CDH20', 'CLDN1', 'CPLX2', 'ELAVL3', 'EOMES', 'FOXA1', 'FOXA2', 'FOXP2', 'GATA4', 'GATA6', 'HHEX', 'HMP19', 'HNF1B', 'HNF4A', 'KLF5', 'LEFTY1', 'LEFTY2', 'NODAL', 'PHOX2B', 'POU3F3', 'PRDM1', 'RXRG', 'SOX17', 'SST', 'CLDN6', 'KRT19', 'FABP1', 'FABP2', 'GSC', 'SOX7', 'AFP', 'CTNNB1', 'GDF1', 'GDF3', 'MIXL1', 'SALL4')

ecto <- c('CDH9', 'COL2A1', 'DMBX1', 'DRD4', 'EN1', 'LMX1A', 'MAP2', 'MYO3B', 'NOS2', 'NR2F1', 'NR2F2', 'OLFM3', 'PAPLN', 'PAX3', 'PAX6', 'POU4F1', 'PRKCA', 'SDC2', 'SOX1', 'TRPM8', 'WNT1', 'ZBTB16', 'BMP4', 'CHRD', 'FGF8', 'FOXJ3', 'GBX1', 'NES', 'NOG', 'OTX2', 'TP63', 'PAX2', 'TUBB3', 'NCAM1', 'VIM')

library(DESeq2)
res <- lfcShrink(dds, contrast=c("condition",'iPS_Y_d60', 'iPS_CEPT_d60'))
res <- reformat.res(res, 'iPS_CEPT_d60', 'iPS_Y_d60')
res <- res[GeneId %in% c(meso, endo, ecto),]
res.meso <- res[GeneId %in% meso,]
res.meso$geneset <- 'mesoderm'
res.ecto <- res[GeneId %in% ecto,]
res.ecto$geneset <- 'ectoderm'
res.endo <- res[GeneId %in% endo,]
res.endo$geneset <- 'endoderm'

res.sig <- rbind(res.endo, res.ecto, res.meso)
res.sig <- res.sig[abs(log2FoldChange)>=1 & padj <= 0.001,]
fwrite(res.sig, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB015/DE/Seungmi_samples/ISB015_CEPT-vs-Y_trilineage_DE.csv')

# iPSC_CEPT <- dds[ , dds$condition %in% c('iPS_CEPT_d60') ]
# iPSC_CEPT_norm <- counts(iPSC_CEPT, normalized=T)
# iPSC_CEPT_norm[c('PAX6'),]
# 
# dds_norm <- counts(dds, normalized=T)
# dds_norm[c('PAX6'),]
#



#
