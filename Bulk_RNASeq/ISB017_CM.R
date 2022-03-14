library(data.table)
library(readxl)
library(stringr)
library(DESeq2)
library(Hmisc)
library(ggplot2)
library(foreach)
library(doParallel)
library(ggrepel)
library(ggthemes)

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/Countfiles_gene_symbol/')

files <- Sys.glob('*.txt')
files <- files[-grep('^CHU', files)]
files <- files[-45]
files

genes_filter <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/35473_Gene_To_Filter_PF.csv')

for (file in files[1]){
  data <- fread(file, header=F)
  samplename <- str_split_fixed(file, '_gene_symbol_counts.txt', n=Inf)[1]
  samplename <- str_split_fixed(samplename, '_htseq', n=Inf)
  samplename <- paste(samplename[2:length(samplename)],collapse="_")
  names(data) <- c('GeneId', paste0(samplename))
  
  genes <- data$GeneId[data$GeneId %nin% genes_filter$GeneSymbol]
  data <- data[GeneId %in% genes,]
}


merged_data <- data.table('GeneId'=data$GeneId)

for (file in files){
  data <- fread(file, header=F)
  samplename <- str_split_fixed(file, '_htseq_counts.txt', n=Inf)[1]
  #samplename <- str_split_fixed(samplename, '_htseq', n=Inf)[1]
  #samplename <- paste(samplename[2:length(samplename)],collapse="_")
  names(data) <- c('GeneId', paste0(samplename))
  
  data <- data[GeneId %in% genes,]
  
  merged_data <- merge(merged_data, data, by='GeneId', all=T)
  
}

fwrite(merged_data, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/ISB017_merged_for_partek.csv', sep=',', col.names = T, row.names = F, quote=T)

# set up DESeq-----

#data <- read.csv('ISB017_merged_for_partek.csv', row.names = 1)
data <- as.data.frame(merged_data)
data[1:10,1:10]
row.names(data) <- data$GeneId
data <- data[-1]
data[1:10,1:10]

cts <- as.data.frame(data)
cts <- as.matrix(cts)
cts <- apply(cts, 2, as.numeric)
row.names(cts) <- row.names(data)

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/')

sampletable <- as.data.table(read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/ISB017_samplesheet_CM.xlsx', sheet=1))
sampletable[,replicate:= paste(condition, replicate, sep="_")]
sampletable
cat(sampletable$replicate, sep='\n')

condition <- data.table(names = sampletable$replicate)
condition[,condition := sampletable$condition]
condition <- condition$condition

coldata <- data.frame(condition = condition, row.names = names(data))
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
dds

save(dds, file='/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/ISB017.DDS.CM.RData')

# DE tests-----
source('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')

tests_table <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/ISB017_samplesheet_CM.xlsx', sheet=2))

testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/CM/'
volcanodir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/CM/volcano'


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

# combine all test results-----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/CM/')
files <- Sys.glob('DE*.csv')
files
ex <- fread(files[1], header=T)
genelist <- ex$GeneId

data.merged <- fread(files[1], header=T)
condition1 <- str_split_fixed(file, '\\.', Inf)[2]
condition2 <- str_split_fixed(file, '\\.', Inf)[4]
data.merged[,test:= paste0(condition1, '.vs.', condition2)]
names(data.merged)[2] <- 'Condition1_mean'
names(data.merged)[3] <- 'Condition2_mean'
data.merged

for (file in files[2:15]){
  data <- fread(file, header=T)
  condition1 <- str_split_fixed(file, '\\.', Inf)[2]
  condition2 <- str_split_fixed(file, '\\.', Inf)[4]
  data[,test:= paste0(condition1, '.vs.', condition2)]
  names(data)[2] <- 'Condition1_mean'
  names(data)[3] <- 'Condition2_mean'
  data.merged <- rbind(data.merged, data)
}

data.merged
data.merged <- data.merged[abs(log2FoldChange) >1 & padj <=0.01, ]
data.merged[,c('Condition1','Condition2'):= tstrsplit(test,'.vs.')]
data.merged


# merged volcano plot----
ggplot(data=data.merged, aes(x=log2FoldChange, y=-log10(padj), color=test)) + geom_point() + theme_bw() +
  labs(title='Combined volcano plot of ISB017 tests', color='DE test')+
  xlim(-12,12)+
  geom_text_repel(data=data.merged.label.7, aes(x=log2FoldChange, y=-log10(padj), label=Gene_short))+
  guides(label=FALSE) + coord_cartesian(ylim=c(0,40))

# try DESeq2 with all 3D cept vs. y-----

data <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/ISB017_merged_for_partek.csv')
data[1:10,1:10]

sampletable <- as.data.table(read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/ISB017_samplesheet_CM.xlsx', sheet=3))
sampletable

row.names(data) <- data$GeneId
data <- data[-1]
data[1:10,1:10]
all(names(data)[2:length(data)] == sampletable$replicate_old)
names(data)[2:length(data)] <- sampletable$replicate_new
data[1:10,1:10]
data <- as.data.frame(data)

cts <- as.data.frame(data[2:45])
cts <- as.matrix(cts)
cts <- apply(cts, 2, as.numeric)
row.names(cts) <- data$GeneId

condition <- data.table(names = sampletable$replicate_new)
condition[,condition := sampletable$condition]
condition <- condition$condition

coldata <- data.frame(condition = condition, row.names = names(data)[2:length(data)])
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
dds

save(dds, file='/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/ISB017.DDS.2Dvs3D.CM.RData')
# DE for 3D Y vs 3D CEPT-----
tests_table <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/ISB017_samplesheet_CM.xlsx', sheet=4))
testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/CM/3DCEPT_vs_Y/'
volcanodir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/DE/CM/3DCEPT_vs_Y/volcano'
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


# output for partek-----

raw.counts <- as.data.frame(counts(dds, normalized=F))
raw.counts
write.csv(raw.counts, file='/Volumes/ncatssctl/NGS_related/BulkRNA/ISB017/Analysis/ISB017_raw_counts_Partek.csv')
