setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/Comparison_with_LancasterD60/')
library(data.table)
library(stringr)


# samples-----
samples.paired <- c('7_H9ES_P31_CEPT_d60_S61','8_H9ES_P31_CEPT_d60_S62','9_H9ES_P31_CEPT_d60_S63')
samples.single <- c('SRR3416403','SRR3416404','SRR3416405','SRR3416415','SRR3416416','SRR3416417')

samples.paired.control <- c('ISB_003_25_S25','ISB_003_26_S26','ISB_003_27_S27', 'ISB_003_28_S28')

# samples <- samples.paired
# samples <- samples.single

# chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
# 
# chunks <- chunk2(samples, 1)

#trimmomatic-----

chunk <- data.table('chunk'=unlist(samples.paired))

chunk[,command:= paste0('java -jar $TRIMMOJAR PE -phred33 /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/',chunk,'_R1_001.fastq.gz /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/',chunk,'_R2_001.fastq.gz -baseout /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/',chunk,'.trim.fastq.gz ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36')]

write(as.character(chunk$command), file=paste0('Swarm_paired_trim.sh'), append=TRUE)

#trimmomatic, paired, for controls----
chunk <- data.table('chunk'=unlist(samples.paired.control))

chunk[,command:= paste0('java -jar $TRIMMOJAR PE -phred33 /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/',chunk,'_R1_001.fastq.gz /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/',chunk,'_R2_001.fastq.gz -baseout /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/',chunk,'.trim.fastq.gz ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36')]

write(as.character(chunk$command), file=paste0('Swarm_paired_control_trim.sh'), append=TRUE)

# swarm -f Swarm_paired_control_trim.sh -g 5 --time=10:00:00 --module trimmomatic

#single
chunk <- data.table('chunk'=unlist(samples.single))

chunk[,command:= paste0('java -jar $TRIMMOJAR SE -phred33 /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/',chunk,'_1.fastq /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/',chunk,'.trim.fastq.gz ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:10 TRAILING:5 MAXINFO:50:0.97 MINLEN:36')]

write(as.character(chunk$command), file=paste0('Swarm_single_trim.sh'), append=TRUE)

# star and htseq, paired-----

chunk <- data.table('chunk'=unlist(samples.paired))

chunk[,command:= paste0('cd ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/ && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/', chunk , '.trim_1P.fastq.gz ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/', chunk , '.trim_2P.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/', chunk, '_hg38 && htseq-count -f bam -r pos -s no -t exon -m union -i gene_name /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/', chunk  ,'_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/htseq/', chunk , '_htseq_counts.txt')]

write(as.character(chunk$command), file=paste0('Swarm_paired_align_htseq.sh'), append=TRUE)

cat(paste0('swarm -f Swarm_paired_align_htseq.sh -g 35 --time=48:00:00 --gres=lscratch:200 --module STAR,htseq'))

# star and htseq, paired, controls-----

chunk <- data.table('chunk'=unlist(samples.paired.control))

chunk[,command:= paste0('cd ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/ && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/', chunk , '.trim_1P.fastq.gz ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/', chunk , '.trim_2P.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/', chunk, '_hg38 && htseq-count -f bam -r pos -s no -t exon -m union -i gene_name /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/', chunk  ,'_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/htseq/', chunk , '_htseq_counts.txt')]

write(as.character(chunk$command), file=paste0('Swarm_paired_control_align_htseq.sh'), append=TRUE)

cat(paste0('swarm -f Swarm_paired_control_align_htseq.sh -g 35 --time=48:00:00 --gres=lscratch:200 --module STAR,htseq'))



# star and htseq, single-----
chunk <- data.table('chunk'=unlist(samples.single))

chunk[,command:= paste0('cd ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/ && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/', chunk , '.trim.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/', chunk, '_hg38 && htseq-count -f bam -r pos -s no -t exon -m union -i gene_name /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/', chunk  ,'_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > ', '/data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/htseq/', chunk , '_htseq_counts.txt')]

write(as.character(chunk$command), file=paste0('Swarm_single_align_htseq.sh'), append=TRUE)

cat(paste0('swarm -f Swarm_single_align_htseq.sh -g 35 --time=48:00:00 --gres=lscratch:200 --module STAR,htseq'))

# start DESeq----
library(Hmisc)
library(data.table)
library(DESeq2)
library(stringr)
library(readxl)
library(RUVSeq)
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/Comparison_with_LancasterD60/htseq/')

files <- Sys.glob('*.txt')
files

genes.exclude <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/28309_gene_IDs_to_exclude_full_list.txt', header=F)
genes.exclude2 <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/Comparison_with_LancasterD60/genes_also_remove.txt', header=F)

for (file in files[1]){
  data <- fread(file, header=F)
  samplename <- str_split_fixed(file, '_gene_symbol_counts.txt', n=Inf)[1]
  samplename <- str_split_fixed(samplename, '_', n=Inf)
  samplename <- paste(samplename[2:length(samplename)],collapse="_")
  names(data) <- c('GeneId', paste0(samplename))

  genes <- data$GeneId
  genes <- genes[-grep('^LINC', genes)]
  genes <- genes[genes %nin% genes.exclude$V1]
  genes <- genes[genes %nin% genes.exclude2$V1]
  
  data <- data[GeneId %in% genes,]
}


merged_data <- data.table('GeneId'=data$GeneId)

for (file in files){
  data <- fread(file, header=F)
  samplename <- str_split_fixed(file, '_htseq_counts.txt', n=Inf)[1]
  samplename <- str_split_fixed(samplename, '_', n=Inf)
  samplename <- paste(samplename[2:length(samplename)],collapse="_")
  names(data) <- c('GeneId', paste0(samplename))
  
  data <- data[GeneId %in% genes,]
  
  merged_data <- merge(merged_data, data, by='GeneId', all=T)
  
}
names(merged_data)

names(merged_data)[5:8] <- c('H9E8_1','H9E8_2','H9E8_3','H9E8_4')
names(merged_data)[9:14] <- gsub('NA_','CerebralOrg_D60_', names(merged_data)[9:14])

fwrite(merged_data, 'ISB013_Lancaster_merged_for_partek.csv', sep=',', col.names = T, row.names = F, quote=T)

data <- as.data.frame(merged_data)
row.names(data) <- data$GeneId
data <- data[-1]

cts <- as.data.frame(data)
cts <- as.matrix(cts)
cts <- apply(cts, 2, as.numeric)
row.names(cts) <- row.names(data)

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/Comparison_with_LancasterD60/')

sampletable <- as.data.table(read_xlsx('sampletable.xlsx', sheet=1))

condition <- data.table(names = sampletable$replicate)
condition[,condition := sampletable$condition]
condition <- condition$condition

coldata <- data.frame(condition = condition, row.names = names(data))
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
dds

save(dds, file='ISB013_Lancaster.DDS.RData')

# batch correction----
library(RColorBrewer)
sampletable

counts <- data
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



# set up new dds with batch correction factor included in the design----
condition <- x

names(pData(set2)) <- c('condition', 'W_1')

dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)

save(dds, file='ISB013_Lancaster.DDS.after_housekp_correction.RData')

# DE test----
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
dds <- DESeq(dds)

teststable <- as.data.table(readxl::read_xlsx('sampletable.xlsx', sheet=2))
condition1 <- teststable$condition1[3]
condition2 <- teststable$condition2[3]

res <- lfcShrink(dds, contrast=c("condition", condition2, condition1)) 

reformat.res <- function(x, condition1, condition2){
  x <- as.data.frame(x)
  x$GeneId <- row.names(x)
  x <- as.data.table(x)
  x[,baseMean:=NULL]
  x[,baseMean1 := rowMeans(as.data.frame(counts(dds,normalized=TRUE)[,dds$condition == condition1]))]
  x[,baseMean2 := rowMeans(as.data.frame(counts(dds,normalized=TRUE)[,dds$condition == condition2]))]
  x <- na.omit(x)
  x <- x[order(rank(padj))]
  x <- x[,c("GeneId", "baseMean1", "baseMean2", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  names(x)[2:3] <- c(condition1, condition2)
  x <- x[order(-abs(log2FoldChange), padj)]
  return(x)
}

res <- reformat.res(res, condition1, condition2)
res <- res[padj < 0.001,]
fwrite(res, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/Comparison_with_LancasterD60/H9E8_vs_CerebralOrg_D60_DE.csv')


# 
fwrite(merged_data, 'ISB013_Lancaster_merged_for_partek.csv', sep=',', col.names = T, row.names = F, quote=T)

dds.corrected <- dds

load('ISB013_Lancaster.DDS.RData')

dds.before.correct <- dds

head(counts(dds.before.correct, normalized=T))

counts.norm.before <- counts(dds.before.correct, normalized=T)
counts.norm.after <- counts(dds.corrected, normalized=T)

all.equal(counts.norm.before, counts.norm.after)

intersect(row.names(counts.norm.after), row.names(counts.norm.before))

counts.norm.after <- subset(counts.norm.after, row.names(counts.norm.after) %in% intersect(row.names(counts.norm.after), row.names(counts.norm.before)))
counts.norm.before <- subset(counts.norm.before, row.names(counts.norm.before) %in% intersect(row.names(counts.norm.after), row.names(counts.norm.before)))

all.equal(counts.norm.before, counts.norm.after)

head(counts.norm.before)
head(counts.norm.after)

#