library(data.table)
library(readxl)
library(DESeq2)
library(stringr)
library(doParallel)
library(foreach)
library(Hmisc)
setwd('/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/')

files <- fread('ISRL010_split_files.txt', header=F)

samplesheet <- as.data.table(readxl::read_xlsx('ISRL010_rerun_samplesheet.xlsx'))
samplesheet
length(unique(samplesheet$wellID))

wells.keep <- samplesheet[condition !='NTC',]

files[,well_S:= tstrsplit(V1, '_L')[1]]
files[,well:= tstrsplit(V1, '_')[1]]

files.keep <- files[well %in% wells.keep$wellID,]

for (well in unique(files.keep$well_S)){
  cat('cat ', well,'_L001_R1_001.fastq.gz ',well,'_L002_R1_001.fastq.gz ',well,'_L003_R1_001.fastq.gz ',well,'_L004_R1_001.fastq.gz > /data/NCATS_ifx/data/rasl/ISRL010/data/',well,'_R1_001.fastq.gz' ,'\n\n', sep='')
}

# subset data by well ID-----
data <- fread('/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/ISRL010_count_tabs.csv', header=T, skip=14L)
data

samplesheet

names(data)

colinfo <-data.table(cols=names(data))
colinfo[,wellID:=tstrsplit(cols, '_S')[1]]
colinfo <- merge(colinfo, wells.keep, by='wellID', all=T)
colinfo

# 2 wells have pct_match < 30%.
#G18_S270_R1_001.fastq.gz
#A18_S198_R1_001.fastq.gz

colinfo <- colinfo[wellID %nin% c('G18', 'A18'),]

data <- subset(data, select=c(colinfo$cols))

data <- as.data.frame(data)
row.names(data) <- data$oligo_id
data <- data[,-127]
data

data.backup <- data
data.backup$RowSums <- raster::rowSums(data.backup)
data.backup <- subset(data.backup, data.backup$RowSums >= 30)

data <- data.backup

exp1 <- as.data.table(readxl::read_xlsx('ISRL010_rerun_samplesheet.xlsx', sheet=2))
exp1

exp2.1 <- as.data.table(readxl::read_xlsx('ISRL010_rerun_samplesheet.xlsx', sheet=3))
exp2.1

exp2.2 <- as.data.table(readxl::read_xlsx('ISRL010_rerun_samplesheet.xlsx', sheet=4))
exp2.2

data <- data[-127]

fwrite(data, 'ISRL010_all_raw_counts.csv', sep=',', col.names = T, row.names = T, quote=T)
fwrite(colinfo, 'colinfo.csv', sep=',', col.names = T, row.names = F, quote=T)

# start DESeq2-----

exp1.wells <- colinfo
exp1.wells <- exp1.wells[wellID %in% exp1$wellID,]

exp2.1.wells <- colinfo
exp2.1.wells <- exp2.1.wells[wellID %in% exp2.1$wellID,]

# exp2.2 will have ISRL009 combined with ISRL010.
# exp2.2.wells <- colinfo
# exp2.2.wells <- exp2.2.wells[wellID %in% exp2.2$wellID,]

exp1.wells

data.exp1 <- subset(data, select=c(names(data) %in% exp1.wells$cols))
all(names(data.exp1) %in% exp1.wells$cols)
names(data.exp1) <- exp1.wells$replicate
cts <- as.data.frame(data.exp1)
cts <- as.matrix(cts)
cts <- apply(cts, 2, as.numeric)
row.names(cts) <- row.names(data.exp1)

condition <- data.table(names = names(data.exp1))
condition[,condition := exp1.wells$condition]
condition <- condition$condition

coldata <- data.frame(condition = condition, row.names = names(data.exp1))
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
dds

save(dds, file='/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/RASL010_exp1.RData') # no threshold for minimum counts here.

dds.data <- as.data.frame(counts(dds, normalized=T))
dds.data.raw <- as.data.frame(counts(dds, normalized=F))

fwrite(dds.data, '/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/RASL010_exp1_normalized_counts.csv', sep=',', quote=F, row.names = T, col.names = T)
fwrite(dds.data.raw, '/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/RASL010_exp1_raw_counts.csv', sep=',', quote=F, row.names = T, col.names = T)
# exp2-----
data.exp2.1 <- subset(data, select=c(names(data) %in% exp2.1.wells$cols))
all(names(data.exp2.1) %in% exp2.1.wells$cols)
names(data.exp2.1) <- exp2.1.wells$replicate
cts <- as.data.frame(data.exp2.1)
cts <- as.matrix(cts)
cts <- apply(cts, 2, as.numeric)
row.names(cts) <- row.names(data.exp2.1)

condition <- data.table(names = names(data.exp2.1))
condition[,condition := exp2.1.wells$condition]
condition <- condition$condition

coldata <- data.frame(condition = condition, row.names = names(data.exp2.1))
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
dds

save(dds, file='/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/RASL010_exp2.1.RData') # no threshold for minimum counts here.

dds.data <- as.data.frame(counts(dds, normalized=T))
dds.data.raw <- as.data.frame(counts(dds, normalized=F))

fwrite(dds.data, '/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/RASL010_exp2.1_normalized_counts.csv', sep=',', quote=F, row.names = T, col.names = T)
fwrite(dds.data.raw, '/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/RASL010_exp2.1_raw_counts.csv', sep=',', quote=F, row.names = T, col.names = T)


# exp2.2-----
# columns should be named as their replicate name not fastq file name.
load('/Volumes/ncatssctl/NGS_related/RASLseq/RASLResults/R_analysis/RASL009_DDS.RData')
rasl9 <- dds

rasl9.raw <- as.data.frame(counts(dds, normalized=F))
#rasl9.raw$GeneId <- row.names(rasl9.raw)

rasl9.wells <- exp2.2[ISRL =='ISRL009',]

rasl9.raw <- subset(rasl9.raw, select=c(names(rasl9.raw) %in% rasl9.wells$replicate))
rasl9.raw$GeneId <- row.names(rasl9.raw)


data.exp2.2 <- data.exp2.1
data.exp2.2$GeneId <- row.names(data.exp2.2)

#data.exp2.2 <- merge(data.exp2.2, rasl9.raw, by='GeneId', all=F)

colinfo.exp2.2 <- as.data.table(readxl::read_xlsx('exp2.2.xlsx', sheet=1))
colinfo.exp2.2

names(rasl9.raw) <- c('H9_ISRL009_1', 'H9_ISRL009_2', 'H9_ISRL009_3', 'H9_ISRL009_4', 'H9_ISRL009_Ecto_1', 'H9_ISRL009_Ecto_2', 'H9_ISRL009_Ecto_3', 'H9_ISRL009_Ecto_4', 'H9_ISRL009_Endo_1', 'H9_ISRL009_Endo_2', 'H9_ISRL009_Endo_3', 'H9_ISRL009_Endo_4', 'H9_ISRL009_Meso_1', 'H9_ISRL009_Meso_2', 'H9_ISRL009_Meso_3', 'H9_ISRL009_Meso_4','GeneId')

names(data.exp2.2) <- c('CEPT-E6D7_1', 'Y27-E6D7_1', 'CEPT-E6D7_2', 'Y27-E6D7_2', 'CEPT-E6D7_3', 'Y27-E6D7_3', 'CEPT-E6D7_4', 'Y27-E6D7_4', 'CEPT-E6D7_5', 'Y27-E6D7_5', 'H9_ISRL010_1', 'H9_ISRL010_2', 'H9_ISRL010_3', 'H9_ISRL010_4', 'CEPT-E6D7_6', 'Y27-E6D7_6', 'CEPT-E6D7_7', 'Y27-E6D7_7', 'CEPT-E6D7_8', 'Y27-E6D7_8', 'CEPT-E6D7_9', 'Y27-E6D7_9', 'CEPT-E6D7_10', 'Y27-E6D7_10', 'CEPT-E6D7_11', 'Y27-E6D7_11', 'CEPT-E6D7_12', 'Y27-E6D7_12', 'CEPT-E6D7_13', 'Y27-E6D7_13', 'H9_ISRL010_Ecto_1', 'H9_ISRL010_Ecto_2', 'H9_ISRL010_Ecto_3', 'H9_ISRL010_Ecto_4', 'CEPT-E6D7_14', 'Y27-E6D7_14', 'H9_ISRL010_Endo_1', 'H9_ISRL010_Endo_2', 'H9_ISRL010_Endo_3', 'H9_ISRL010_Endo_4', 'CEPT-E6D7_15', 'Y27-E6D7_15', 'H9_ISRL010_Meso_1', 'H9_ISRL010_Meso_2', 'H9_ISRL010_Meso_3', 'H9_ISRL010_Meso_4', 'CEPT-E6D7_16', 'Y27-E6D7_16','GeneId')

data.exp2.2 <- merge(data.exp2.2, rasl9.raw, by='GeneId', all=F)

all(names(data.exp2.2) %in% colinfo.exp2.2$replicate)
data.exp2.2<- as.data.frame(data.exp2.2)
row.names(data.exp2.2) <- data.exp2.2$GeneId
data.exp2.2 <- data.exp2.2[-1]

cts <- as.data.frame(data.exp2.2)
cts <- as.matrix(cts)
cts <- apply(cts, 2, as.numeric)
row.names(cts) <- row.names(data.exp2.2)

condition <- data.table(names = names(data.exp2.2))
condition[,condition := colinfo.exp2.2$condition[2:65]]
condition <- condition$condition

coldata <- data.frame(condition = condition, row.names = names(data.exp2.2))
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
dds

save(dds, file='/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/RASL010_exp2.2_before_correction.RData') # no threshold for minimum counts here.

dds.data <- as.data.frame(counts(dds, normalized=T))
dds.data.raw <- as.data.frame(counts(dds, normalized=F))

fwrite(dds.data, '/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/RASL010_exp2.2_normalized_counts.csv', sep=',', quote=F, row.names = T, col.names = T)
fwrite(dds.data.raw, '/Volumes/ncatssctl/NGS_related/RASLseq/ISRL010/RASL010_exp2.2_raw_counts.csv', sep=',', quote=F, row.names = T, col.names = T)

# triangle plot, 05/24/21----
setwd('//128.231.11.251/ncatssctl/NGS_related/RASLseq/ISRL010')
#load('RASL010_exp2.1.RData')
data <- fread("//128.231.11.251/ncatssctl/NGS_related/RASLseq/RASLResults/RASLSeq010_triangleplot/20210518_ISRL010_EB_analysis/Filtered_L2FC_EB_genelevel_433gene_ISRL010.txt")
names(data)[1] <- 'gene'


#View(as.data.frame(dds@colData))
#data <- as.data.frame(counts(dds, normalized=T))

# genes.dt <- data.table(oligo=row.names(data))
# genes.dt[,GeneId:=tstrsplit(oligo, '_')[2]]
# genes.dt[,sum:=rowSums(data, na.rm=T)]
# genes.dt <- genes.dt[order(GeneId, -sum)]
# 
# genes.keep <- genes.dt[, head(.SD, 1), keyby = GeneId]
# oligos.keep <- genes.keep$oligo
# 
# data.keep <- subset(data, row.names(data) %in% oligos.keep)
# length(row.names(data.keep))
# length(genes.keep$GeneId)
# 
# row.names(data.keep) <- genes.keep$GeneId
# 
# data.keep$gene <- row.names(data.keep)
# data.keep <- as.data.table(data.keep)

marker.lists <- as.data.table(readxl::read_xlsx('Trilineage_markers_validated_probe_RASLSeq.xlsx'))
meso.genes <- unique(na.omit(marker.lists$meso))
endo.genes <- unique(na.omit(marker.lists$endo))
ecto.genes <- unique(na.omit(marker.lists$ecto))

meso.dt <- data[gene %in% meso.genes,]
endo.dt <- data[gene %in% endo.genes,]
ecto.dt <- data[gene %in% ecto.genes,]
meso.means <- colMeans(meso.dt[,-1])
endo.means <- colMeans(endo.dt[,-1])
ecto.means <- colMeans(ecto.dt[,-1])

data.tern <- rbind(meso.means, endo.means, ecto.means)
data.tern <- as.data.table(t(data.tern))
data.tern

names(data.tern) <- c('Mesoderm','Endoderm','Ectoderm')

# metadata <- as.data.frame(dds@colData)
# metadata$replicate <- row.names(metadata)
# metadata <- as.data.table(metadata)
# metadata

metadata <- data.table(replicate=names(data)[-1])
metadata$replicate <- gsub('H9_', 'H9-', metadata$replicate)
metadata[,condition:= tstrsplit(replicate, '_')[1]]

data.tern[,Sample:= metadata$condition]

data.tern$Mesoderm <- data.tern$Mesoderm+2
data.tern$Endoderm <- data.tern$Endoderm+2
data.tern$Ectoderm <- data.tern$Ectoderm+2

data.tern$Sample <- factor(data.tern$Sample, levels=c('H9-Ecto', 'H9-Endo','H9-Meso',
                                               'CEPT-E6D7','Y27-E6D7'))

ggtern(data.tern, aes(Mesoderm, Endoderm, Ectoderm, color=Sample)) + 
  theme_bw() +
  geom_point(size=1)+
  limit_tern(0.6,0.7,0.78)

#save(data.tern, file='Data.ternaryplot.RData')

# get means for all Y27 or CEPT

data.tern.means <- data.table(meso.means=c(mean(data.tern[Sample=='CEPT-E6D7',Mesoderm]),
                                           mean(data.tern[Sample=='Y27-E6D7',Mesoderm])),
                              endo.means=c(mean(data.tern[Sample=='CEPT-E6D7',Endoderm]),
                                           mean(data.tern[Sample=='Y27-E6D7',Endoderm])),
                              ecto.means=c(mean(data.tern[Sample=='CEPT-E6D7',Ectoderm]),
                                           mean(data.tern[Sample=='Y27-E6D7',Ectoderm])),
                              Sample=c('Average CEPT',
                                       'Average Y27'))
data.tern.means
names(data.tern.means)[1:3] <-c('Mesoderm','Endoderm','Ectoderm')

data.tern.both <- rbind(data.tern, data.tern.means)

ggtern(data.tern.both, aes(Mesoderm, Endoderm, Ectoderm, color=Sample)) + 
  theme_bw() +
  geom_point(size=c(rep(1, 44), rep(4, 2)))+
  limit_tern(0.6,0.7,0.78)

save(data.tern.both, file='Data.ternaryplot.RData')
