library(data.table)
library(readxl)
library(stringr)
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/')
samples <- fread('samples.txt', header=F)

indir <- '/data/NCATS_ifx/data/mRNASeq/ISB013/'
outdir <- '/data/NCATS_ifx/data/mRNASeq/ISB013/htseq'
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/')


for (sample in samples){
  sample <- gsub('_star.genome.sorted.bam', '', sample)
  
  cat(paste0('htseq-count -f bam -r pos -s no -t exon -m union ', indir , 'bam/' , sample  , '_star.genome.sorted.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > ', outdir, '/', sample , '_htseq_counts.txt' , "\n\n"))
}

# swarm -f htseq_swarm.sh -t 4 -g 4 --module htseq --time=10:00:00

# convert to gene symbols----
for (file in files[24:27]){
  
  dt <- fread(file)
  dt <- dt[1:(nrow(dt)-5),]
  names(dt) <- c("ENSG.full", "Counts")
  
  ensg.genes <- data.table("ENSG.full" = dt$ENSG.full)
  
  ensg.genes[,ENSG.short := tstrsplit(ENSG.full, "\\.")[1]]
  
  genes <- ensg.genes$ENSG.short
  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest")
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
  G_list <- as.data.table(G_list)
  ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG.short", by.y="ensembl_gene_id")
  ensg.genes <- na.omit(ensg.genes)
  dt <- subset(dt, dt$ENSG.full %in% ensg.genes$ENSG.full)
  dt <- merge(dt, ensg.genes, by="ENSG.full")
  dt <- dt[,c("external_gene_name", "Counts")]
  dt <- dt[order(external_gene_name, -Counts)]
  dt <- dt[!duplicated(external_gene_name),]
  
  sample_id <- str_split_fixed(file, "_htseq_counts.txt",2)[1]
  
  fwrite(dt, paste0("/data/NCATS_ifx/data/mRNASeq/ISB013/htseq/Countfiles_gene_symbol/", sample_id, "_gene_symbol_counts.txt"), col.names = F, row.names = F, quote=F, sep="\t")
  
}

# 013120 make DDS for all samples----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/Countfiles_gene_symbols/')

sampletable <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/ISB013_sampletable.xlsx'))
sampletable
directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/Countfiles_gene_symbols/"
dds <- DESeqDataSetFromHTSeqCount(sampletable,
                                  directory = directory,
                                  design = ~ condition)

dds <- DESeq(dds)
save(dds, file='/Volumes/ncatssctl/NGS_related/BulkRNA/ISB013/ISB013.DDS.RData')

