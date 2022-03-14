library(data.table)
setwd('//128.231.11.251/ncatssctl/NGS_related/RASLseq/ISRL010/')
rasl <- fread('RASL010_merged_raw_counts.csv')
rasl

setwd('Separate_count_files_ISRL010/')
for (i in names(rasl)[2:49]){
  #data <- subset(rasl, select=c('V1', i))
  #names(data)[1] <- 'oligoID'
  #fwrite(data, paste0('RASL010_',i,'_raw_counts.csv'), sep = ',', quote=F, col.names = T)
  cat(paste0('RASL010_',i,'_raw_counts.csv.md5sum'))
  cat('\n')
}
