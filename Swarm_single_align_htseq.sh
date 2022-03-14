cd /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/ && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/SRR3416403.trim.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416403_hg38 && htseq-count -f bam -r pos -s no -t exon -m union -i gene_name /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416403_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/htseq/SRR3416403_htseq_counts.txt
cd /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/ && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/SRR3416404.trim.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416404_hg38 && htseq-count -f bam -r pos -s no -t exon -m union -i gene_name /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416404_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/htseq/SRR3416404_htseq_counts.txt
cd /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/ && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/SRR3416405.trim.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416405_hg38 && htseq-count -f bam -r pos -s no -t exon -m union -i gene_name /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416405_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/htseq/SRR3416405_htseq_counts.txt
cd /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/ && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/SRR3416415.trim.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416415_hg38 && htseq-count -f bam -r pos -s no -t exon -m union -i gene_name /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416415_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/htseq/SRR3416415_htseq_counts.txt
cd /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/ && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/SRR3416416.trim.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416416_hg38 && htseq-count -f bam -r pos -s no -t exon -m union -i gene_name /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416416_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/htseq/SRR3416416_htseq_counts.txt
cd /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/ && STAR --runThreadN $SLURM_CPUS_PER_TASK --genomeDir /fdb/STAR_current/GENCODE/Gencode_human/release_27/genes-75 --sjdbOverhang 75 --outSAMunmapped Within --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesIn /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Trimmed_FASTQs/SRR3416417.trim.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416417_hg38 && htseq-count -f bam -r pos -s no -t exon -m union -i gene_name /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/Aligned_BAMs/SRR3416417_hg38Aligned.sortedByCoord.out.bam /fdb/GENCODE/Gencode_human/release_27/gencode.v27.annotation.gtf > /data/NCATS_ifx/data/mRNASeq/ISB013/sra_merge/htseq/SRR3416417_htseq_counts.txt
