module load trimgalore
module load samtools
module load fastqc
module load bowtie

R1_fastq="Raw fastq with read 1"
R2_fastq="Raw fastq with read 2"
samp_name="Sample name"

# Trimming with TrimGalore
trim_galore --paired $R1_fastq $R2_fastq -o Trim_galore_output/

prefix_R1="$(echo $R1_fastq | cut -d '.' -f1)"
prefix_R2="$(echo $R2_fastq | cut -d '.' -f1)"

# Assessing Quality Check with FastQC after Trimming
fastqc -t 10 -o FastQC_after_trim Trim_galore_output/${prefix_R1}*.fq.gz
fastqc -t 10 -o FastQC_after_trim Trim_galore_output/${prefix_R2}*.fq.gz

# Alignment of reads using Bowtie2
bowtie2  --very-sensitive --end-to-end --no-mixed --no-discordant --phred33 -I 10 -X 700 -x hg38_drosophila_combined_genome -1 Trim_galore_output/${prefix_R1}*.fq.gz -2 Trim_galore_output/${prefix_R2}*.fq.gz -S ${samp_name}.sam 2> ${samp_name}_bowtie2_log.txt

# Filter reads from SAM file with mapping quality > 10
samtools view -h -F 4 -q 10 -bS ${samp_name}.sam > ${samp_name}_filtered.bam

# Sort filtered BAM file
java -jar picard.jar SortSam -I ${samp_name}_filtered.bam \
  -O ${samp_name}_filtered_sorted.bam -SORT_ORDER coordinate

# Removing reads that mapped >1 times and reads that did not map
sambamba-0.6.8 view -h -t 2 -f bam -F "[XS] == null and not unmapped and not duplicate" \
  ${samp_name}_filtered_sorted.bam > \
  ${samp_name}_filtered_sorted_unique.bam

# Removing duplicate reads
java -jar picard.jar MarkDuplicates I=${samp_name}_filtered_sorted_unique.bam \
  O=${samp_name}_filtered_noDUP_unique.bam REMOVE_DUPLICATES=true \
  METRICS_FILE=${samp_name}_filtered_noDUP.txt

# Filtering regions from the ENCODE black list using bedtools
bedtools intersect -v -a ${samp_name}_filtered_noDUP_unique.bam -b hg38_encode_blacklist_chip-seq.bed > ${samp_name}_filtered_noDUP_unique.noB_unique.bam

# Index of the filtered BAM file with all the reads (Mouse and Drosophila)
samtools index ${samp_name}_filtered_noDUP_unique.noB_unique.bam

# Count total reads in BAM file
samtools view -c ${samp_name}_filtered_noDUP_unique.noB_unique.bam > ${samp_name}_count_reads_filtered_TOTAL.txt

# Count mouse reads in BAM file
samtools view ${samp_name}_filtered_noDUP_unique.noB_unique.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY | wc -l > ${samp_name}_count_reads_filtered_mouse.txt

# Count reads from Drosophila chromosomes (spikein reads) in the filtered BAM. Note that chromosomes 4, X, Y 
# in Drosophila differ from chr4, chrX and chrY from mouse (they lack the 'chr' prefix)
samtools view ${samp_name}_filtered_noDUP_unique.noB_unique.bam 2L 2R 3L 3R 4 X Y | wc -l > ${samp_name}_count_reads_filtered_spikein.txt

# Subset the filtered BAM file for mouse reads and Drosophila reads
samtools view -b ${samp_name}_filtered_noDUP_unique.noB_unique.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ${samp_name}_filtered_mouse_reads.bam
samtools view -b ${samp_name}_filtered_noDUP_unique.noB_unique.bam 2L 2R 3L 3R 4 X Y > ${samp_name}_filtered_Drosophila_reads.bam

# Make index for filtered BAM file with mouse reads
samtools index ${samp_name}_filtered_human_reads_final.bam

# Creating a bigwig file from the BAM file with mouse reads using Deeptools
bamCoverage --bam ${samp_name}_filtered_human_reads_final.bam -o ${samp_name}_filtered_human_reads.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --ignoreDuplicates \
    --ignoreForNormalization chrX chrY \
    --numberOfProcessors 10 \
    --extendReads \
    --centerReads
