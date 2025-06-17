# Required packages:
# - Trim Galore
# - Fastqc
# - Bowtie2
# - Samtools
# - BEDtools
# - Sambamba
# - deepTools


R1_fastq="Raw fastq file with read 1"
R2_fastq="Raw fastq file with read 2"
samp_name="Sample name"

# Trimming adapter sequences with Trim Galore from the raw fastq files
trim_galore --paired $R1_fastq $R2_fastq

prefix_R1="$(echo $R1_fastq | cut -d '.' -f1)"
prefix_R2="$(echo $R2_fastq | cut -d '.' -f1)"

# Performing quality control check of the sequencing data
fastqc -t 10 ${prefix_R1}*.fq.gz
fastqc -t 10 ${prefix_R2}*.fq.gz

# Alignment of reads against mouse (mm10) and Drosophila/spikein (dm6) reference genomes using Bowtie2. The genomes of both species are concatenated in a single fasta file.
bowtie2  --very-sensitive --end-to-end --no-mixed --no-discordant --phred33 -I 10 -X 700 -x mm10_and_dm6_genomes.fasta -1 ${prefix_R1}*.fq.gz -2 ${prefix_R2}*.fq.gz -S ${samp_name}.sam 2> ${samp_name}_bowtie2_log.txt

# Filtering reads from SAM file with mapping quality > 10
samtools view -h -F 4 -q 10 -bS ${samp_name}.sam > ${samp_name}_filtered.bam

# Sorting filtered BAM file by coordinates
java -jar picard.jar SortSam -I ${samp_name}_filtered.bam \
  -O ${samp_name}_filtered_sorted.bam -SORT_ORDER coordinate

# Removing reads that did not map and reads that mapped >1 times
sambamba-0.6.8 view -h -t 2 -f bam -F "[XS] == null and not unmapped and not duplicate" \
  ${samp_name}_filtered_sorted.bam > \
  ${samp_name}_filtered_sorted_unique.bam

# Removing duplicate reads
java -jar picard.jar MarkDuplicates I=${samp_name}_filtered_sorted_unique.bam \
  O=${samp_name}_filtered_noDUP_unique.bam REMOVE_DUPLICATES=true \
  METRICS_FILE=${samp_name}_filtered_noDUP.txt

# Filtering regions from the ENCODE black list using bedtools. Black list for mm10: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz
bedtools intersect -v -a ${samp_name}_filtered_noDUP_unique.bam -b mm10_ENCODE_blacklist_chip-seq.bed > ${samp_name}_filtered_noDUP_unique.noB_unique.bam

# Indexing BAM file with total reads
samtools index ${samp_name}_filtered_noDUP_unique.noB_unique.bam

# Counting total reads in BAM file
samtools view -c ${samp_name}_filtered_noDUP_unique.noB_unique.bam > ${samp_name}_count_reads_filtered_TOTAL.txt

# Counting reads aligned to the mouse chromosomes in BAM file
samtools view ${samp_name}_filtered_noDUP_unique.noB_unique.bam $(echo chr{1..19} chrX chrY) | wc -l > ${samp_name}_count_reads_filtered_mouse.txt

# Counting reads aligned to the Drosophila chromosomes (spikein reads) in BAM file 
# Note that chromosomes names '4', 'X', 'Y' are from Drosophila and differ from chomosome names 'chr4', 'chrX' and 'chrY' from mouse
samtools view ${samp_name}_filtered_noDUP_unique.noB_unique.bam 2L 2R 3L 3R 4 X Y | wc -l > ${samp_name}_count_reads_filtered_Drosophila.txt

# Subset the BAM file to extract reads uniquely aligned to either the mouse or Drosophila genomes
samtools view -b ${samp_name}_filtered_noDUP_unique.noB_unique.bam $(echo chr{1..19} chrX chrY) > ${samp_name}_filtered_mouse_reads.bam
samtools view -b ${samp_name}_filtered_noDUP_unique.noB_unique.bam 2L 2R 3L 3R 4 X Y > ${samp_name}_filtered_Drosophila_reads.bam

# Make index for BAM file containing reads uniquely aligned to the mouse genome
samtools index ${samp_name}_filtered_mouse_reads_final.bam

# Creating a bigwig file (used for visualization with genome browser) from the BAM file containing reads uniquely aligned to the mouse genome
bamCoverage --bam ${samp_name}_filtered_mouse_reads_final.bam -o ${samp_name}_filtered_mouse_reads_final.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --ignoreDuplicates \
    --ignoreForNormalization chrX chrY \
    --extendReads \
    --centerReads
