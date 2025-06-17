# Required packages:
# - samtools

# Sicer2 requires bed files for both IP and Input samples. bed files are obtained from the previously processed bam files containing reads uniquely aligned to the mouse genome.
# Sorting bam files by read name
samtools sort -n -o ${samp_name}_filtered_mouse_reads_final_nsorted.bam ${samp_name}_filtered_mouse_reads_final.bam
# Converting bam to bed format
bedtools bamtobed -i ${samp_name}_filtered_mouse_reads_final_nsorted.bam  > ${samp_name}_filtered_mouse_reads_final_nsorted.bed

# Running Sicer2
# -t expects the bed file containing the data from the IP sample
# -c expects the bed file containing the data from the Input sample
sicer -t ${samp_name_ChIP}_filtered_mouse_reads_final_nsorted.bed \
        -c ${samp_name_Input}_filtered_mouse_reads_final_nsorted.bed \
        --species mm10 \
	-fdr 0.05 \
	-w 200 \
	-g 600
