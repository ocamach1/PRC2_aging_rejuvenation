module load samtools

##sicer2 requires bed files for both ChIP and Input samples
#Sort bam files by name
samtools sort -n -o ${samp_name}_filtered_mouse_reads_final_nsorted.bam $outdir ${samp_name}_filtered_mouse_reads_final.bam
#Convert bam to bed
bedtools bamtobed -i ${samp_name}_filtered_mouse_reads_final_nsorted.bam  > ${samp_name}_filtered_mouse_reads_final_nsorted.bed

# Run sicer2 specifying ChIP and Input samples
sicer -t ${samp_name_ChIP}_filtered_mouse_reads_final_nsorted.bed \
        -c ${samp_name_Input}_filtered_mouse_reads_final_nsorted.bed \
        --species mm10 \
	-fdr 0.05 \
	-w 200 \
	-g 600
