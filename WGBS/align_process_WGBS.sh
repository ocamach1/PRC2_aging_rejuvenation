
cd $outdir/$samp_name/Trim_galore_output

########### EXTRA
mkdir fastqc_before_trimming
fastqc -t 8 -o $outdir/$samp_name/Trim_galore_output/fastqc_before_trimming $R1_fastq
fastqc -t 8 -o $outdir/$samp_name/Trim_galore_output/fastqc_before_trimming $R2_fastq
############

echo "Running Trim galore"
trim_galore -j 8 --paired $R1_fastq $R2_fastq -o $outdir/$samp_name/Trim_galore_output/
#trim_galore -j 16 --paired --three_prime_clip_R1 90 $R1_fastq $R2_fastq -o $outdir/$samp_name/Trim_galore_output/

# In case we use SWIFT library:
#trim_galore -j 4  --paired --clip_R1 10 --clip_R2 15 $R1_fastq $R2_fastq -o $outdir/$samp_name/Trim_galore_output/

#echo "Finished trimgalore"

cd $outdir/$samp_name/Bismark_output

echo "Running bismark alignment"

trimmed_R1="$(echo $R1_fastq | cut -d '/' -f8)"
trimmed_R2="$(echo $R2_fastq | cut -d '/' -f8)"

echo "Trimmed R1 file is ${trimmed_R1}"
echo "Trimmed R2 file is ${trimmed_R2}"

prefix_R1="$(echo $trimmed_R1 | cut -d '.' -f1)"
prefix_R2="$(echo $trimmed_R2 | cut -d '.' -f1)"

echo "Trimmed R1 prefix is ${prefix_R1}"
echo "Trimmed R2 prefix is ${prefix_R2}"
 
###########EXTRA
#cd $outdir/$samp_name/Trim_galore_output
#mkdir fastqc_after_trimming
#fastqc -t 4 -o $outdir/$samp_name/Trim_galore_output/fastqc_after_trimming $outdir/$samp_name/Trim_galore_output/${prefix_R1}*.fq.gz
#fastqc -t 4 -o $outdir/$samp_name/Trim_galore_output/fastqc_after_trimming $outdir/$samp_name/Trim_galore_output/${prefix_R2}*.fq.gz
###########

# CHANGE THE GENOME HERE

bismark --bowtie2 --bam --parallel 8 --output_dir $outdir/$samp_name/Bismark_output --temp_dir $outdir/$samp_name/Bismark_temp /dcs07/afeinber/data/personal/ocamacho/mm10_and_EM-Seq_control -1 $outdir/$samp_name/Trim_galore_output/${prefix_R1}*.fq.gz -2 $outdir/$samp_name/Trim_galore_output/${prefix_R2}*.fq.gz
