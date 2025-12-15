# Required packages:
# - EDD

# EDD requires:
- Filtered BAM files for both IP (DNA bound to the target of interest) and Input (total chromatin DNA that serves as background control) samples. 
- Configuration file: available from https://github.com/CollasLab/edd including the default parameters. Default parameters were used in this study.
- mm10 chromosome sizes:
	- Available by using the command 'fetchChromSizes mm10 > mm10.sizes' from the the UCSC Genome Browser utilities software 
	(https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads).

edd mm10.sizes --config-file EDD.config.txt --write-log-ratios --write-bin-scores \
${samp_name_ChIP}_filtered_mouse_reads_final_nsorted.bam ${samp_name_INPUT}_filtered_mouse_reads_final_nsorted.bam
