library(RnBeads)
library(DecompPipeline)


# Importing the data.
# Upload data specifying the sample annotation file ("skin_aging_sample_annotation.tsv") as input, which contains the bedmethyl 
# (Bismark cov) files and the corresponding sample names included in the analysis. In this analysis, the content of the annotation file is:

#sampleID,filename
#Y4F-1,Y4F-1.bismark.cov
#Y4F-3,Y4F-3.bismark.cov
#Y4F-4,Y4F-4.bismark.cov
#5760-con,5760-con.bismark.cov
#5761-con,5761-con.bismark.cov
#5762-con,5762-con.bismark.cov
#5764-con,5764-con.bismark.cov
#253-dox,253-dox.bismark.cov
#284-dox,284-dox.bismark.cov
#285-dox,285-dox.bismark.cov
#288-dox,288-dox.bismark.cov
#904-dox,904-dox.bismark.cov
#O-005,O-005.bismark.cov
#O-989,O-989.bismark.cov
#O-990,O-990.bismark.cov
#O-991,O-991.bismark.cov
#O-992,O-992.bismark.cov
#Y-458,Y-458.bismark.cov
#Y-459,Y-459.bismark.cov
#Y-461,Y-461.bismark.cov
#Y-463,Y-463.bismark.cov
#Y-464,Y-464.bismark.cov

rnb.options(
  assembly = "mm10",
  identifiers.column = "submitter_id",
  import = TRUE,
  import.default.data.type = "bed.dir",
  import.table.separator = "\t",
  import.sex.prediction = TRUE,
  qc = TRUE,
  preprocessing = FALSE,
  exploratory = FALSE,
  inference = FALSE,
  differential = FALSE,
  export.to.bed = FALSE,
  export.to.trackhub = NULL,
  export.to.csv = FALSE
)

sample.anno <- "annotation/skin_aging_sample_annotation.tsv"
idat.folder <- "idat/"
dir.report <- paste0("report",Sys.Date(),"/")
temp.dir <- "/tmp"
options(fftempdir = temp.dir)
rnb.set <- rnb.run.analysis(
  dir.reports = dir.report,
  sample.sheet = sample.anno,
  data.dir = idat.folder)

# Preprocessing and filtering
# All sites that have a read coverage lower than 5 in any of the samples are removed.
# High and low coverage outliers, sites with missing values, annotated SNPs, and sites on the sex chromosomes are also removed.
data.prep <- prepare.data(rnb.set = rnb.set,
                             analysis.name = "skin_aging",
                             filter.coverage = TRUE,
                             min.coverage = 5,
                             min.covg.quant = 0.001,
                             max.covg.quant = 0.999,
                             filter.na = TRUE,
                             filter.snp = TRUE,
                             filter.sex.chromosomes = TRUE,
                             execute.lump = TRUE)


# Selection of CpG subsets
# We select the 5000 most variably methylated CpGs across the samples for downstream analysis.
cg_subset <- prepare.CG.subsets(rnb.set=data.prep$rnb.set.filtered,
                                marker.selection = "var",
                                n.markers = 5000)

# Methylome deconvolution
# We use MeDeCom with a grid of values for the number of components (K) ranging from 2 to 15, which covers homogeneous to heterogeneous
# samples. We also specify a grid for the regularization parameter (λ) from strong (0.01) to no regularization (0).
md.res <- start.medecom.analysis(
  rnb.set = data.prep$rnb.set.filtered,
  cg.groups = cg_subset,
  Ks = 2:15,
  lambda.grid = c(0,0.01,0.001,0.0001,0.00001),
  factorviz.outputs = TRUE,
  analysis.name = "skin_aging",
  cores = 15)

# Interpretation of the results is done with the FactorViz to visualize and ineractively explore the deconvolution results.









