# This is an R function that ranks all Mouse genes in the Bioconductor library TxDb.Mmusculus.UCSC.mm10.refGene using 
# Jensen-Shannon distance (JSD) and the average mutual information based on the method described in https://doi.org/10.1186/s12859-019-2777-6. 
# It should be run within an R session.

# The following R libraries must be installed:
#   - GenomicFeatures
#   - GenomicRanges
#   - rtracklayer
#   - TxDb.Mmusculus.UCSC.mm10.refGene


  default usage (replicate reference data is available):

   source("jsGrank.R")
   rankGenes(refVrefFiles,testVrefFiles,inFolder,outFolder,
             tName,rName)

   # refVrefFiles is a vector of BIGWIG files that contain the
   # JSD values of a test/reference comparison. 
   # For example: if
   #
   # JSD-lungnormal-1-VS-lungnormal-2.bw 
   # JSD-lungcancer-3-VS-lungnormal-1.bw 
   # JSD-lungnormal-3-VS-lungnormal-2.bw
   # 
   # are available, then set 
   # 
   # textVrefFiles <- c("JSD-lungnormal-1-VS-lungnormal-2.bw",
   #                    "JSD-lungnormal-3-VS-lungnormal-1.bw",
   #                    "JSD-lungnormal-3-VS-lungnormal-2.bw")
   #
   # testVrefFiles is a vector of BIGWIG files that contain the  
   # JSD values of available test/reference comparisons. 
   # For example: if 
   #
   # JSD-lungcancer-1-VS-lungnormal-1.bw  
   # JSD-lungcancer-2-VS-lungnormal-2.bw 
   # JSD-lungcancer-3-VS-lungnormal-3.bw 
   # 
   # are available, then set 
   # 
   # textVrefFiles <- c("JSD-lungcancer-1-VS-lungnormal-1.bw",
   #                    "JSD-lungcancer-2-VS-lungnormal-2.bw",
   #                    "JSD-lungcancer-3-VS-lungnormal-3.bw")
   #
   # inFolder is the directory that contains the JSD files
   # outFolder is the directory used to write the result  
   # (a .tsv file).
   # 
   # For example:
   # 
   # inFolder  <- "/path/to/in-folder/"
   # outFolder <- "/path/to/out-folder/"
   #
   # tName and rName are strings providing names for the 
   # test and reference phenotypes.
   #
   # For example: 
   #
   # tName <- "lungcancer"
   # rName <- "lungnormal"

  default usage (no replicate reference data is available):  

   source("jsGrank.R")
   rankGenes(c(),testVrefFiles,inFolder,outFolder,
             tName,rName)

   # testVrefFiles is a vector of BIGWIG files that contain the  
   # JSD values of available test/reference comparisons. 
   # For example: if 
   #
   # JSD-lungcancer-1-VS-lungnormal-1.bw 
   # JSD-lungcancer-2-VS-lungnormal-2.bw 
   # JSD-lungcancer-3-VS-lungnormal-3.bw 
   # 
   # are available, then set 
   # 
   # textVrefFiles <- c("JSD-lungcancer-1-VS-lungnormal-1.bw",
   #                    "JSD-lungcancer-2-VS-lungnormal-2.bw",
   #                    "JSD-lungcancer-3-VS-lungnormal-3.bw")
   #
   # inFolder is the directory that contains the JSD files
   # outFolder is the directory used to write the result 
   # (a .tsv file).
   # 
   # For example:
   # 
   # inFolder  <- "/path/to/in-folder/"
   # outFolder <- "/path/to/out-folder/"
   #
   # tName and rName are strings providing names for the 
   # test and reference phenotypes.
   #
   # For example: 
   #
   # tName <- "lungcancer"
   # rName <- "lungnormal"
   
  requirements:


