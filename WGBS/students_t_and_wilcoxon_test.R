library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(grid)
library(Matrix)


####################
# Perform Student's t-test
# Every sample is treated as an independent biological replicate
####################

# Import data from bigwig files generated with the informME pipeline
MML_1<-import("MML-Y4F-1.bw")
genome(MML_1)<-"mm10"
colnames(mcols(MML_1))[1]<-"MML"
MML_1$phenoname="Y-1"
MML_1$group="Y"

MML_2<-import("MML-Y4F-3.bw")
genome(MML_2)<-"mm10"
colnames(mcols(MML_2))[1]<-"MML"
MML_2$phenoname="Y-2"
MML_2$group="Y"

MML_3<-import("MML-Y4F-4.bw")
genome(MML_3)<-"mm10"
colnames(mcols(MML_3))[1]<-"MML"
MML_3$phenoname="Y-3"
MML_3$group="Y"

MML_4<-import("MML-5760-con.bw")
genome(MML_4)<-"mm10"
colnames(mcols(MML_4))[1]<-"MML"
MML_4$phenoname="O-1"
MML_4$group="O"

MML_5<-import("MML-5761-con.bw")
genome(MML_5)<-"mm10"
colnames(mcols(MML_5))[1]<-"MML"
MML_5$phenoname="O-2"
MML_5$group="O"

MML_6<-import("MML-5762-con.bw")
genome(MML_6)<-"mm10"
colnames(mcols(MML_6))[1]<-"MML"
MML_6$phenoname="O-3"
MML_6$group="O"

MML_7<-import("MML-5764-con.bw")
genome(MML_7)<-"mm10"
colnames(mcols(MML_7))[1]<-"MML"
MML_7$phenoname="O-4"
MML_7$group="O"

MML_8<-import("MML-253-dox.bw")
genome(MML_8)<-"mm10"
colnames(mcols(MML_8))[1]<-"MML"
MML_8$phenoname="O+OSKM-1"
MML_8$group="O+OSKM"

MML_9<-import("MML-284-dox.bw")
genome(MML_9)<-"mm10"
colnames(mcols(MML_9))[1]<-"MML"
MML_9$phenoname="O+OSKM-2"
MML_9$group="O+OSKM"

MML_10<-import("MML-285-dox.bw")
genome(MML_10)<-"mm10"
colnames(mcols(MML_10))[1]<-"MML"
MML_10$phenoname="O+OSKM-3"
MML_10$group="O+OSKM"

MML_11<-import("MML-288-dox.bw")
genome(MML_11)<-"mm10"
colnames(mcols(MML_11))[1]<-"MML"
MML_11$phenoname="O+OSKM-4"
MML_11$group="O+OSKM"

MML_12<-import("MML-904-dox.bw")
genome(MML_12)<-"mm10"
colnames(mcols(MML_12))[1]<-"MML"
MML_12$phenoname="O+OSKM-5"
MML_12$group="O+OSKM"

# Make combined object for boxplots and density plots
MML = c(MML_1,MML_2,MML_3,MML_4,MML_5,MML_6,MML_7,MML_8,MML_9,MML_10,MML_11,MML_12)
MML$phenoname<-factor(MML$phenoname,levels=c("Y-1","Y-2","Y-3","O-1","O-2","O-3","O-4","O+OSKM-1","O+OSKM-2","O+OSKM-3","O+OSKM-4","O+OSKM-5"))

# Get average MML/NME value for each sample
means_df <- as.data.frame(mcols(MML)) %>%
  group_by(phenoname) %>%
  summarize(mean_MML = mean(MML), .groups = "drop") %>%
  arrange(phenoname)

# Specify the samples that need to be compared
group1_ids <- c("O-1", "O-2", "O-3", "O-4")       # Young
group2_ids <- c("Y-1", "Y-2", "Y-3") # Old

group1_vals <- means_df %>%
  filter(phenoname %in% group1_ids) %>%
  pull(mean_MML)

group2_vals <- means_df %>%
  filter(phenoname %in% group2_ids) %>%
  pull(mean_MML)

# Perform Student's t-test
res <- t.test(group1_vals, group2_vals, var.equal = TRUE)
print(res)



####################
# Perform Wilcoxon signed rank test
####################
# The genome is tiled into 150-bp genomic units, corresponding to the analytical units for which the informME pipeline calculated MML and MML values.
# Per-group averages are computed for each 150 bp region.
# A paired Wilcoxon signed-rank test across the 150 bp regions is performed.

# Make a 150bp "tiling" of the mouse genome over chromosomes 1-19
chrsOfInterest <- paste("chr",1:19,sep="")
tiles <- tileGenome(seqinfo(Mmusculus), tilewidth=150,
                    cut.last.tile.in.chrom=TRUE)
tiles <- tiles[seqnames(tiles) %in% chrsOfInterest]

# Import bigwig files generated with the informME pipeline, specifying the 2 groups needed to be compared
inFolder<-"aging_bw_files"
outFolder<-"wilcoxon_test_results"

# Specify control (young) samples
MML1 <- "MML-Y4F-1"
MML2 <- "MML-Y4F-3"
MML3 <- "MML-Y4F-4"
control <- c(MML1,MML3,MML4)

# Specify treated (old) samples
MML4 <- "MML-5760-con" 
MML5 <- "MML-5761-con"
MML6 <- "MML-5762-con"
MML7 <- "MML-5764-con"
treated <- c(MML4,MML5,MML6,MML7)

# Make function to add data onto tiles
addTrack <- function(tiles,trackFileName,inFolder){
  GR <- import.bw(paste(inFolder,trackFileName,".bw",sep=""))
  mcols(tiles)[[trackFileName]] <- NaN
  olaps <- findOverlaps(GR,tiles)
  mcols(tiles)[[trackFileName]][subjectHits(olaps)] <- GR$score[queryHits(olaps)]
  tiles
}
for (ctrl in controls){
  tiles <- addTrack(tiles,ctrl,inFolder)
}
for (trtd in treated){
  tiles <- addTrack(tiles,trtd,inFolder)
}

# Calculate averages over groups in each 150 bp genomic unit
tiles$ctrlAve <- rowMeans(as.matrix(mcols(tiles)[,controls]),na.rm=TRUE)
tiles$treatAve <- rowMeans(as.matrix(mcols(tiles)[,treated]),na.rm=TRUE)

# Find genomic units with data in both groups and do a (paired) wilcoxon signed rank test
notNAinds <- !(is.na(tiles$ctrlAve)|is.na(tiles$treatAve))

# Perform paired Wilcoxon signed-rank test across the 150 bp regions
testResults <- wilcox.test(tiles$ctrlAve[notNAinds],
                           tiles$treatAve[notNAinds],
                           alternative = "two.sided",paired = TRUE,conf.int = TRUE)
print(testResults)
