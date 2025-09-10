library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggplot2)#for plotting
library(grid)
library(ggpubr)
library(rstatix)
library(svglite)
library(rtracklayer)


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


######################
######################
# Plot average MML and NME levels on Ezh2 binding sites in mESC (as defined in available ChIP-Seq data downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4774490)
# Boxplots represent the distribution of the average values of MML/NME of each sample

ezh2_peaks<-import("Ezh2_binding_sites_mESC.bed")
MML$chrom_state = "blank"
olaps <- findOverlaps(ezh2_peaks,MML)
MML$chrom_state[subjectHits(olaps)]="marked"
MML_subset <- MML[MML$chrom_state == "marked"]

# Compute average MML/NME per sample
mean_per_sample <- tapply(
  MML_subset$MML,
  MML_subset$phenoname,
  mean,
  na.rm = TRUE
)

df <- data.frame(
  Sample = names(mean_per_sample),
  Value = as.numeric(mean_per_sample),
  stringsAsFactors = FALSE
)
sample_to_group <- data.frame(
  Sample = names(mean_per_sample),
  Group = c(rep("Young", 3), rep("Old", 4), rep("Old+OSKM", 5))
)

df <- merge(df, sample_to_group, by = "Sample")
df$Group <- factor(df$Group, levels = c("Young", "Old", "Old+OSKM"))
manPalette <- c("#619CFF", "#F8766D", "#00BA38")

p <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.35, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(y = "Average MML", x = "") +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0.11, 0.15)) +
  scale_y_continuous(breaks = seq(0.11, 0.15, by = 0.01)) +
  scale_fill_manual(values = manPalette)
ggsave(p, file="average_MML_average_boxplots_Ezh2.svg", width=5.2, height=9.0,path="plots/")


######################
######################
# Plot average MML and NME levels on H3K27me3-marked regions in mESC (as defined in available ChIP-Seq data downloaded from https://www.encodeproject.org/experiments/ENCSR000CFN/)
# Boxplots represent the distribution of the average values of MML/NME of each sample

H3K27me3_peaks<-import("H3K27me3_peaks_mESC.bed")
MML$chrom_state = "blank"
olaps <- findOverlaps(H3K27me3_peaks,MML)
MML$chrom_state[subjectHits(olaps)]="marked"
MML_subset <- MML[MML$chrom_state == "marked"]

# Compute average MML/NME per sample
mean_per_sample <- tapply(
  MML_subset$MML,
  MML_subset$phenoname,
  mean,
  na.rm = TRUE
)

df <- data.frame(
  Sample = names(mean_per_sample),
  Value = as.numeric(mean_per_sample),
  stringsAsFactors = FALSE
)
sample_to_group <- data.frame(
  Sample = names(mean_per_sample),
  Group = c(rep("Young", 3), rep("Old", 4), rep("Old+OSKM", 5))
)

df <- merge(df, sample_to_group, by = "Sample")
df$Group <- factor(df$Group, levels = c("Young", "Old", "Old+OSKM"))
manPalette <- c("#619CFF", "#F8766D", "#00BA38")

p <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.35, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(y = "Average MML", x = "") +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0.11, 0.15)) +
  scale_y_continuous(breaks = seq(0.11, 0.15, by = 0.01)) +
  scale_fill_manual(values = manPalette)
ggsave(p, file="average_MML_average_boxplots_H3K27me3.svg", width=5.2, height=9.0,path="plots/")




 


