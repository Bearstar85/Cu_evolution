#Houskeeping the loading packages####
#Lets clear old objects
rm(list=ls())
getwd()
#and set the work directory in case we have moved around or opened another project
setwd("~/Documents/R/Cu_evolution/Barcodes/Barcoding_C12W1/C12W1_unknowns")
dir()

#Bash commands for formation input data
#cd /Users/xanbjg/Documents/R/Cu_evolution/Barcodes/C12W1_unknowns
#rm *.Unknown_*
#for i in *.Unknown_Sequences*; do mv "$i" "${i/.Unknown_Sequences/}"; done

#Load packages and data (not adapted)
library(reshape2)
library(multtest) 
library(lattice)
library(devtools)
library(ggbiplot)
library(plyr)
library(scales)
library(ggplot2)
library(devtools)
library(grid)
#library(gplots)
library(tidyverse)
library(lubridate) # useful for working with dates
library(cowplot) # useul for combining multiple plots
library(ggthemes)
library(broom)
library(ggpubr)
library(gridExtra)
#library(ggpmisc)
library(gplots)

#Process data#####

#Make a vector list of all input file names
InputFiles <- list.files(path = "~/Documents/R/Cu_evolution/Barcodes/Barcoding_C12W1/C12W1_unknowns/Input")
InputFiles
#And read in the get Indexing file
#Allels <- read.delim("../Allele_indexing.txt", header=TRUE)
Indexing <- read.delim("../Indexing.txt", header=TRUE)

#Indexing and summarizing the input files (note that 1_header.tmp needs to be created by loop below before this works, if its not been run before run loop once and come back here)

for (i in InputFiles) {
  #Import data (too be looped)
  input=paste("Input/",i,"", sep="")
  myData <- read.table(file = input, sep = '\t', header = TRUE)
  sapply(myData, class)
  
  
  #Generate sample ID for indexing experimental samples
  
  Name <- i
  ID <- gsub('.tsv', '', Name)
  N <- nrow(myData)
  Sample <- rep(ID, each = N)
  Sample <- as.data.frame(Sample)
  colnames(Sample) <- c("ID")
  Index <- join(Sample, Indexing, by = "ID")
  
  #Fuse the sample indexing with data
  myData1 <- cbind(Index, myData)
  
  #And count reads
  Nreads <- sum(myData1$Quantity)
  Read_counts <- cbind(ID, Nreads)
  
  #Well save the header and export data
  output <- paste("Indexed/",ID,".csv", sep="")
  output2 <- paste("Read_counts/",ID,".csv", sep="")
  output3 <- paste("Plots/",ID,".pdf", sep="")
  
  header <- colnames(myData1)
  
  write.table(header, "Indexed/1_header.tmp", sep=",",  col.names=FALSE, row.names =FALSE)
  write.table(myData1, output, sep=",",  col.names=FALSE, row.names =FALSE)
  write.table(Read_counts, output2, sep=",",  col.names=FALSE, row.names =FALSE)
  
}

#Read_count_summary#####

#This section requires that the below files have been made by 
#Okey lets merge the count files so we see where were at
Unknown_counts <- read.csv("Read_counts/Unknown_counts.txt", header=FALSE)
Total_hits <- read.csv("../C12W1_abundances/Read_counts/ReadCount.txt", header=FALSE)
#VG_hits <- read.csv("../C12W1_abundances/Read_counts/VG_adapted/ReadCount.txt", header=FALSE)
#Total_hits <- rbind(GP_hits, VG_hits)
colnames(Total_hits) <- c("ID", "Hits")
colnames(Unknown_counts) <- c("ID", "Unknown")
Hits_vs_Unknowns <- join(Total_hits, Unknown_counts, by = "ID")

#Compute fraction hits
Hits_vs_Unknowns$Fraction <- Hits_vs_Unknowns$Hits/(Hits_vs_Unknowns$Hits+Hits_vs_Unknowns$Unknown)
Hits_vs_Unknowns$Total <- (Hits_vs_Unknowns$Hits+Hits_vs_Unknowns$Unknown)

#and amend the raw read counts per sample
RawReads <- read.table(file = "RawReads.txt", sep = '\t', header = TRUE)
Hits_vs_Unknowns2 <- join(Hits_vs_Unknowns, RawReads, by = "ID")
Hits_vs_Unknowns2$Fraction2 <- Hits_vs_Unknowns2$Hits/(Hits_vs_Unknowns2$Raw.Reads)
head(Hits_vs_Unknowns2)

#And also add all the mapped reads counts including allels that where not uniqe
All_allels_total <- read.csv(file = "../C12W1_abundances/Read_counts/Counts_all_knownMapped.txt", header = FALSE)
colnames(All_allels_total) <- c("ID", "All_allels_total")
Hits_vs_Unknowns3 <- join(All_allels_total, Hits_vs_Unknowns2, by = "ID")
Hits_vs_Unknowns3$NonUniqe_Fraction <- (Hits_vs_Unknowns3$All_allels_total-Hits_vs_Unknowns3$Hits)/Hits_vs_Unknowns3$All_allels_total

#and save results
write.table(Hits_vs_Unknowns3, "Results/Hits_vs_Unknowns.txt", sep='\t',  col.names=TRUE, row.names =FALSE)
unique(Hits_vs_Unknowns3$ID)
mean(Hits_vs_Unknowns3$Fraction)
max(Hits_vs_Unknowns3$Fraction)
min(Hits_vs_Unknowns3$Fraction)

mean(Hits_vs_Unknowns3$Fraction2)
max(Hits_vs_Unknowns3$Fraction2)
min(Hits_vs_Unknowns3$Fraction2)

sapply(Hits_vs_Unknowns3, class)

mean(Hits_vs_Unknowns3$Hits)
mean(Hits_vs_Unknowns3$Unknown)

Fig1 <- ggplot(Hits_vs_Unknowns, stats = "identity", aes(Fraction)) +
  geom_col(data=Hits_vs_Unknowns, aes(ID, Fraction)) +
  #geom_hline(yintercept = 0.02) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = "Fraction mapped Reads of assembeled amplicons") +
  #coord_cartesian(ylim=c(0.0001, 0.01), expand = F) + #ylim=c(-0, 12.5)
  #scale_y_discrete(limits=c(0.0001,0.001,0.01,0.1,1)) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #labels = trans_format("log10", math_format(10^.x))) +
  #coord_flip() +
  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=5))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  #stat_cor(aes(color = scientific_name), label.x = 3)
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(legend.position ="top",
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #scale_fill_manual(breaks=c(3)),
        legend.margin = margin(2, 2, 2, 2),
        #guides(color = FALSE), #,color = FALSE, size = FALSE,
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig1)

dev.copy(pdf, "Plots/MappedReads.pdf")
dev.off()

Fig2 <- ggplot(Hits_vs_Unknowns, stats = "identity", aes(Total)) +
  geom_col(data=Hits_vs_Unknowns, aes(ID, Total)) +
  #geom_hline(yintercept = 0.02) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = "assembeled amplicons") +
  #coord_cartesian(ylim=c(0.0001, 0.01), expand = F) + #ylim=c(-0, 12.5)
  #scale_y_discrete(limits=c(0.0001,0.001,0.01,0.1,1)) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #labels = trans_format("log10", math_format(10^.x))) +
  #coord_flip() +
  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=5))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  #stat_cor(aes(color = scientific_name), label.x = 3)
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(legend.position ="top",
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #scale_fill_manual(breaks=c(3)),
        legend.margin = margin(2, 2, 2, 2),
        #guides(color = FALSE), #,color = FALSE, size = FALSE,
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig2)

dev.copy(pdf, "Plots/TotalReads.pdf")
dev.off()

Fig3 <- ggplot(Hits_vs_Unknowns3, stats = "identity", aes(Raw.Reads)) +
  geom_col(data=Hits_vs_Unknowns3, aes(ID, Raw.Reads)) +
  #geom_hline(yintercept = 0.02) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = "Raw reads all PIPT C12W1") +
  #coord_cartesian(ylim=c(0.0001, 0.01), expand = F) + #ylim=c(-0, 12.5)
  #scale_y_discrete(limits=c(0.0001,0.001,0.01,0.1,1)) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #labels = trans_format("log10", math_format(10^.x))) +
  #coord_flip() +
  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=5))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  #stat_cor(aes(color = scientific_name), label.x = 3)
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(legend.position ="top",
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #scale_fill_manual(breaks=c(3)),
        legend.margin = margin(2, 2, 2, 2),
        #guides(color = FALSE), #,color = FALSE, size = FALSE,
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig3)

dev.copy(pdf, "Plots/RawReads.pdf")
dev.off()

Fig4 <- ggplot(Hits_vs_Unknowns3, stats = "identity", aes(Fraction2)) +
  geom_col(data=Hits_vs_Unknowns3, aes(ID, Fraction2)) +
  #geom_hline(yintercept = 0.02) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = "Mapped fraction Raw reads PIPT C12W1") +
  #coord_cartesian(ylim=c(0.0001, 0.01), expand = F) + #ylim=c(-0, 12.5)
  #scale_y_discrete(limits=c(0.0001,0.001,0.01,0.1,1)) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #labels = trans_format("log10", math_format(10^.x))) +
  #coord_flip() +
  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=5))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  #stat_cor(aes(color = scientific_name), label.x = 3)
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(legend.position ="top",
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #scale_fill_manual(breaks=c(3)),
        legend.margin = margin(2, 2, 2, 2),
        #guides(color = FALSE), #,color = FALSE, size = FALSE,
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig4)

dev.copy(pdf, "Plots/FractionRawReadsUsefull.pdf")
dev.off()

Fig5 <- ggplot(Hits_vs_Unknowns3, stats = "identity", aes(NonUniqe_Fraction)) +
  geom_col(data=Hits_vs_Unknowns3, aes(ID, NonUniqe_Fraction)) +
  #geom_hline(yintercept = 0.02) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = "Fraction non-uniqe allels filtered out") +
  #coord_cartesian(ylim=c(0.0001, 0.01), expand = F) + #ylim=c(-0, 12.5)
  #scale_y_discrete(limits=c(0.0001,0.001,0.01,0.1,1)) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #labels = trans_format("log10", math_format(10^.x))) +
  #coord_flip() +
  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=5))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  #stat_cor(aes(color = scientific_name), label.x = 3)
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(legend.position ="top",
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #scale_fill_manual(breaks=c(3)),
        legend.margin = margin(2, 2, 2, 2),
        #guides(color = FALSE), #,color = FALSE, size = FALSE,
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig5)

dev.copy(pdf, "Plots/FractionNon-uniqe.pdf")
dev.off()

#Find true unknown sequences#####

#Lets start with the controls and see what reads are present  in all MM samples
#Start with GP
GP_MM1 <- read.table(file = "Input/P21502_250.tsv", sep = '\t', header = TRUE)
GP_MM2 <- read.table(file = "Input/P21502_251.tsv", sep = '\t', header = TRUE)
GP_MM3 <- read.table(file = "Input/P21502_252.tsv", sep = '\t', header = TRUE)

colnames(GP_MM1) <- c("GP_MM1", "Unknown.barcode")
colnames(GP_MM2) <- c("GP_MM2", "Unknown.barcode")
colnames(GP_MM3) <- c("GP_MM3", "Unknown.barcode")


sapply(GP_MM1, class)
GP_MMall <- join(GP_MM1, GP_MM2, by = "Unknown.barcode")
GP_MMall <- join(GP_MMall, GP_MM3, by = "Unknown.barcode")

#Okey lets remove NAs to get reads present in all samples

GP_MMall2 <- na.omit(GP_MMall)
#SEARCH for  a sequence####
DNA <- "CCTCAAACCCCATCGAATACAAACTCTCCATATGCATCAAAAGCTGTGCCGCCTTATCACTCGCCCCCTTTTCCAATGAATTGGCATAGGCATTAATAATCGCACCAAATGTCAATGCGTCGGGTAGTAAATCATCATCCCCCTCATTGTACAGACGGCCCATTTCATCTAGTAAAAACTTGGCACGATCGGCGCAATTCTTGTCACGGGAACGTGCTACTGCGTTGATTACCGCGTTGAAGGATCGTGCATCGGGTTTGCAGTTTATGTCTAGTCCACTATCGTTATTGTAAAAGTCGTACATGTTACGGAGCAATTCTTCTGCTTTACGTGCCGATCCTCGTTCATTGCTACGTGCCCATGCCGTGATGACGGAGTTGTAACTAATTTTGTTTGGTTTGACATGATCCCATTTTCCATCGATATTATGTTTGGGAGAGGAGTAGAGGTGCAGCATTTGATCGAGGAGGGCTTGTGCGTGCTTTGCGGAATCTTTGCGTGTGCATCGTGCC"
DNApresenet <- subset.data.frame(GP_MMall2, grepl(print(DNA), GP_MMall2$Unknown.barcode))
DNApresenet

#Lets see which are most numerous and least variable

GP_MMall2$mean <- (GP_MMall2$GP_MM1+GP_MMall2$GP_MM2+GP_MMall2$GP_MM3)/3
plot(GP_MMall2$mean)
dev.copy(pdf, "Plots/GP_Unknows_plot.pdf")
dev.off()

#Computes SD value for replicates
GP_MMall2 <- GP_MMall2 %>% 
  mutate(ID = row_number()) %>%
  group_by(ID) %>%
  do(data.frame(., SD = sd(unlist(.[c("GP_MM1", "GP_MM3")]), na.rm=T)))

GP_MMall2$CV <- (GP_MMall2$SD/GP_MMall2$mean)

#I want to contrast with what is normal in true allels of these samples
GP_MM1_allels <- read.table(file = "../C12W1_abundances/Indexed/P21502_250.csv", sep = ',', header = F)
GP_MM2_allels <- read.table(file = "../C12W1_abundances/Indexed/P21502_251.csv", sep = ',', header = F)
GP_MM3_allels <- read.table(file = "../C12W1_abundances/Indexed/P21502_252.csv", sep = ',', header = F)

head(GP_MM1_allels)
header <- read.table(file = "../C12W1_abundances/Indexed/1_header.tmp", sep = ',', header = F)
header <- t(header)
GP_MMall_allels <- rbind(GP_MM1_allels, GP_MM2_allels, GP_MM3_allels, GP_MM3_allels)
colnames(GP_MMall_allels) <- (header)
GP_MMall_allels <- subset.data.frame(GP_MMall_allels, grepl("GP", GP_MMall_allels$Strain))
Allel_count_GP <- mean(GP_MMall_allels$Total)

Allel_count_GP <- sd(GP_MMall_allels$Total)
hist(GP_MMall_allels$Total, xlim=c(0,500), breaks=80, main="GP x 3 MM")
dev.copy(pdf, "Plots/GP_Readcountallels_absolute.pdf")
dev.off()


#So most alleles are present at >20 reads (only 0 are not)

Allel_count_GP
plot(GP_MMall2$mean~GP_MMall2$CV, log="xy")
abline(h=20, col="green")
log="xy"
scale_y_log10()
dev.copy(pdf, "Plots/GP_Unknows_plot.pdf")
dev.off()

#Lets do same for VG
VG_MM1 <- read.table(file = "Input/P21502_272.tsv", sep = '\t', header = TRUE)
VG_MM2 <- read.table(file = "Input/P21502_273.tsv", sep = '\t', header = TRUE)
VG_MM3 <- read.table(file = "Input/P21502_274.tsv", sep = '\t', header = TRUE)
VG_MM4 <- read.table(file = "Input/P21502_275.tsv", sep = '\t', header = TRUE)

colnames(VG_MM1) <- c("VG_MM1", "Unknown.barcode")
colnames(VG_MM2) <- c("VG_MM2", "Unknown.barcode")
colnames(VG_MM3) <- c("VG_MM3", "Unknown.barcode")
colnames(VG_MM4) <- c("VG_MM4", "Unknown.barcode")

sapply(VG_MM1, class)
VG_MMall <- join(VG_MM1, VG_MM2, by = "Unknown.barcode")
VG_MMall <- join(VG_MMall, VG_MM3, by = "Unknown.barcode")
VG_MMall <- join(VG_MMall, VG_MM4, by = "Unknown.barcode")

#Okey lets remove NAs to get reads present in all samples

VG_MMall2 <- na.omit(VG_MMall)

#Lets see which are most numerous and least variable

VG_MMall2$mean <- (VG_MMall2$VG_MM1+VG_MMall2$VG_MM2+VG_MMall2$VG_MM3+VG_MMall2$VG_MM4)/4
plot(VG_MMall2$mean)
dev.copy(pdf, "Plots/VG_Unknows_plot.pdf")
dev.off()

sapply(VG_MMall2, class)
VG_MMall2$VG_MM1 <- as.numeric(as.integer(VG_MMall2$VG_MM1))
VG_MMall2$VG_MM2 <- as.numeric(as.integer(VG_MMall2$VG_MM2))
VG_MMall2$VG_MM3 <- as.numeric(as.integer(VG_MMall2$VG_MM3))
VG_MMall2$VG_MM4 <- as.numeric(as.integer(VG_MMall2$VG_MM4))
colnames(Unknown_counts) <- c("ID", "Unknown")

#Computes SD value for replicates
VG_MMall2 <- VG_MMall2 %>% 
  mutate(ID = row_number()) %>%
  group_by(ID) %>%
  do(data.frame(., SD = sd(unlist(.[c("VG_MM1", "VG_MM3", "VG_MM4")]), na.rm=T)))

VG_MMall2$CV <- (VG_MMall2$SD/VG_MMall2$mean)

#I want to contrast with what is normal in true allels of these samples
VG_MM1_allels <- read.table(file = "../C12W1_abundances/Indexed/P21502_272.csv", sep = ',', header = F)
VG_MM2_allels <- read.table(file = "../C12W1_abundances/Indexed/P21502_273.csv", sep = ',', header = F)
VG_MM3_allels <- read.table(file = "../C12W1_abundances/Indexed/P21502_274.csv", sep = ',', header = F)
VG_MM4_allels <- read.table(file = "../C12W1_abundances/Indexed/P21502_275.csv", sep = ',', header = F)
head(VG_MM1_allels)
header <- read.table(file = "../C12W1_abundances/Indexed/1_header.csv", sep = ',', header = F)
header <- t(header)
VG_MMall_allels <- rbind(VG_MM1_allels, VG_MM2_allels, VG_MM3_allels, VG_MM3_allels)
colnames(VG_MMall_allels) <- (header)
VG_MMall_allels <- subset.data.frame(VG_MMall_allels, grepl("VG", VG_MMall_allels$Strain))
Allel_count_VG <- mean(VG_MMall_allels$Total)

Allel_count_VG <- sd(VG_MMall_allels$Total)
hist(VG_MMall_allels$Total, xlim=c(0,500), breaks=100, main="VG x 4 MM")
dev.copy(pdf, "Plots/VG_Readcountallels_absolute.pdf")
dev.off()

#So most alleles are present at >20 reads here (only 7 out of  244 are not). 

Allel_count_VG
plot(VG_MMall2$mean~VG_MMall2$CV, log="xy")
abline(h=20, col="green")
log="xy"
scale_y_log10()
dev.copy(pdf, "Plots/VG_Unknows_plot.pdf")
dev.off()

##Lets push out unknown that are above 20 reads per sample but not annotated####
unknowns <- rbind(GP_MMall2, VG_MMall2)
Unknowns2 <- unknowns %>% filter(mean > 10)

#Are they all uniqe?
#No 25 overlapp
nrow(unknowns)
df_uniq <- unique(unknowns$Unknown.barcode)
length(df_uniq)

#But none of the high abundance ones do
nrow(Unknowns2)
df_uniq <- unique(Unknowns2$Unknown.barcode)
length(df_uniq)

write.table(Unknowns2, "Results/unknow_allels.txt", sep='\t',  col.names=TRUE, row.names =FALSE)

#High abundant reads across all samples####

High_RA <- data.frame(Unknowns2$Unknown.barcode)
colnames(High_RA) <- c("Unknown.barcode")
#i <- "P21502_109.tsv"

#This takes about 20 min to run, OPTIONAL to skip and READ in file below.
for (i in InputFiles) {
  #Import data (too be looped)
  input=paste("Input/",i,"", sep="")
  myData <- read.table(file = input, sep = '\t', header = TRUE)
  ID <- gsub(".tsv", "", i)
  colnames(myData) <- c(ID, "Unknown.barcode")
 #Normalization to RA
  Nreads <- sum(myData[,1])
  myData[,1] <- (myData[,1]/Nreads)


  High_RA <- join(High_RA, myData, by = "Unknown.barcode", type = "full")
}

#Write out the results.
write.table(High_RA, "Results/Summary_all_unknow_allels_RA.txt", sep='\t',  col.names=TRUE, row.names =FALSE)
High_RA2 <- High_RA

#Advisable to just read in this file instead of repeating above loop.
#Results/Summary_all_unknow_allels.txt contains raw counts
#Results/Summary_all_unknow_allels_RA.txt contains Relative abundaces normalized per file (i.e. fraction of unknow amplicons not total)
High_RA2 <- read.table(file = "Results/Summary_all_unknow_allels_RA.txt", sep = '\t', header = T)

sapply(High_RA2, class)
dim(High_RA2)
#sum(High_RA2)

#Lets add up the RA to get a sense of which amplicons are present at high concentrations in multiple samples
High_RA_means <- High_RA2 %>% 
  mutate(ID = row_number()) %>%
  group_by(ID) %>%
  do(data.frame(., Sum = sum(unlist(.[c(2:102)]), na.rm=T)))

#Maybe counting the number of observations is a better filter? (it is not)
#High_RA_means <- High_RA_means %>% 
 # mutate(ID) %>%
  #group_by(ID) %>%
  #do(data.frame(., N = sum(!is.na(.[c(2:107)]), na.rm=T)))

#GP_MMall_allels <- subset.data.frame(GP_MMall_allels, grepl("GP", GP_MMall_allels$Strain))
#Allel_count_GP <- mean(GP_MMall_allels$Total)

#So we have quite a few amplicons showing up at >100 observations, fewer over 500, and max 5000.
hist(High_RA_means$Sum, xlim=c(0,0.4), ylim = c(0,100), breaks=80, main="Sum RA unknown")
dev.copy(pdf, "Plots/All_unknowns_total_observations.pdf")
dev.off()

#In how many samples are amplicons encountered?
#any true allels would be expected in at least 10 samples (MM, genotype, t9)
#alleles associated with positive selection in a treatment in >15.
#This is probablt not a good filter, since its prone to pick up reads associated with common strains
#hist(High_RA_means$N, xlim=c(0,107), ylim = c(0,100), breaks=107, main="N RA unknown")
#dev.copy(pdf, "Plots/All_unknowns_N_observations.pdf")
#dev.off()

#Lets remove everything below 100 (or simillar number for the RA, this get abput 50% reads overrepresented in MM, which is good)
High_RA_means100 <- High_RA_means %>% filter(Sum > 0.025)
High_RA_means500 <- High_RA_means %>% filter(Sum > 0.1)
High_RA_means_less100 <- High_RA_means %>% filter(Sum < 100)

#lets quantify how many observations we have in each category

#hahah ca. 175,000 amplicons are only seen once (only works on raw reads)
hist(High_RA_means_less100$Sum, xlim=c(0,1), breaks=100, main="Sum low abundance unknown")
#dev.copy(pdf, "Plots/All_unknowns_less100observations.pdf")
#dev.off()

#903 amplicons seen more than 100 times (combined 335787)
#sum(High_RA_means100$Sum)

#154 amplicons (76 accounting for primer variants) seen more than 500 time (combined 187170)
#sum(High_RA_means500$Sum)

#38847 amplicons seen less than 100 times (combined 758662)
#sum(High_RA_means_less100$Sum)

#Now lets focus in the >500 group. Where are these reads?
#How many observations do we have?
High_RA_means500 <- High_RA_means500 %>% 
  mutate(ID) %>%
  group_by(ID) %>%
  do(data.frame(., N = sum(!is.na(.[c(2:102)]), na.rm=T)))

#Lets index the samples

#IDs <- unique(High_RA_means500$ID)
#IDs <- as.character(as.integer(IDs))
Indexed500 <- data.frame(ID=character(),
                         Index=character(),
                           Total=numeric(),
                           stringsAsFactors=FALSE)

#Lets change structure so we can index samples and plot them.
for (i in InputFiles) {
  Name <- gsub(".tsv", "", i) 
  N <- nrow(High_RA_means500)
  ID <- rep(Name, each = N)
  Index <- High_RA_means500$"ID"
  Total <- High_RA_means500[grep(Name, names(High_RA_means500))]
  colnames(Total) <- c("Total")
  Column <- cbind.data.frame(ID, Index, Total)
  Indexed500 <- rbind(Indexed500, Column)
}

Indexing <- read.delim("../Indexing.txt", header=TRUE)
Indexed500_2 <- join(Indexed500, Indexing, by = "ID")

#Lets plot them
df_uniq <- unique(Indexed500_2$Index)
colourCount = length(unique(Indexed500_2$Replicates))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
df_uniq
Indexed500_2$Replicate <- as.character(as.integer(Indexed500_2$Replicate))
#i <- "6"

for (i in df_uniq) {
  OneAmplicon <- subset.data.frame(Indexed500_2, grepl(paste("\\b",i,"\\b", sep = ""), Indexed500_2$Index))
  StrainAmplicon <- subset.data.frame(OneAmplicon, grepl('Strain', OneAmplicon$Experiment))
  ExperimentAmplicon <- subset.data.frame(OneAmplicon, grepl(pattern = 'GP|VG', OneAmplicon$Experiment))
  output <- paste("Plots/HighAbundance/",i,".pdf", sep="")
  output2 <- paste("Plots/HighAbundance/",i,"_Strains.pdf", sep="")
  
  
  #This needs to be here for some annoying reason, loop terminates otherwise
  #quartz()
  
  suppressWarnings(FigTotal <- ggplot(ExperimentAmplicon, aes(Timepoint, Total), message=FALSE) +
                     geom_point(mapping = aes(color = Replicate), size = 2, shape = 1) + #shape = Barcode,
                     facet_wrap(~Experiment + Treatment) +
                     theme(strip.text.x = element_text(size=12, angle=0, face = "italic"),
                           strip.background = element_rect(colour="white", fill="white")) +
                     geom_path(mapping = aes(color = Replicate), size = 0.5) +
                     #scale_color_brewer(palette = "Paired") +
                     #geom_errorbar(aes(xmin = (Low), xmax = (High), width=.4), position = position_dodge(0)) +
                     #geom_vline(xintercept = 0.86, color="red") +
                     #geom_vline(xintercept = 8,6, color="red", linetype = "dashed") +
                     #scale_shape_manual(values=c(19, 0, 1, 2)) +
                     #scale_color_manual(values = getPalette(colourCount)) +
scale_color_manual(values = c("firebrick", "orange", "gold", "blue", "cyan", "magenta")) +
                     #scale_color_brewer(values = getPalette(colourCount)) +
                     #geom_line(mapping = aes(t_SM ~ Ratio, color = Model_SM)) +
                     #scale_linetype_manual(values = 2) +
                     coord_cartesian(xlim=c(-1, 44), expand = F) + #ylim=c(-0, 12.5)
                     #scale_y_discrete(limits=c(0,-3,-6,-9,-12,-15,-18)) +
                     #scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                                   #labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
                     #annotation_logticks(sides = "b", size = 0.1) + # adds linier tickmarks
                     #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
                     background_grid(major = "none", minor = "none") + # add thin horizontal lines
                     panel_border(colour = "black", size = 1) + # and a border around each panel
                     labs (x=("Time (days)"), y=("Counts"), title = i) +
                     #theme(plot.title = element_text(vjust = - 13, hjust = 0.04)) +
                     theme(panel.spacing = unit(0.1, "lines")) +
                     theme(legend.title=element_text(size=12)) +
                     theme(legend.text=element_text(size=6)) +
                     theme(text=(element_text(size=12))) +
                     theme(axis.text=(element_text(size=12))) +
                     theme(panel.background = element_blank()) +
                     theme(legend.text = element_text(face = "italic")) +
                     theme(aspect.ratio=1))
  
suppressWarnings(print(FigTotal))
  #Sys.sleep(1)
  
  suppressWarnings(dev.copy(pdf, output))
  suppressWarnings(dev.off())
  
  
  #This needs to be here for some annoying reason, loop terminates otherwise
  #quartz()
  
  suppressWarnings(FigStrains <- ggplot(StrainAmplicon, aes(Sample, Total)) + #, message=FALSE
                     geom_point(mapping = aes(), size = 2, shape = 1) +
                     coord_cartesian(xlim = c(), expand = F) + #ylim=c(-0, 12.5)
                     #scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                                   #labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
                     background_grid(major = "none", minor = "none") + # add thin horizontal lines
                     panel_border(colour = "black", size = 1) + # and a border around each panel
                     labs (x=("Strain Genotype sample"), y=("Counts"), title = i) +
                     theme(panel.spacing = unit(0.1, "lines")) +
                     theme(legend.title=element_text(size=12)) +
                     theme(legend.text=element_text(size=6)) +
                     theme(text=(element_text(size=10))) +
                     theme(axis.text=(element_text(size=6))) +
                     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                     theme(panel.background = element_blank()) +
                     theme(legend.text = element_text(face = "italic")) +
                     theme(aspect.ratio=0.25) +
                     theme(legend.position ="top"))
  
  suppressWarnings(print(FigStrains))
  #Sys.sleep(1)
  
  suppressWarnings(dev.copy(pdf, output2))
  suppressWarnings(dev.off())
  
}

#Lets make a heatmap
V <- as.vector(High_RA_means500$ID)
High_RA_means500M <- subset(High_RA_means500, select = -c(Unknown.barcode, Sum, ID, N))
#Reformat to Matrix
MatrixAll2 <- as.matrix(High_RA_means500M)
class(MatrixAll2) <-"numeric"
#MatrixAll2 <- MatrixAll2[,-1]
row.names(MatrixAll2)=V
MatrixAll2 <- t(MatrixAll2)
#Lets change INDEX of samples
INDEX = as.data.frame(rownames(MatrixAll2))
colnames(INDEX) <- c("ID")
V <- join(INDEX, Indexing)
V <- as.vector(V$Sample)
row.names(MatrixAll2)=V

#Logtransform matrix
MatrixAll2_log <- log10(MatrixAll2+1)

#Make a heatmap without clustering
quartz()
colRamp <- colorRampPalette(c("white", "red", "black", "green"), space="rgb")(64)
heatmap(MatrixAll2, col = colRamp, Rowv = NA, Colv = NA, cexRow=0.2, cexCol=0.2, na.rm = TRUE, distfun = dist)
dev.copy(pdf, "Plots/UnknownHeat_raw.pdf")
dev.off()

MatrixAll2_cF <- scale(MatrixAll2, scale = T, center = F) #center = F, scale = T)
MatrixAll2_cT <- scale(MatrixAll2, scale = T, center = T) #center = T, scale = T)

#Change NAs to 0 for clustering
MatrixAll2_cT[is.na(MatrixAll2_cT)] <- 0
MatrixAll2_cF[is.na(MatrixAll2_cF)] <- 0
MatrixAll2[is.na(MatrixAll2)] <- 0
MatrixAll2_log[is.na(MatrixAll2_log)] <- 0

#2#Apply clustering (Eucladian distances)
quartz()
heatmap.2(MatrixAll2_cF, col = colRamp, scale="none", revC = T, margins = c(8,6),
          xlab = "UnknowAmplicon_ID", ylab = "Sample", trace="none", cexRow=0.2, cexCol=0.1,
          symm=F,symkey=T,symbreaks=F)
dev.copy(pdf, "Plots/UnknownHeat_cF.pdf")
dev.off()

heatmap.2(MatrixAll2, col = colRamp, scale="none", revC = T, margins = c(8,6),
          xlab = "UnknowAmplicon", ylab = "Sample", trace="none", cexRow=0.2, cexCol=0.1,
          symm=F,symkey=T,symbreaks=F)
dev.copy(pdf, "Plots/UnknownHeat_Raw2.pdf")
dev.off()

#By looking them over i think most are artefacts associate with one strain, here are the once that are not:
Real <- c("6", "9", "10", "11", "62", "63") 
Real <- subset.data.frame(High_RA_means500, grepl(pattern = '\\b6\\b|\\b9\\b|\\b10\\b|\\b11\\b|\\b62\\b|\\b63\\b', High_RA_means500$ID))
Real

#SEARCH for  a sequence####
DNA <- "CCTCAAACCCCATCGAATACAAACTCTCCATATGCATCAAAAGCTGTGCCGCCTTATCACTCGCCCCCTTTTCCAATGAATTGGCATAGGCATTAATAATCGCACCAAATGTCAATGCGTCGGGTAGTAAATCATCATCCCCCTCATTGTACAGACGGCCCATTTCATCTAGTAAAAACTTGGCACGATCGGCGCAATTCTTGTCACGGGAACGTGCTACTGCGTTGATTACCGCGTTGAAGGATCGTGCATCGGGTTTGCAGTTTATGTCTAGTCCACTATCGTTATTGTAAAAGTCGTACATGTTACGGAGCAATTCTTCTGCTTTACGTGCCGATCCTCGTTCATTGCTACGTGCCCATGCCGTGATGACGGAGTTGTAACTAATTTTGTTTGGTTTGACATGATCCCATTTTCCATCGATATTATGTTTGGGAGAGGAGTAGAGGTGCAGCATTTGATCGAGGAGGGCTTGTGCGTGCTTTGCGGAATCTTTGCGTGTGCATCGTGCC"
DNApresenet <- subset.data.frame(High_RA_means100, grepl(print(DNA), High_RA_means100$Unknown.barcode))
DNApresenet

write.table(High_RA_means500, "Results/unknow_allels_AbuntantAllS.txt", sep='\t',  col.names=TRUE, row.names =FALSE)

#Adapth the plots to the abundant Mastermix samples instead####
MM_all <- subset.data.frame(High_RA_means, grepl(pattern = '\\b[8-9][0-3]\\b|\\b[0-8][0-9]\\b|\\b[0-9]\\b', High_RA_means$ID))
#to select 124 firs use pattern = '\\b12[0-4]\\b|\\b1[0-1][0-9]\\b|\\b[0-9][0-9]\\b|\\b[0-9]\\b'

Indexed_MM_all <- data.frame(ID=character(),
                             Index=character(),
                             Total=numeric(),
                             stringsAsFactors=FALSE)
#Lets change structure so we can index samples and plot them.
for (i in InputFiles) {
  Name <- gsub(".tsv", "", i) 
  N <- nrow(MM_all)
  ID <- rep(Name, each = N)
  Index <- MM_all$"ID"
  Total <- MM_all[grep(Name, names(MM_all))]
  colnames(Total) <- c("Total")
  Column <- cbind.data.frame(ID, Index, Total)
  Indexed_MM_all <- rbind(Indexed_MM_all, Column)
}

Indexing <- read.delim("../Indexing.txt", header=TRUE)
Indexed_MM_all2 <- join(Indexed_MM_all, Indexing, by = "ID")

#Lets plot them
df_uniq <- unique(Indexed_MM_all2$Index)
colourCount = length(unique(Indexed_MM_all2$Replicates))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
df_uniq
Indexed_MM_all2$Replicate <- as.character(as.integer(Indexed_MM_all2$Replicate))
#i <- "6"

for (i in df_uniq) {
  OneAmplicon <- subset.data.frame(Indexed_MM_all2, grepl(paste("\\b",i,"\\b", sep = ""), Indexed_MM_all2$Index))
  StrainAmplicon <- subset.data.frame(OneAmplicon, grepl('Strain', OneAmplicon$Experiment))
  ExperimentAmplicon <- subset.data.frame(OneAmplicon, grepl(pattern = 'GP|VG', OneAmplicon$Experiment))
  output <- paste("Plots/HighAbundance/",i,".pdf", sep="")
  output2 <- paste("Plots/HighAbundance/",i,"_Strains.pdf", sep="")
  
  
  #This needs to be here for some annoying reason, loop terminates otherwise
  #quartz()
  
  suppressWarnings(FigTotal <- ggplot(ExperimentAmplicon, aes(Timepoint, Total), message=FALSE) +
                     geom_point(mapping = aes(color = Replicate), size = 2, shape = 1) + #shape = Barcode,
                     facet_wrap(~Experiment + Treatment) +
                     theme(strip.text.x = element_text(size=12, angle=0, face = "italic"),
                           strip.background = element_rect(colour="white", fill="white")) +
                     geom_path(mapping = aes(color = Replicate), size = 0.5) +
                     #scale_color_brewer(palette = "Paired") +
                     #geom_errorbar(aes(xmin = (Low), xmax = (High), width=.4), position = position_dodge(0)) +
                     #geom_vline(xintercept = 0.86, color="red") +
                     #geom_vline(xintercept = 8,6, color="red", linetype = "dashed") +
                     #scale_shape_manual(values=c(19, 0, 1, 2)) +
                     #scale_color_manual(values = getPalette(colourCount)) +
                     scale_color_manual(values = c("firebrick", "orange", "gold", "blue", "cyan", "magenta")) +
                     #scale_color_brewer(values = getPalette(colourCount)) +
                     #geom_line(mapping = aes(t_SM ~ Ratio, color = Model_SM)) +
                     #scale_linetype_manual(values = 2) +
                     coord_cartesian(xlim=c(-1, 44), expand = F) + #ylim=c(-0, 12.5)
                     #scale_y_discrete(limits=c(0,-3,-6,-9,-12,-15,-18)) +
                     #scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                     #labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
                     #annotation_logticks(sides = "b", size = 0.1) + # adds linier tickmarks
                     #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
                     background_grid(major = "none", minor = "none") + # add thin horizontal lines
                     panel_border(colour = "black", size = 1) + # and a border around each panel
                     labs (x=("Time (days)"), y=("Counts"), title = i) +
                     #theme(plot.title = element_text(vjust = - 13, hjust = 0.04)) +
                     theme(panel.spacing = unit(0.1, "lines")) +
                     theme(legend.title=element_text(size=12)) +
                     theme(legend.text=element_text(size=6)) +
                     theme(text=(element_text(size=12))) +
                     theme(axis.text=(element_text(size=12))) +
                     theme(panel.background = element_blank()) +
                     theme(legend.text = element_text(face = "italic")) +
                     theme(aspect.ratio=1))
  
  suppressWarnings(print(FigTotal))
  #Sys.sleep(1)
  
  suppressWarnings(dev.copy(pdf, output))
  suppressWarnings(dev.off())
  
  
  #This needs to be here for some annoying reason, loop terminates otherwise
  #quartz()
  
  suppressWarnings(FigStrains <- ggplot(StrainAmplicon, aes(Sample, Total)) + #, message=FALSE
                     geom_point(mapping = aes(), size = 2, shape = 1) +
                     coord_cartesian(xlim = c(), expand = F) + #ylim=c(-0, 12.5)
                     #scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                     #labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
                     background_grid(major = "none", minor = "none") + # add thin horizontal lines
                     panel_border(colour = "black", size = 1) + # and a border around each panel
                     labs (x=("Strain Genotype sample"), y=("Counts"), title = i) +
                     theme(panel.spacing = unit(0.1, "lines")) +
                     theme(legend.title=element_text(size=12)) +
                     theme(legend.text=element_text(size=6)) +
                     theme(text=(element_text(size=10))) +
                     theme(axis.text=(element_text(size=6))) +
                     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                     theme(panel.background = element_blank()) +
                     theme(legend.text = element_text(face = "italic")) +
                     theme(aspect.ratio=0.25) +
                     theme(legend.position ="top"))
  
  suppressWarnings(print(FigStrains))
  #Sys.sleep(1)
  
  suppressWarnings(dev.copy(pdf, output2))
  suppressWarnings(dev.off())
  
}

#Lets make a heatmapp
V <- as.vector(MM_all$ID)
MM_allM <- subset(MM_all, select = -c(Unknown.barcode, Sum, ID))
#Reformat to Matrix
MatrixAll2 <- as.matrix(MM_allM)
class(MatrixAll2) <-"numeric"
#MatrixAll2 <- MatrixAll2[,-1]
row.names(MatrixAll2)=V
MatrixAll2 <- t(MatrixAll2)
#Lets change INDEX of samples
INDEX = as.data.frame(rownames(MatrixAll2))
colnames(INDEX) <- c("ID")
V <- join(INDEX, Indexing)
V <- as.vector(V$Sample)
row.names(MatrixAll2)=V

MatrixAll2_log <- log10(MatrixAll2+1)

#Make a heatmap without clustering
quartz()
colRamp <- colorRampPalette(c("white", "red", "black", "green"), space="rgb")(64)
heatmap(MatrixAll2, col = colRamp, Rowv = NA, Colv = NA, cexRow=0.2, cexCol=0.2, na.rm = TRUE, distfun = dist)
dev.copy(pdf, "Plots/MM_Heat_raw.pdf")
dev.off()

MatrixAll2_cF <- scale(MatrixAll2, scale = T, center = F) #center = F, scale = T)
MatrixAll2_cT <- scale(MatrixAll2, scale = T, center = T) #center = T, scale = T)

#Change NAs to 0 for clustering
MatrixAll2_cT[is.na(MatrixAll2_cT)] <- 0
MatrixAll2_cF[is.na(MatrixAll2_cF)] <- 0
MatrixAll2[is.na(MatrixAll2)] <- 0
MatrixAll2_log[is.na(MatrixAll2_log)] <- 0

#2#Apply clustering (Eucladian distances)
quartz()
heatmap.2(MatrixAll2_cF, col = colRamp, scale="none", revC = T, margins = c(8,6),
          xlab = "UnknowAmplicon_ID", ylab = "Sample", trace="none", cexRow=0.2, cexCol=0.2,
          symm=F,symkey=T,symbreaks=F)
dev.copy(pdf, "Plots/MM_Heat_cF.pdf")
dev.off()

heatmap.2(MatrixAll2, col = colRamp, scale="none", revC = T, margins = c(8,6),
          xlab = "UnknowAmplicon", ylab = "Sample", trace="none", cexRow=0.2, cexCol=0.2,
          symm=F,symkey=T,symbreaks=F)
dev.copy(pdf, "Plots/MM_Heat_Raw2.pdf")
dev.off()

write.table(MM_all, "Results/unknow_allels_MM.txt", sep='\t',  col.names=TRUE, row.names =FALSE)

#amplicon artefacts of VG90
#109, 113, 160, 155, 171, 149, 157, 110, 156, 166, 91, 105, 154, 98, 153, 78, 82, 92, 76, 130, 67



