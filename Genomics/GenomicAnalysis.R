#Houskeeping the loading packages####
#Lets clear old objects
rm(list=ls())
getwd()
#and set the work directory in case we have moved around or opened another project
#setwd("~/Users/xanbjg/Documents/R/Cu_evolution/Genomics")
dir()

# loading packages containing functions. Access the functions that these contain.
# check if you have these under packages
# you probably don't have it so you have to install it before. Everytime I start R studio
# I have to run these librarys. 

#packages & functions
library(drc)
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
library(tidyverse)
library(lubridate) # useful for working with dates
library(cowplot) # useul for combining multiple plots
library(ggthemes)
library(broom)
library(ggpubr)
library(gridExtra)
#library(ggpmisc)
library(gplots)
#library(staplr)
library(staplr)
library(dplyr)
library(naniar)
#install.packages("naniar")
#Read in data and look at it####

#The strain-per-strain fitness estimates and the barcoded once
Barcode_Fitness <- read.table(file = "Input/Fitness.txt", sep = '\t', header = TRUE)
sapply(Barcode_Fitness, class)

Fitness_Cu <- subset.data.frame(Barcode_Fitness, grepl('Copper', Barcode_Fitness$Treatment))
Fitness_Cu <- subset.data.frame(Fitness_Cu, select = c("Strain", "Barcode_mean"))
colnames(Fitness_Cu) <- c("Strain", "Barcode_Cu")
Fitness_C <- subset.data.frame(Barcode_Fitness, grepl('Control', Barcode_Fitness$Treatment))
Fitness_C <- subset.data.frame(Fitness_C, select = c("Strain", "Barcode_mean"))
colnames(Fitness_C) <- c("Strain", "Barcode_C")
Barcode_Fitness2 <- join(Fitness_C, Fitness_Cu, by = "Strain")

#Lets also comute the relative inhibition of growth via barcode observations
Barcode_Fitness2$Barcoded_CuInhibition <- Barcode_Fitness2$Barcode_Cu/Barcode_Fitness2$Barcode_C
plot(Barcode_Fitness2$Barcode_Cu~Barcode_Fitness2$Barcoded_CuInhibition)
#super correlated so this will not affect analysis much

#Metadata for linking genomic ID and adding other phenotypes
Fitness <- read.table(file = "Input/DRCpredictions.txt", sep = '\t', header = TRUE)
sapply(Fitness, class)

#Lets also compute the responsrange between EC05 and EC95
Fitness$ResponseRange_Cu <- Fitness$EC95/Fitness$EC05
plot(Fitness$ResponseRange_Cu~Fitness$EC50)
#So the slope/responsrange and EC50 is not correlated so they may capture diffrent modes of tolerance

#Lets reduce to data to Cu phenotype relevant traits that are not entirely autocorrelated
Fitness2 <- subset.data.frame(Fitness, select = c("Strain", "EC05", "EC50", 
                                                         "EC95", "ResponseRange_Cu", "Growth_Rate", "FvFm",
                                                         "Surface_Vol_ratio", "CuGrowth"))
Meta <- join(Barcode_Fitness2, Fitness2, by = "Strain")

Index <- read.table(file = "Input/StrainMeta.txt", sep = '\t', header = TRUE)

#Read in information about the copper metabolism genes
Annotations <- read.table(file = "Input/CuMetabolism_SM_short.txt", sep = '\t', header = TRUE)
sapply(Annotations, class)

#Make a list of input files of gene coverage across Cu metabolism genes
InputFiles <- list.files(path = "~/Documents/R/Cu_evolution/Genomics/Input/CuGeneCoverage")
InputFiles

#myData <- read.table(file = "Input/CuGeneCoverage/CuGenesOfInterest.lst", sep = '\t', header = TRUE)
#P14203_101_CuCoverage.tsv

#sapply(Indexed, class)
#library(RColorBrewer)
#n <- 97
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))

WGScoverage  <- NULL;
#i <-"P14203_105_CuCoverage.tsv"
for (i in InputFiles) {
  #Import data (too be looped)
  #i <-"P14203_105_CuCoverage.tsv"
  input=paste("Input/CuGeneCoverage/",i,"", sep="")
  myData <- read.table(file = input, sep = '\t', header = TRUE)
  
  #Generate sample ID for indexing experimental samples and add metadata
  Name <- i
  ID <- gsub('_CuCoverage.tsv', '', Name)
  N <- nrow(myData)
  Sample <- rep(ID, each = N)
  Sample <- as.data.frame(Sample)
  colnames(Sample) <- c("ID")
  myData2 <- join(Sample, Index, by = "ID")

  #and fuse with data
  Indexed <- cbind(myData2, myData)
  
  #Lets normalize in two ways
  #1 versus the genome
  MeanGenomeCoverage <- subset.data.frame(Indexed, grepl('Genome', Indexed$Gene))
  MeanGenomeCoverage <- MeanGenomeCoverage$Model.mean
  
  Indexed$Coverage_vs_genome <- Indexed$Model.mean/MeanGenomeCoverage
  
  #2 contig based
  Indexed$Contig.mean <- as.numeric(as.character(Indexed$Contig.mean))
  Indexed$Coverage_vs_contig <- Indexed$Model.mean/Indexed$Contig.mean
  
  #Add the data to growing dataframe
  WGScoverage <- rbind(WGScoverage, Indexed)
  
  #Okey lets also make some plots
  Titel <- unique(Indexed$Strain)
  output3 <- paste("Plots/",ID,".pdf", sep="")
  output3
  
  Fig1 <- ggplot(Indexed, stats = "identity", aes(Gene, fill = Gene)) +
    geom_col(data=Indexed, aes(Gene, Coverage_vs_genome)) +
    geom_hline(yintercept = 1) +
    background_grid(major = "none", minor = "none") + # add thin horizontal lines
    #theme_classic() +
    panel_border(colour = "black", size = 1) + # and a border around each panel
    labs (title = Titel) +
    #coord_cartesian(ylim=c(0.0001, 0.01), expand = F) + #ylim=c(-0, 12.5)
    #scale_y_discrete(limits=c(0.0001,0.001,0.01,0.1,1)) +
    #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #labels = trans_format("log10", math_format(10^.x))) +
    #coord_flip() +
    theme(plot.title = element_text(vjust = -6, hjust = 0.1)) +
    theme(panel.spacing = unit(0.1, "lines")) +
    theme(legend.title=element_blank()) +
    theme(legend.text=element_text(size=10)) +
    theme(text=(element_text(size=10))) +
    theme(axis.text=(element_text(size=3))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.background = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(legend.text = element_text(face = "italic")) +
    theme(aspect.ratio=0.25) +
    #stat_cor(aes(color = scientific_name), label.x = 3)
    theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
    theme(legend.position ="none")
          #legend.justification = c("right", "top"),
          #legend.box.just = "right",
          #scale_fill_manual(breaks=c(3)),
          #legend.margin = margin(2, 2, 2, 2),
          #guides(color = FALSE), #,color = FALSE, size = FALSE,
          #legend.box.background = element_rect(fill='white'),
          #legend.background = element_blank(),
          #legend.key = element_rect(fill = NA, color = NA),
          #legend.spacing.x=unit(0, "cm"),
          #legend.spacing.y=unit(0, "cm"))
  
 # print(Fig1)
  
 # dev.copy(pdf, output3)
 # dev.off()
  
}

#Heres the avarage genome coverage per strain
WGScoverage2 <- subset.data.frame(WGScoverage, grepl('Genome', WGScoverage$Gene))
mean(WGScoverage2$Model.mean)
#Lets see what the per strain avarage coverage is, we probably want to normalize to that
#Remove the genome data
WGScoverage3 <- WGScoverage[ !(WGScoverage$Gene %in% c('Genome')), ]

#Lets also change 0 to NA  and remove samples with No coverage from analysis 
#(lots of this is from GP_44 with very low genomic coverage 0.6)
WGScoverage3$Coverage_vs_genome[WGScoverage3$Coverage_vs_genome==0] <- NA
WGScoverage3 <- WGScoverage3[complete.cases(WGScoverage3$Coverage_vs_genome), ]


#OPTIONAL change to pool copies of genes in close proximity on contig
#Basicly go through and change WGScoverage4 to WGScoverage3 and "gene" to "Putative_identity_INDEX"

WGScoverage4 <- join(WGScoverage3, Annotations, by = "Gene", type="right")
WGScoverage4 <- ddply(WGScoverage4, c("Putative_identity_INDEX", "Strain"), summarise,
                         Sum = sum(Coverage_vs_genome), N = n(), sd = sd(Coverage_vs_genome), 
                         Min = min(Coverage_vs_genome), 
                         Max = max(Coverage_vs_genome))

#Here it is possible to change
WGScoverageMean <- ddply(WGScoverage4, c("Putative_identity_INDEX"), summarise,
                 Mean = mean(Sum), N = n(), sd = sd(Sum), 
                 Min = min(Sum), 
                 Max = max(Sum))

WGScoverageMean$Range <- WGScoverageMean$Max-WGScoverageMean$Min

#Here it is possible to change 

#Lets round the numbers of for a tabel
WGScoverageMean_export <- WGScoverageMean
WGScoverageMean_export$Mean <- signif(WGScoverageMean_export$Mean, digits = 2)
WGScoverageMean_export$sd <- signif(WGScoverageMean_export$sd, digits = 2)
WGScoverageMean_export$Min <- signif(WGScoverageMean_export$Min, digits = 2)
WGScoverageMean_export$Max <- signif(WGScoverageMean_export$Max, digits = 2)
WGScoverageMean_export$Range <- signif(WGScoverageMean_export$Range, digits = 2)
write.table(WGScoverageMean_export, file = "Results/WGScoverageMean.txt", sep = '\t', col.names = TRUE)
#Add back the annotations (only for all genes)
#WGScoverageMean2 <- join(WGScoverageMean, Annotations, by = "Putative_identity_INDEX", type="left")

Fig1 <- ggplot(WGScoverageMean, stats = "identity", aes(Putative_identity_INDEX, fill = Putative_identity_INDEX)) +
  geom_col(data=WGScoverageMean, aes(Putative_identity_INDEX, Mean)) +
  geom_hline(yintercept = 1) +
  geom_errorbar(aes(ymin = Min, ymax = Max, width=1), color="black") +
  scale_y_continuous(trans = "log2") + #change the scale on y axis
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = "Avarage_min_max",) +
  theme(plot.title = element_text(vjust = -6, hjust = 0.1)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=3))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  #stat_cor(aes(color = scientific_name), label.x = 3)
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(legend.position ="none")

print(Fig1)

dev.copy(pdf, "Plots/AvarageCoverage_Summed.pdf")
dev.off()

#USe WGScoverage4-Putative_identity_INDEX-Sum for summed orthologs, or WGScoverage3-gene-coverage_vs_genome for genemodels
Fig1B <- ggplot(WGScoverage4, stats = "identity", aes(Putative_identity_INDEX, fill = Putative_identity_INDEX)) +
  geom_boxplot(data=WGScoverage4, aes(Putative_identity_INDEX, Sum), lwd=0.2, outlier.size=0.2) +
  geom_hline(yintercept = 1) +
  #geom_errorbar(aes(ymin = Min, ymax = Max, width=1), color="black") +
  scale_y_continuous(trans = "log2") + #change the scale on y axis
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = "AllStrains") +
  theme(plot.title = element_text(vjust = -6, hjust = 0.1)) +
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
  theme(legend.position ="none")

print(Fig1B)

dev.copy(pdf, "Plots/BoxplotCoverage_summed.pdf")
dev.off()

#Changes in coverage filter####

#Okey we should sort away any genes that have less than 1 in min-max range
#since this suggests no change in copy number
WGScoverageMean_changing <- WGScoverageMean
WGScoverageMean_changing$Range[WGScoverageMean_changing$Range<1] <- NA
WGScoverageMean_changing <- WGScoverageMean_changing[complete.cases(WGScoverageMean_changing$Range), ]
#Drop all data save the mean values 
WGScoverageMean_changing <- WGScoverageMean_changing[,c("Putative_identity_INDEX", "Mean", "Range")]
#Lets sort out these genes using join
Changing_genes <- join(WGScoverage4, WGScoverageMean_changing, by = "Putative_identity_INDEX", type="right")
names(Changing_genes)[names(Changing_genes) == 'Sum'] <- "Coverage_vs_genome"

Fig1C <- ggplot(Changing_genes, stats = "identity", aes(Putative_identity_INDEX, fill = Putative_identity_INDEX)) +
  geom_boxplot(data=Changing_genes, aes(Putative_identity_INDEX, Coverage_vs_genome), lwd=0.2, outlier.size=0.2) +
  geom_hline(yintercept = 1) +
  #geom_errorbar(aes(ymin = Min, ymax = Max, width=1), color="black") +
  scale_y_continuous(trans = "log2") + #change the scale on y axis
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = "", x="Gene family summed by contig possition", y="Normalized Coverage") +
  theme(plot.title = element_text(vjust = -6, hjust = 0.1)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=7))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  #stat_cor(aes(color = scientific_name), label.x = 3)
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(legend.position ="none")

print(Fig1C)

dev.copy(pdf, "Plots/BoxplotChanging_Indexed.pdf")
dev.off()

#Genotype-Phenotype correlation#####

#First add relevant metadata
#Since this works the strain index is right here
sapply(Barcode_Fitness, class)
IndexShort <- subset.data.frame(Fitness, select = c("Population", "Strain"))
Changing_genes2 <- join(Changing_genes, Meta, by = "Strain", type = "left")
Changing_genes2 <- join(Changing_genes2, IndexShort, by = "Strain", type = "left")
sapply(Changing_genes2, class)

#There is some duplicated colums so i will extract data for graphs
#Changing_genes2 <- subset.data.frame(Changing_genes2, select = c("Putative_identity_INDEX", "Population", "Strain", "CuGrowth", "EC50", "Barcode_Cu", "Coverage_vs_genome"))
plot(Changing_genes2$Coverage_vs_genome~Changing_genes2$EC50)
unique(Changing_genes2$Putative_identity_INDEX)
unique(Changing_genes2$Strain)
sapply(Changing_genes2, class)
library(RColorBrewer)
n <- 58
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

FigReg <- ggplot(Changing_genes2, aes(EC50, Coverage_vs_genome)) +
  geom_point(mapping = aes(color = Strain, shape = Population), size = 1.2) + #color = Population, shape = Population
  facet_wrap(~Putative_identity_INDEX) +
  #labs(colour="Population-Timepoint") +
  stat_smooth(mapping = aes(), size = 0.5, method = 'lm', formula = y ~ x, se= TRUE) +
  geom_hline(yintercept=1) +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values=c(col_vector)) + #aesthetics = c("colour", "fill")) +
  scale_shape_manual(values=c(1, 2)) +
  #coord_cartesian(ylim=c(0.5, 28000), xlim=c(0.5, 28000), expand = F) + #ylim=c(-0, 12.5)
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                #labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim=c(6, 10), ylim=c(0.25,26), expand = F) + #ylim=c(-10000,+10000)
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
               labels = trans_format("log2", math_format(2^.x))) + #change the scale on y axis
  #annotation_logticks(sides = "b", size = 0.1) + # adds linier tickmarks
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x=expression("EC50 (" *mu ~ "M Cu)"), y=("Normalized Coverage")) + #title = i
  #theme(plot.title = element_text(vjust = - 13, hjust = 0.04)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  #theme(legend.title=element_text(size=3)) +
  theme(legend.text=element_text(size=6)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=6))) +
  theme(strip.text = element_text(size=3)) +
  theme(panel.background = element_blank()) +
  theme(strip.background =element_rect(fill='white')) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position="right")

print(FigReg)

dev.copy(pdf, "Plots/Reg_EC50.pdf")
dev.off()

FigReg2 <- ggplot(Changing_genes2, aes(Barcode_Cu, Coverage_vs_genome)) +
  geom_point(mapping = aes(color = Strain, shape = Population), size = 1.2) + #color = Population, shape = Population
  facet_wrap(~Putative_identity_INDEX) +
  #labs(colour="Population-Timepoint") +
  stat_smooth(mapping = aes(), size = 0.5, method = 'lm', formula = y ~ x, se= TRUE) + #Fits liniear regression to data
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values=c(col_vector)) + #aesthetics = c("colour", "fill")) +
  scale_shape_manual(values=c(1, 2)) +
  geom_hline(yintercept=1) +
  #coord_cartesian(ylim=c(0.5, 28000), xlim=c(0.5, 28000), expand = F) + #ylim=c(-0, 12.5)
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim=c(-0.2, 1), ylim=c(0.25,26), expand = F) + #ylim=c(-10000,+10000)
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) + #change the scale on y axis
  scale_x_discrete(limits=c(0,0.3,0.6,0.9)) +
  #annotation_logticks(sides = "b", size = 0.1) + # adds linier tickmarks
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x=expression("Barcoded copper growth rate"~(day^{-1})), y=("Normalized Coverage")) + #title = i
  #theme(plot.title = element_text(vjust = - 13, hjust = 0.04)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  #theme(legend.title=element_text(size=3)) +
  theme(legend.text=element_text(size=6)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=6))) +
  theme(strip.text = element_text(size=3)) +
  theme(panel.background = element_blank()) +
  theme(strip.background =element_rect(fill='white')) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position="right")

print(FigReg2)

dev.copy(pdf, "Plots/Reg_Barcoded_Cu.pdf")
dev.off()
#Okey lets try and build a matrix to correlate phenotype-genotype
#We also need a dataframe for cbinding Matrix data
MatrixAll <- data.frame(unique(Changing_genes$Putative_identity_INDEX))
colnames(MatrixAll) <- (c('Putative_identity_INDEX'))
ID <- unique(Changing_genes$Strain)
ID
#df_uniq2 <- subset(df_uniq, grepl("GP",df_uniq))
sapply(Changing_genes2, class)

#SOME problem with the Strain names in Changing_genes2, they wont match the old names... Creats issues downstream
for (i in ID) {
  #Lets also make a matrix format for further analysis
  #i <- "GP_2-4_26"
Coverage <- subset.data.frame(Changing_genes, grepl(paste("\\b",i,"\\b", sep = ""), Changing_genes$Strain), select = c("Putative_identity_INDEX", "Strain", "Coverage_vs_genome"))
Coverage <- subset(Coverage, select = -c(Strain))
colnames(Coverage) <- (c("Putative_identity_INDEX", i))

  #Then we just join this with growing matrix
  MatrixAll <- join(MatrixAll, Coverage, by = "Putative_identity_INDEX")
  
}

sapply(MatrixAll, class)

#Lets change the fucking header now
#MatrixHead <- colnames(MatrixAll)
#MatrixHead <- gsub("_P21502_[0-9]{3}_C12W1", "", MatrixHead)
#MatrixAll2 <- MatrixAll
#colnames(MatrixAll2) <- MatrixHead


#Reformat to Matrix
V <- as.vector(MatrixAll$Putative_identity_INDEX)
MatrixAll2 <- as.matrix(MatrixAll)
class(MatrixAll2) <-"numeric"
MatrixAll2 <- MatrixAll2[,-1]
row.names(MatrixAll2)=V
#MatrixAll2 <- t(MatrixAll2)

#And add other phenotypes
W <- as.vector(Meta$Strain)
Meta2 <- as.matrix(Meta)
class(Meta2) <-"numeric"
Meta2 <- Meta2[,-1]
row.names(Meta2)=W
Meta2 <- t(Meta2)

#We need to remove strains lacking data in either database
x <- colnames(Meta2)
y <- colnames(MatrixAll2)
#These samples are not in both datasets
x[!(x %in% y)]
y[!(y %in% x)]

#Delete them
Meta2 <- Meta2[,colnames(Meta2)!="GP2-4_44"]
Meta2 <- Meta2[,colnames(Meta2)!="VG1-2_65or99"]
Meta2 <- Meta2[,colnames(Meta2)!="VG1-2_99or65"]
Meta2 <- Meta2[,colnames(Meta2)!="GP2-4_45and46"]
MatrixAll2 <- MatrixAll2[,colnames(MatrixAll2)!="GP2-4_45"]
MatrixAll2 <- MatrixAll2[,colnames(MatrixAll2)!="GP2-4_46"]
colnames(Meta2)
colnames(MatrixAll2)

#Check that they are gone for real
x <- colnames(Meta2)
y <- colnames(MatrixAll2)
#These samples are not in both datasets
x[!(x %in% y)]
y[!(y %in% x)]

MatrixAll3 <- rbind(MatrixAll2, Meta2)

MatrixAll_cF <- scale(MatrixAll3, scale = F, center = F) #center = F, scale = T)
MatrixAll_cT <- scale(MatrixAll3, scale = T, center = T) #center = T, scale = T)

#Change NA to 0
MatrixAll_cT[is.na(MatrixAll_cT)] <- 0
MatrixAll_cF[is.na(MatrixAll_cF)] <- 0
MatrixAll3[is.na(MatrixAll3)] <- 0

#Heatmaps####
#Make a heatmap of raw data without clustering
quartz()
colRamp <- colorRampPalette(c("white", "red", "black", "green"), space="rgb")(64)
heatmap(MatrixAll3, col = colRamp, Rowv = NA, Colv = NA, cexRow=0.5, cexCol=0.5, na.rm = TRUE, distfun = dist)
dev.copy(pdf, "Plots/Heat_raw.pdf")
dev.off()

#2#Apply clustering (Euclidean distances)
quartz()
heatmap.2(MatrixAll3, col = colRamp, scale="row", revC = T, margins = c(4,10),
          xlab = "Strain", ylab = "Phenotype/genotype", trace="none", cexRow=0.5, cexCol=0.4,
          symm=F,symkey=T,symbreaks=F)
dev.copy(pdf, "Plots/MM_Heat_Scaled.pdf")
dev.off()

write.table(MatrixAll3, "Results/CorrMAtrix_genopheno.txt", sep='\t',  col.names=TRUE, row.names =FALSE)
#symm=F, symkey=F,symbreaks=F

#CorrelationMatrix#####
#This is the R values between all parameters
quartz()
CorrelationMAtrix<- cor(t(MatrixAll3), use = "complete.obs")
heatmap.2(CorrelationMAtrix, col = colRamp, scale="none", revC = T, margins = c(10,10),
          xlab = "Phenotype/genotype", ylab = "Phenotype/genotype", trace="none", cexRow=0.4, cexCol=0.4,
          symm=F,symkey=T,symbreaks=F)
dev.copy(pdf, "Plots/Heat_Corr.pdf")
dev.off()

#install.packages("Hmisc")
library("Hmisc")
#Pvalues adjusted based on Holmes correction
Pvalues <- rcorr(as.matrix(t(MatrixAll3)), type = "pearson")
P <- print(Pvalues$P)
#n = NROW(P)*NCOL(P)

Variable <- colnames(P)
Variable.pair <- paste(Variable[row(P)], Variable[col(P)], sep="_vs_")
i <- lower.tri(P)
P_holms <- data.frame(Variable.pair[i], p.value=P[i])
P_holms$Holm <- p.adjust(P, method="holm")
P_holms

P_holms$p.value[P_holms$p.value>0.05] <- NA
P_holms_sig <- P_holms[complete.cases(P_holms$p.value), ]
write.table(P_holms_sig, "Results/P_holms_sig.txt", sep='\t',  col.names=TRUE, row.names =FALSE)

#Okey lets export list of significant correlations

#PCAs####
#Add Population/strain vector
IndexShort2 <- IndexShort[ !(IndexShort$Strain %in% c('GP2-4_44', 'GP2-4_45', 'GP2-4_46', 'GP2-4_46',"VG1-2_65", "VG1-2_99")), ]

#Vectors for coloring PCA
P <- IndexShort2$Population
ID <- IndexShort2$Strain
#Run PCA on many factors
PCA <- prcomp(na.omit(t(MatrixAll3)), scale = TRUE)
summary(PCA)
dev.off()
dev.new(width=5, height=5)
PCA_all <-ggbiplot(PCA, ellipse=TRUE, labels=ID, groups=P, alpha=0, varname.size = 2, labels.size=2) +
  #coord_cartesian(ylim=c(-0.01,0.31), expand = F) +
  scale_color_manual(values=c("gray14", "#B85633")) +
  #geom_point(aes(colour=P), size = 5) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines 
  panel_border(colour = "black", size = 1) + # and a border around each panel
  #labs (x="PC1 (32.9%)", y="PC2 (18.2%)") +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  #theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=16))) +
  theme(axis.text=(element_text(size=16))) +
  #theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position = "none") +
  xlim(-2, 3) +
  ylim(-2.25, 2)

print(PCA_all)
dev.copy(pdf, "Plots/PCAall.pdf")
dev.off()


#The PCA is just a messier form of data then the heatmap of correlation coefficents. Also poor resolution in the two first axises.

#But as a final test how about we resolve the patterns in Cu genes to PCA1, and plot this against the PCA1 of the phenotypes?

Phenotypes <- MatrixAll3[-c(1:18),]
Genotypes <- MatrixAll3[c(1:18),]
PCA_Geno <- prcomp(na.omit(t(Genotypes)), scale = TRUE)
PCA_Pheno <-prcomp(na.omit(t(Phenotypes)), scale = TRUE)
summary(PCA_Geno)
summary(PCA_Pheno)

dev.off()
dev.new(width=5, height=5)
PCA1 <-ggbiplot(PCA_Geno, ellipse=TRUE, labels=ID, groups=P, alpha=0, varname.size = 2, labels.size=2) +
  #coord_cartesian(ylim=c(-0.01,0.31), expand = F) +
  scale_color_manual(values=c("gray14", "#B85633")) +
  #geom_point(aes(colour=P), size = 5) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines 
  panel_border(colour = "black", size = 1) + # and a border around each panel
  #labs (x="PC1 (32.9%)", y="PC2 (18.2%)") +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  #theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=16))) +
  theme(axis.text=(element_text(size=16))) +
  #theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position = "none") +
  xlim(-2, 3) +
  ylim(-2.25, 2)

print(PCA1)

PCA2 <-ggbiplot(PCA_Pheno, ellipse=TRUE, labels=ID, groups=P, alpha=0, varname.size = 2, labels.size=2) +
  #coord_cartesian(ylim=c(-0.01,0.31), expand = F) +
  scale_color_manual(values=c("gray14", "#B85633")) +
  #geom_point(aes(colour=P), size = 5) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines 
  panel_border(colour = "black", size = 1) + # and a border around each panel
  #labs (x="PC1 (32.9%)", y="PC2 (18.2%)") +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  #theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=16))) +
  theme(axis.text=(element_text(size=16))) +
  #theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position = "none") +
  xlim(-2, 3) +
  ylim(-2.25, 2)

print(PCA2)


#Okey lets tru and correlate first component of each dataset
PCA_Geno <- prcomp(na.omit((Genotypes)), scale = TRUE)
PCA_Pheno <-prcomp(na.omit((Phenotypes)), scale = TRUE)

#Get the principal compunents and change header to plot them
PCA_P <- as.data.frame(PCA_Pheno$rotation)
colnames(PCA_P) <- gsub('PC', 'P_PC', colnames(PCA_P))
PCA_G <- as.data.frame(PCA_Geno$rotation)
colnames(PCA_G) <- gsub('PC', 'G_PC', colnames(PCA_G))

#Merge data
PCA_GvsP<- cbind(PCA_P, PCA_G, IndexShort2)

plot(PCA_GvsP$G_PC1, PCA_GvsP$P_PC1)
text(PCA_GvsP$G_PC1, PCA_GvsP$P_PC1-1, labels=PCA_GvsP$Strain)

library(ggrepel)
ggplot(PCA_GvsP, aes(G_PC1, P_PC1)) +
  geom_point(aes(color = Population)) +
  geom_text_repel(aes(label = Strain))

ggplot(PCA_GvsP, aes(P_PC1, G_PC2)) +
  geom_point(aes(color = Population)) +
  geom_text_repel(aes(label = Strain))

ggplot(PCA_GvsP, aes(P_PC1, G_PC2)) +
  geom_point(aes(color = Population)) +
  geom_text_repel(aes(label = Strain))

ggplot(PCA_GvsP, aes(P_PC2, G_PC1)) +
  geom_point(aes(color = Population)) +
  geom_text_repel(aes(label = Strain))

ggplot(PCA_GvsP, aes(P_PC2, G_PC2)) +
  geom_point(aes(color = Population)) +
  geom_text_repel(aes(label = Strain))

#Cannot say that that was enlightening in any way. PCA are to crude

#Multiple regression####
#I am almost against even doing this given the dataspread and lack of correlation
FitData <- as.data.frame(t(MatrixAll3))
Variables <- colnames(FitData[c(1:18),])
FitData[, c(1:18)]
sapply(FitData, class)
fit <- lm(FitData$EC50 ~ FitData[, c(1)]+FitData[, c(2)]+FitData[, c(3)]+
          FitData[, c(4)]+FitData[, c(5)]+FitData[, c(6)]+
          FitData[, c(7)]+FitData[, c(8)]+FitData[, c(9)]+
          FitData[, c(10)]+FitData[, c(11)]+FitData[, c(12)]+
          FitData[, c(13)]+FitData[, c(14)]+FitData[, c(15)]+
          FitData[, c(16)]+FitData[, c(17)]+FitData[, c(18)],
          data=FitData)

summary(fit) # show results

library(MASS)
step <- stepAIC(fit, direction="both")
step$anova # display results

#So we can use for example the 5 most predictive genetic factors
Iterative_fit <- lm(FitData$EC50 ~ FitData[, c(13)]+FitData[, c(7)]+FitData[, c(1)]+
            +FitData[, c(16)]+FitData[, c(3)],
          data=FitData)

summary(Iterative_fit) # show results
print
colnames(FitData[, c(13)])
colnames(FitData)
#[13] "Metallothionein-like protein_C54"
# [7] "Cytochrome c6 isoform A_C29"
#[1] "Cu-transporting P1B-type ATPases_C15"
# 16 "Multimetal transporter (ZIP)_C111"  
# [3] "Cue5-4_homolog_C3"


#That was a pointless, and statistically unsound exploration. Even cheeroicking data does not make a multiple regression significant
#Clearly, across all strains, gene-copy number does not explain Cu tolerance, however the data is twisted and turned
#Still possible that some aspect of very tolerant individuals can be partially explained by copy-number variance
#Like the duplication of the ZIP multimetal transporter in VG1-2_81
###########