#Houskeeping the loading packages####
#Lets clear old objects
rm(list=ls())
getwd()
#and set the work directory in case we have moved around or opened another project
setwd("~/Users/xanbjg/Documents/R/Cu_evolution/Genomics")
dir()

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
library("Hmisc")
#install.packages("naniar")
#Read in data and look at it####

#The strain-per-strain fitness estimates and the barcoded once
Barcode_Fitness <- read.table(file = "Input/Fitness_DL.txt", sep = '\t', header = TRUE)
sapply(Barcode_Fitness, class)

#Correlation of GP-Cu
Fitness_Cu <- subset.data.frame(Barcode_Fitness, grepl('Copper', Barcode_Fitness$Treatment))
Fitness_Cu_GP <- subset.data.frame(Fitness_Cu, grepl('GP', Fitness_Cu$Population))
ml <- lm(Fitness_Cu_GP$Growth_Rate~Fitness_Cu_GP$Barcode_mean)
summary(ml)
summary(ml)$r.squared

#Correlation of VG-Cu
Fitness_Cu_VG <- subset.data.frame(Fitness_Cu, grepl('VG', Fitness_Cu$Population))
ml <- lm(Fitness_Cu_VG$Growth_Rate~Fitness_Cu_VG$Barcode_mean)
summary(ml)$r.squared
summary(ml)

#Correlation of GP-C
Fitness_C <- subset.data.frame(Barcode_Fitness, grepl('Control', Barcode_Fitness$Treatment))
#Removing strain GP2-4_42 which is an outlier with very low growth rate
Fitness_C_GP <- subset.data.frame(Fitness_C, grepl('GP', Fitness_C$Population))
Fitness_C_GP <- Fitness_C_GP[-grep("GP2-4_42",Fitness_C_GP$Strain),]
ml <- lm(Fitness_C_GP$Growth_Rate~Fitness_C_GP$Barcode_mean)
summary(ml)$r.squared
summary(ml)
#Correlation of VG-C
Fitness_C_VG <- subset.data.frame(Fitness_C, grepl('VG', Fitness_C$Population))
#Removing strain VG1-2_63 which is an outlier with very low growth rate
Fitness_C_VG <- Fitness_C_VG[-grep("VG1-2_63",Fitness_C_VG$Strain),]
ml <- lm(Fitness_C_VG$Growth_Rate~Fitness_C_VG$Barcode_mean)
summary(ml)$r.squared
summary(ml)
#Then we can move on
Fitness_C <- subset.data.frame(Fitness_C, select = c("Strain", "Barcode_mean"))
Fitness_Cu <- subset.data.frame(Fitness_Cu, select = c("Strain", "Barcode_mean"))
colnames(Fitness_Cu) <- c("Strain", "Barcode_Cu")
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

#Make a matrix for correlative stats and plots
W <- as.vector(Meta$Strain)
Meta2 <- as.matrix(Meta)
class(Meta2) <-"numeric"
Meta2 <- Meta2[,-1]
row.names(Meta2)=W
Meta2 <- t(Meta2)


#Pvalues adjusted based on Holmes correction ()
Pvalues <- rcorr(as.matrix(t(Meta2)), type = "pearson")
P <- print(Pvalues$P)
r <- print(Pvalues$r)
#n = NROW(P)*NCOL(P)

Variable <- colnames(P)
Variable.pair <- paste(Variable[row(P)], Variable[col(P)], sep="_vs_")
i <- lower.tri(P)
P_holms <- data.frame(Variable.pair[i], p.value=P[i], r.value=r[i])
P_holms$Holm <- p.adjust(P, method="holm")
P_holms
#Okey lets export list of significant correlations
write.table(P_holms, "Results/P_holms_pheno.txt", sep='\t',  col.names=TRUE, row.names =FALSE)

#PCAs####
#Add Population/strain vector
IndexShort <- subset.data.frame(Fitness, select = c("Population", "Strain"))
IndexShort2 <- IndexShort[!(IndexShort$Strain %in% c('GP2-4_45', 'GP2-4_46',"VG1-2_65", "VG1-2_99")), ]
#IndexShort2 <- sort(IndexShort2$Strain)
IndexShort2 <- IndexShort2 %>% arrange(Strain)
#unique(Meta$Strain)
#Vectors for coloring PCA
P <- IndexShort2$Population
ID <- IndexShort2$Strain
ID
colnames(Meta2)
#Run PCA on many factors
PCA <- prcomp(na.omit(t(Meta2)), scale = TRUE)
summary(PCA)

dev.off()
dev.new(width=5, height=5)
PCA_all <-ggbiplot(PCA, ellipse=TRUE, labels=ID, groups=P, alpha=0, varname.size = 2, labels.size=2) + #labels=ID, groups=P
  #coord_cartesian(ylim=c(-0.01,0.31), expand = F) +
  scale_color_manual(values=c("gray14", "#B85633")) +
  #geom_point(aes(colour=P), size = 5) +
  background_grid(major = "none", minor = "none", col) + # add thin horizontal lines 
  theme(panel.background = element_rect(fill = 'white', color = 'white')) +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x="PC1 (30.5%)", y="PC2 (23.3%)") +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=16)) +
  theme(text=(element_text(size=16))) +
  theme(axis.text=(element_text(size=16))) +
  #theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  #theme(legend.position = "none")
  xlim(-3.9, 2.2) +
  ylim(-2.7, 2.7) +
  theme(legend.position = c(0.995, 0.995),
        legend.justification = c("right", "top"),
        #legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"),
        legend.key.size = unit(0.3, "cm"))

print(PCA_all)
dev.copy(pdf, "Plots/PCApheno.pdf")
dev.off()

#Scale data
Meta2_cF <- scale(Meta2, scale = F, center = F) #center = F, scale = T)
Meta2_cT <- scale(Meta2, scale = T, center = T) #center = T, scale = T)
Meta2_all <- Meta2
#Change NA to 0
Meta2_cT[is.na(Meta2_cF)] <- 0
Meta2_cF[is.na(Meta2_cT)] <- 0
Meta2_all[is.na(Meta2_all)] <- 0

#Heatmaps####
#Make a heatmap of raw data without clustering
quartz()
colRamp <- colorRampPalette(c("white", "red", "black", "green"), space="rgb")(64)
heatmap(Meta2_all, col = colRamp, Rowv = NA, Colv = NA, cexRow=0.5, cexCol=0.5, na.rm = TRUE, distfun = dist)
dev.copy(pdf, "Plots/Heat_raw.pdf")
dev.off()

#2#Apply clustering (Euclidean distances)
quartz()
heatmap.2(Meta2_all, col = colRamp, scale="row", revC = T, margins = c(4,10),
          xlab = "Strain", ylab = "Phenotype", trace="none", cexRow=0.5, cexCol=0.4,
          symm=F,symkey=T,symbreaks=F)
dev.copy(pdf, "Plots/MM_Heat_absolut.pdf")
dev.off()


quartz()
heatmap.2(Meta2_cF, col = colRamp, scale="row", revC = T, margins = c(4,10),
          xlab = "Strain", ylab = "Phenotype", trace="none", cexRow=0.5, cexCol=0.4,
          symm=F,symkey=T,symbreaks=F)
dev.copy(pdf, "Plots/MM_Heat_cF_Scaled.pdf")
dev.off()

#CorrelationMatrix#####
#This is the R values between all parameters
quartz()
CorrelationMAtrix<- cor(t(Meta2), use = "complete.obs")
heatmap.2(CorrelationMAtrix, col = colRamp, scale="none", revC = T, margins = c(10,10),
          xlab = "Phenotype", ylab = "Phenotype", trace="none", cexRow=1, cexCol=1,
          symm=F,symkey=T,symbreaks=F)
dev.copy(pdf, "Plots/Heat_Corr_pheno.pdf")
dev.off()
