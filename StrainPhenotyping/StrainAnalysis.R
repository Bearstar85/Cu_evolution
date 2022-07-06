
#Lets clear the pbjects from day1
rm(list=ls())
getwd()
#and set the work directory in case we have moved around or opened another project
setwd("/Users/xanbjg/Documents/R/PIPT/Strains/")
dir()

#Load packages (based on Bengts recommendation)
#install.packages("gridExtra")
library(ggplot2)
library(tidyverse)
library(lubridate) # useful for working with dates
library(cowplot) # useul for combining multiple plots
library(scales)
library(ggthemes)
library(dplyr)
library(broom)
library(ggpubr)
library(gridExtra)
library(car)

#PCA packages
library(lattice)
library(devtools)
library(plyr)
library(scales)
library(grid)
library(ggbiplot)
#library(ggpmisc)

#Read and Tranform data####
#Lets read in the data
All <- read.csv2("Strain_summary.csv")
head(All)
dim(All)

#Change factors to numeric on approperiate data
cols <- c(7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)
All[,cols] = apply(All[,cols], 2, function(x) as.numeric(as.character(x)))
#All[,cols] <- as.numeric(as.character(unlist(All[,cols])))
#All[,cols] <- lapply(All[,cols], as.numeric)
head(All)

#PCA on QC and traits####

#Lets also subset only the continous variables (remove error colums too) and run a PCA
#Conver dataframe to matrix for heatmap
cols2 <- c(7, 8, 10, 12, 13, 15, 17, 19, 21, 27, 28, 29, 30, 31) #add 25, 26 when data available
AllM <- as.matrix(All)
class(AllM) <-"numeric"
AllMc <- subset.matrix(AllM, select = cols2)
dim(AllMc)
head(AllMc)

#Add strain header
ID <- All$Running.ID 
row.names(AllMc)=ID
head(AllMc)

#Vector for coloring PCA
P <- All$Population

#Run PCA on many factors
PCA <- prcomp((AllMc), scale = TRUE)
dev.off()
dev.new(width=5, height=5)
PCA_all <-ggbiplot(PCA, ellipse=TRUE, labels=ID, groups=P, alpha=0, varname.size = 3, labels.size=2) +
  #coord_cartesian(ylim=c(-0.01,0.31), expand = F) +
  #scale_color_manual(values=c("#000000", "#008000", "#99FF99", "#2846FF", "#99CCFF")) +
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
#theme(legend.position = c(0.995, 0.995),
#legend.justification = c("right", "top"),
#legend.box.just = "right",
#legend.margin = margin(2, 2, 2, 2),
#legend.box.background = element_rect(fill='white'),
#legend.background = element_blank(),
#legend.spacing.x=unit(0, "cm"),
#legend.spacing.y=unit(0, "cm"))

print(PCA_all)
dev.copy(pdf, "PCAall.pdf")
dev.off()

#The QC parameters seems to seperate the two populations
#Dpes Fv/Fm, pH, pre-density or wait time correlates with EC50, or is this only batch effect? 

#Does Fv/Fm correlate with EC50?

dev.off()
dev.new(width=5, height=5)
par(mfrow=c(2,2))

plot(All$EC50~All$FvFm, col = All$Population)
S1 <- lm(All$EC50~All$Biovolume)
summary(S1)
abline(S1)

plot(All$EC50~All$pH, col = All$Population)
S1 <- lm(All$EC50~All$pH)
summary(S1)
abline(S1)

plot(All$EC50~All$Pre_Density, col = All$Population)
S1 <- lm(All$EC50~All$Pre_Density)
summary(S1)
abline(S1)

plot(All$EC50~All$Wait_time, col = All$Population)
S1 <- lm(All$EC50~All$Wait_time)
summary(S1)
abline(S1)

dev.copy(pdf, "QCregressions.pdf")
dev.off()

#The regression between Wait-time and EC50 may be lost when poolinf populations
#Lets run it on seperated data using ggplot

dev.off()
dev.new(width=5, height=5)

WaitTime <- ggplot(data = All) +
  geom_point(mapping = aes(x = Wait_time, y = EC50, color = Population, fill = Population))+
  geom_smooth(mapping = aes(x = Wait_time, y = EC50, color = Population), method = "lm") +
  #coord_cartesian(ylim=c(-0.01,4), expand = F) +
  scale_color_manual(values=c("Black", "Red", "Green")) +
  scale_fill_manual(values=c("Black", "Red", "Green")) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines 
  panel_border(colour = "black", size = 1) + # and a border around each panel
  #labs (x="Year", y=expression(Max~bio-volume~(mm^{3}~L^{-1}))) +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=12))) +
  theme(axis.text=(element_text(size=12))) +
  theme(legend.text = element_text(face = "italic")) +
  #theme(aspect.ratio=1) +
  #stat_cor(aes(color = scientific_name), label.x = 3) +
  theme(legend.position = c(0.995, 0.995),
        legend.justification = c("right", "top"),
        #legend.box.just = "right",
        legend.margin = margin(2, 2, 2, 2),
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(WaitTime)
dev.copy(pdf, "Wait_time.pdf")
dev.off()

#So there is still a possitive slope but low statistical significance

#Trait based analysis####
#Lets do another PCA, reducing complexity to phenotypic trains

#Lets also subset only the continous variables (remove error colums too) and run a PCA
#Conver dataframe to matrix for heatmap
cols3 <- c(10, 17, 19, 21, 27, 28, 29, 30, 31) #add 25, 26 when data available
#AllM <- as.matrix(All)
class(AllM) <-"numeric"
AllMc2 <- subset.matrix(AllM, select = cols3)
dim(AllMc2)
head(AllMc2)

#Run PCA on traits
PCA2 <- prcomp((AllMc2), scale = TRUE)
dev.off()
dev.new(width=5, height=5)
PCA_traits <-ggbiplot(PCA2, ellipse=TRUE, labels=ID, groups=P, alpha=0, varname.size = 3, labels.size=2) +
  #coord_cartesian(ylim=c(-0.01,0.31), expand = F) +
  #scale_color_manual(values=c("#000000", "#008000", "#99FF99", "#2846FF", "#99CCFF")) +
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
  xlim(-3, 4) +
  ylim(-3, 4)
#theme(legend.position = c(0.995, 0.995),
#legend.justification = c("right", "top"),
#legend.box.just = "right",
#legend.margin = margin(2, 2, 2, 2),
#legend.box.background = element_rect(fill='white'),
#legend.background = element_blank(),
#legend.spacing.x=unit(0, "cm"),
#legend.spacing.y=unit(0, "cm"))

print(PCA_traits)
dev.copy(pdf, "PCAtraits.pdf")
dev.off()

#Okey so none of the morphological traits coencide with tolerance. 
#also the tolerant strain Vg1-2#15, #20, and #28 are all over the morpology spectrum (PC1)
#lets plot a few anyway


dev.off()
dev.new(width=5, height=5)
par(mfrow=c(2,2))

plot(All$EC50~All$Growth_Rate, col = All$Population)
S1 <- lm(All$EC50~All$Growth_Rate)
summary(S1)
abline(S1)

plot(All$EC50~All$Biovolume, col = All$Population)
S1 <- lm(All$EC50~All$Biovolume)
summary(S1)
abline(S1)

plot(All$EC50~All$Diamter_Length_ratio, col = All$Population)
S1 <- lm(All$EC50~All$Pre_Density)
summary(S1)
abline(S1)

plot(All$EC50~All$Surface_Vol_ratio, col = All$Population)
S1 <- lm(All$EC50~All$Surface_Vol_ratio)
summary(S1)
abline(S1)

dev.copy(pdf, "Traitregressions.pdf")
dev.off()

#No significant correlations. 
#Whats interesting here is that Vg1-2#15 has an odd morphology
#its wide and short and consequently has an abnormaly small surface to volume ratio. 
#But overall, this does not have an effect on Cu tolerance amongst the populations 

#Effects on growth rate####

#Does light intensity boost RO5s growth? No and its the same between experiments!
dev.off()
dev.new(width=5, height=5)
par(mfrow=c(2,2))
    
plot(All$RO5_Growth_Rate~All$Light_Intensity, col = All$Population)
S1 <- lm(All$RO5_Growth_Rate~All$Light_Intensity)
summary(S1)
abline(S1)

#Neither is all other strains affected
plot(All$Growth_Rate~All$Light_Intensity, col = All$Population)
S1 <- lm(All$Growth_Rate~All$Light_Intensity)
summary(S1)
abline(S1)

#There may be some Fv/Fm loss in some slow growing strains but relationship is full of outliers
plot(All$Growth_Rate~All$FvFm, col = All$Population)
S1 <- lm(All$Growth_Rate~All$FvFm)
summary(S1)
abline(S1)

#Interestingly its mainly small cells (<200 um3) that have low growth!
plot(All$Growth_Rate~All$Biovolume, col = All$Population)

dev.copy(pdf, "GrowthRate.pdf")
dev.off()

#Lets make a complex multiple regression trying to explain Growth rate (Fv/Fm important)
S2 <- lm(All$Growth_Rate~All$FvFm + All$Light_Intensity + All$EC50 + All$Biovolume)
summary(S2)

#Respons range of Cu####

#Compute respons range normalized to EC50
ResponsRange <- as.vector(All$EC95-All$EC05)/All$EC50
plot(ResponsRange)
plot(ResponsRange~All$EC05)
plot(ResponsRange~All$EC50)
plot(ResponsRange~All$Biovolume)
plot(ResponsRange~All$Growth_Rate)
plot(ResponsRange~All$Wait_time)
plot(ResponsRange~All$Surface_Vol_ratio)
plot(ResponsRange~All$FvFm)

#A long respons range seems to coencide with a low EC05 (i.e early onset of tox respons) 
#but nothing else really

#Population descriptors####

