#Houskeeping the loading packages####
#Lets clear old objects
rm(list=ls())
getwd()
#and set the work directory in case we have moved around or opened another project
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
#Read in data####
Barcoding <- read.table(file = "Input/All_strains.txt", sep = '\t', header = TRUE)
sapply(Barcoding, class)

#Remove the strain samples
unique(Barcoding$Strain)
Barcoding2 <- subset.data.frame(Barcoding, grepl('VG|GP', Barcoding$Experiment))
Barcoding_Avarages <- ddply(Barcoding2, c("Population", "Treatment", "Timepoint", "Strain"), summarise,
                     mean = mean(Relative_abundance), sd = sd(Relative_abundance ))

#Remove RO5 samples
Barcoding_Avarages2 <- subset.data.frame(Barcoding_Avarages, grepl('VG|GP', Barcoding_Avarages$Population))

#lets also Remove RO5AC from individual samples
unique(Barcoding_Avarages2$Treatment)
Barcoding_Avarages2 <- subset.data.frame(Barcoding_Avarages2, !grepl('RO5AC', Barcoding_Avarages2$Strain))

#Need to clone start values to Cu
Clone0 <- subset.data.frame(Barcoding_Avarages2, grepl('0', Barcoding_Avarages2$Timepoint))
Clone0$Treatment <- recode_factor(Clone0$Treatment, "Control" = "Copper")
Barcoding_Avarages2 <- rbind(Barcoding_Avarages2, Clone0)

#Need to add GP GP2-4_45dummy label to make color match with model
#DummyStrain <- subset.data.frame(Barcoding_Avarages2, grepl('GP2-4_45and46', Barcoding_Avarages2$Strain))
#DummyStrain$Strain <- recode_factor(DummyStrain$Strain, "GP2-4_45and46" = "GP2-4_46dummy")
#DummyStrain$mean <- DummyStrain$mean*0
#DummyStrain$sd <- DummyStrain$sd*0
#Barcoding_Avarages2 <- rbind(Barcoding_Avarages2, DummyStrain)
##Strain selection graphs#####
#Mimicing the models result

#Lets first pedict pop evolution in growth rate
#StrainsO <- cbind(GP_density2, VG_density2)
#longDataO <-melt(StrainsO)

library(RColorBrewer)
n <- 58
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

#Lets make it all in one graph
sapply(Barcoding_Avarages2, class)
Barcoding_Avarages2$Timepoint <- as.numeric(as.character(Barcoding_Avarages2$Timepoint))

Fig2.D <- ggplot(data = Barcoding_Avarages2, aes(x = Timepoint, y = mean , fill=Strain)) + geom_area(position='fill') +
  #facet_wrap(~pop) +
  facet_grid(rows = vars(Population), cols = vars(Treatment)) +
  labs (x="Time (days)", y=("Relative biomass"), title = "") +
  theme(plot.title = element_text(vjust = - 8, hjust = 0.04)) +
  coord_cartesian(xlim=c(0, 42), ylim=c(0, 1), expand = F) + #ylim=c(-10000,+10000)
  scale_color_manual(values=c(col_vector), aesthetics = c("colour", "fill")) +
  scale_x_discrete(limits=c(10,20,30,40)) +
  scale_y_discrete(limits=c(0.2,0.4,0.6,0.8,1)) +
  panel_border(colour = "black", size = 1) +
  background_grid(major = "none", minor = "none") +# and a border around each panel
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=9)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=10, colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(strip.background =element_rect(fill='white'))+
  theme(aspect.ratio=1.1) +
  theme(legend.position = "right", 
        legend.spacing.x=unit(0.1, "cm"),
        legend.spacing.y=unit(0, "cm"),
        legend.box.background = element_rect(fill='white')) +
  guides(fill=guide_legend(ncol=3), color = guide_legend(override.aes = list(size = 0.1)))

Fig2.D

dev.copy(pdf, "Plots/Fig2.D.pdf")
dev.off()

#PIPT Density, Fv/Fm, growth rate, RFU####
PIPT_whole <- read.table(file = "Input/PIPT_wholeExp.txt", sep = '\t', header = TRUE)
sapply(PIPT_whole, class)

#Remove NA values
PIPT_whole <- PIPT_whole[complete.cases(PIPT_whole), ]

#Break-up the data and compute avarages and conf. intervals
P <- unique(PIPT_whole$Parameter)
P
Avarages_PIPT  <- NULL;
for (i in P) {
  #i <- "FvFm"
  #Select data
  mydata <- subset.data.frame(PIPT_whole, grepl(paste("\\b",i,"\\b", sep = ""), PIPT_whole$Parameter))
  #Compute means and conf. intervals
  mydata2 <- ddply(mydata, c("Experiment","Population", "Treatment", "Timepoint"), summarise,
                   Mean = mean(Value), N = n(), 
                   High = mean(Value)+2*sd(Value)/sqrt(n()), 
                   Low = mean(Value)-2*sd(Value)/sqrt(n()))
  length(mydata)
  #Add index back
  n <- nrow(mydata2)
  Parameters <- rep(i, each = n)
  mydata3 <- cbind(mydata2, Parameters)
  Avarages_PIPT <- rbind(mydata3, Avarages_PIPT)
}
unique(Avarages_PIPT$Parameter)
Rate <- subset.data.frame(Avarages_PIPT, grepl("Rate", Avarages_PIPT$Parameter))
RFU <- subset.data.frame(Avarages_PIPT, grepl("RFU", Avarages_PIPT$Parameter))
FvFm <- subset.data.frame(Avarages_PIPT, grepl("FvFm", Avarages_PIPT$Parameter))
pH <- subset.data.frame(Avarages_PIPT, grepl("pH", Avarages_PIPT$Parameter))

Fig2A <- ggplot(Rate, aes(Timepoint, Mean)) +
  facet_grid(cols = vars(Experiment)) + # , rows = vars(Parameter),  +
  geom_point(mapping = aes(color = interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 2, shape=16, stroke = 0.7) + #shape = Trea
  geom_errorbar(aes(ymin = Low, ymax = High, width=1, color = interaction(Treatment,Population,sep="-",lex.order=TRUE))) +
  geom_line(mapping = aes(colour = interaction(Treatment,Population,sep="-",lex.order=TRUE), linetype=interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 1) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  scale_color_manual(values=c("gray14", "#008000", "#B85633", "Blue", "dodgerblue1", "Blue")) + #aesthetics = c("colour", "fill")) +
  scale_linetype_manual(values = c(1,2,1,1,2,1)) +
  coord_cartesian(xlim=c(-1, 46), ylim=c(-0.5, 2.2), expand = F) + #ylim=c(-10000,+10000)
  #scale_x_discrete(limits=c(-0.5, 0 ,0.5, 1, 1.5, 2)) +
  scale_y_discrete(limits=c(-0.5, 0, 0.5, 1, 1.5, 2)) +
  #annotation_logticks(sides = "l") + # adds linier tickmarks
  #coord_cartesian(ylim=c(-4,4), expand = F) + #changes the y axis
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x="Time (days)", y=expression("Growth rate "~(day^{-1})), title = "") +
  theme(plot.title = element_text(vjust = - 35, hjust = 0.04)) +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=15)) +
  theme(text=(element_text(size=15))) +
  theme(axis.text=(element_text(size=15, colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(strip.background =element_rect(fill='white')) +
  theme(legend.position = "top", 
        legend.spacing.x=unit(0.1, "cm"),
        legend.spacing.y=unit(0, "cm"),
        legend.key.width = unit(1.5,"cm"),
        legend.box.background = element_rect(fill='white'))
#guides(fill=guide_legend(nrow = 3), color = guide_legend(override.aes = list(size = 1)))
#guides(fill=guide_legend(ncol=2), color = guide_legend(override.aes = list(length=1))) #color = guide_legend(override.aes = list(length=1))

print(Fig2A)
dev.copy(pdf, "Plots/Fig2A.pdf")
dev.off()

FigSx2 <- ggplot(FvFm, aes(Timepoint, Mean)) +
  facet_grid(cols = vars(Experiment)) + # , rows = vars(Parameter),  +
  geom_point(mapping = aes(color = interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 2, shape=16, stroke = 0.7) + #shape = Trea
  geom_errorbar(aes(ymin = Low, ymax = High, width=1, color = interaction(Treatment,Population,sep="-",lex.order=TRUE))) +
  geom_line(mapping = aes(colour = interaction(Treatment,Population,sep="-",lex.order=TRUE), linetype=interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 1) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  scale_color_manual(values=c("gray14", "#008000", "#B85633", "Blue", "dodgerblue1", "Blue")) + #aesthetics = c("colour", "fill")) +
  scale_linetype_manual(values = c(1,2,1,1,2,1)) +
  coord_cartesian(xlim=c(-1, 43), ylim=c(0, 1), expand = F) + #ylim=c(-10000,+10000)
  #scale_x_discrete(limits=c(-0.5, 0 ,0.5, 1, 1.5, 2)) +
  scale_y_discrete(limits=c(0, 0.2, 0.4, 0.6, 0.8)) +
  #annotation_logticks(sides = "l") + # adds linier tickmarks
  #coord_cartesian(ylim=c(-4,4), expand = F) + #changes the y axis
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x="Time (days)", y=expression("Fv/Fm"), title = "") +
  theme(plot.title = element_text(vjust = - 30, hjust = 0.04)) +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=12)) +
  theme(text=(element_text(size=12))) +
  theme(axis.text=(element_text(size=12, colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.5) +
  theme(strip.background =element_rect(fill='white')) +
  theme(legend.position = "top", 
        legend.spacing.x=unit(0.1, "cm"),
        legend.spacing.y=unit(0, "cm"),
        legend.key.width = unit(1.5,"cm"),
        legend.box.background = element_rect(fill='white'))
#guides(fill=guide_legend(nrow = 3), color = guide_legend(override.aes = list(size = 1)))
#guides(fill=guide_legend(ncol=2), color = guide_legend(override.aes = list(length=1))) #color = guide_legend(override.aes = list(length=1))

print(FigSx2)
dev.copy(pdf, "Plots/FigSx2.pdf")
dev.off()

FigSx3 <- ggplot(pH, aes(Timepoint, Mean)) +
  facet_grid(cols = vars(Experiment)) + # , rows = vars(Parameter),  +
  geom_point(mapping = aes(color = interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 2, shape=16, stroke = 0.7) + #shape = Trea
  geom_errorbar(aes(ymin = Low, ymax = High, width=1, color = interaction(Treatment,Population,sep="-",lex.order=TRUE))) +
  geom_line(mapping = aes(colour = interaction(Treatment,Population,sep="-",lex.order=TRUE), linetype=interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 1) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  scale_color_manual(values=c("gray14", "#008000", "#B85633", "Blue", "dodgerblue1", "Blue")) + #aesthetics = c("colour", "fill")) +
  scale_linetype_manual(values = c(1,2,1,1,2,1)) +
  coord_cartesian(xlim=c(-1, 43), ylim=c(7, 9), expand = F) + #ylim=c(-10000,+10000)
  #scale_x_discrete(limits=c(-0.5, 0 ,0.5, 1, 1.5, 2)) +
  scale_y_discrete(limits=c(7, 7.5, 8, 8.5, 9)) +
  #annotation_logticks(sides = "l") + # adds linier tickmarks
  #coord_cartesian(ylim=c(-4,4), expand = F) + #changes the y axis
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x="Time (days)", y=expression("pH"), title = "") +
  theme(plot.title = element_text(vjust = - 30, hjust = 0.04)) +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=12)) +
  theme(text=(element_text(size=12))) +
  theme(axis.text=(element_text(size=12, colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.5) +
  theme(strip.background =element_rect(fill='white')) +
  theme(legend.position = "top", 
        legend.spacing.x=unit(0.1, "cm"),
        legend.spacing.y=unit(0, "cm"),
        legend.key.width = unit(1.5,"cm"),
        legend.box.background = element_rect(fill='white'))
#guides(fill=guide_legend(nrow = 3), color = guide_legend(override.aes = list(size = 1)))
#guides(fill=guide_legend(ncol=2), color = guide_legend(override.aes = list(length=1))) #color = guide_legend(override.aes = list(length=1))

print(FigSx3)
dev.copy(pdf, "Plots/FigSx3.pdf")
dev.off()

FigSx4 <- ggplot(RFU, aes(Timepoint, Mean)) +
  facet_grid(cols = vars(Experiment)) + # , rows = vars(Parameter),  +
  geom_point(mapping = aes(color = interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 2, shape=16, stroke = 0.7) + #shape = Trea
  geom_errorbar(aes(ymin = Low, ymax = High, width=1, color = interaction(Treatment,Population,sep="-",lex.order=TRUE))) +
  geom_line(mapping = aes(colour = interaction(Treatment,Population,sep="-",lex.order=TRUE), linetype=interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 1) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  scale_color_manual(values=c("gray14", "#008000", "#B85633", "Blue", "dodgerblue1", "Blue")) + #aesthetics = c("colour", "fill")) +
  scale_linetype_manual(values = c(1,2,1,1,2,1)) +
  coord_cartesian(ylim=c(0.0003, 1), expand = F) + #ylim=c(-10000,+10000)
  #scale_x_discrete(limits=c(0,2,4,6,8,10)) +
  scale_y_continuous(trans = "log10") + #change the scale on y axis
  #annotation_logticks(sides = "l") + # adds linier tickmarks
  #coord_cartesian(ylim=c(-4,4), expand = F) + #changes the y axis
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x="Time (days)", y=expression("Density (RFU)"), title = "") +
  theme(plot.title = element_text(vjust = - 30, hjust = 0.04)) +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=12)) +
  theme(text=(element_text(size=12))) +
  theme(axis.text=(element_text(size=12, colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.5) +
  theme(strip.background =element_rect(fill='white')) +
  theme(legend.position = "top", 
        legend.spacing.x=unit(0.1, "cm"),
        legend.spacing.y=unit(0, "cm"),
        legend.key.width = unit(1.5,"cm"),
        legend.box.background = element_rect(fill='white'))
#guides(fill=guide_legend(nrow = 3), color = guide_legend(override.aes = list(size = 1)))
#guides(fill=guide_legend(ncol=2), color = guide_legend(override.aes = list(length=1))) #color = guide_legend(override.aes = list(length=1))

print(FigSx4)
dev.copy(pdf, "Plots/FigSx4.pdf")
dev.off()


##################
