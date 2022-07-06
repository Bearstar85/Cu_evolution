#Houskeeping the loading packages####
#Lets clear old objects
rm(list=ls())
getwd()
#and set the work directory in case we have moved around or opened another project
setwd("~/Users/xanbjg/Documents/R/Cu_evolution/Genomics")
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

#lets also Remove RO5AC from samples
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

#Next step now is to do the same but for the barcode quantification of strains

#Predicting evolution of growth rate####

#Okey now we want to compute fitness for all strains using the barcoding quantification on day 9 of experiment

#Phenotyping using metabarcoding####

#Okey now we want to compute fitness for all strains using the barcoding quantification on day 9 of experiment

#We read in this data and join with barcodes data
PIPT_summary <- read.table(file = "Input/PIPT_1value.txt", sep = '\t', header = TRUE)
sapply(PIPT_summary, class)
sapply(Barcoding2, class)
unique(Barcoding2$Bottel)
unique(PIPT_summary$Bottel)

Barcoding_fitness <- join(Barcoding2, PIPT_summary, by = "Bottel", type = "left")
unique(Barcoding_fitness$Bottel)

Barcoding_fitness_GP <- subset.data.frame(Barcoding_fitness, grepl('GP', Barcoding_fitness$Population))
Barcoding_fitness_GP <- subset.data.frame(Barcoding_fitness_GP, grepl("GP2-4*", Barcoding_fitness_GP$Strain))
Barcoding_fitness_VG <- subset.data.frame(Barcoding_fitness, grepl('VG', Barcoding_fitness$Population))
Barcoding_fitness_VG <- subset.data.frame(Barcoding_fitness_VG, grepl("VG1-2*", Barcoding_fitness_VG$Strain))
Barcoding_fitness <- rbind(Barcoding_fitness_GP, Barcoding_fitness_VG)
#We loost the t0 data and need to get this back into all samples as the start conc.
#Lets use the avarages from the MM reps and selecting the right pop.
Time0 <- subset.data.frame(Barcoding_Avarages2, grepl('0', Barcoding_Avarages2$Timepoint))
Time0_GP <- subset.data.frame(Time0, grepl('GP', Time0$Population))
Time0_GP <- subset.data.frame(Time0_GP, grepl('Control', Time0_GP$Treatment))
Time0_GP <- subset.data.frame(Time0_GP, grepl("GP2-4*", Time0_GP$Strain))
Time0_VG <- subset.data.frame(Time0, grepl('VG', Time0$Population))
Time0_VG <- subset.data.frame(Time0_VG, grepl('Control', Time0_VG$Treatment))
Time0_VG <- subset.data.frame(Time0_VG, grepl("VG1-2*", Time0_VG$Strain))
Time0_2 <- rbind(Time0_GP, Time0_VG)

#Okey lets pick out the data we want to add
A <- Time0_2$mean
sd <- Time0_2$sd
S <- Time0_2$Strain
Time0 <- cbind(S, A, sd)
Time0 <- as.data.frame(Time0)
colnames(Time0) <- c("Strain", "t0Density", "t0Density_sd")

Barcoding_fitness2 <- join(Barcoding_fitness, Time0, by = "Strain", type = "left")

Barcoding_fitness3 <- subset.data.frame(Barcoding_fitness2, grepl("9", Barcoding_fitness2$Timepoint))
sapply(Barcoding_fitness3, class)


#Lets compute the growth rate
min(Barcoding_fitness3$Relative_abundance)
max(Barcoding_fitness3$Relative_abundance)
median(Barcoding_fitness3$Relative_abundance)
min(Barcoding_fitness3$t0Density)
max(Barcoding_fitness3$t0Density)
median(Barcoding_fitness3$t0Density)

#Firsts compute the relative growth rate per strain
Barcoding_fitness3$t0Density <- as.numeric(as.character(Barcoding_fitness3$t0Density))
Barcoding_fitness3$Relative_rate <- (log(Barcoding_fitness3$Relative_abundance/Barcoding_fitness3$t0Density)/9)
min(Barcoding_fitness3$Relative_rate)
max(Barcoding_fitness3$Relative_rate)
median(Barcoding_fitness3$Relative_rate)

#Then correct the observed bottle growth rate
Barcoding_fitness3$Barcode_rate <- Barcoding_fitness3$Growth_rate9days + Barcoding_fitness3$Relative_rate  
histogram(Barcoding_fitness3$Barcode_rate)

#Now we loop through and replace 0 with minimal value per replicate, to account for small sample sizes
#I use the NGS sample ID so there is no risk of mixing up samples/timepoints
f <- unique(Barcoding_fitness3$ID)
Header <- colnames(Barcoding_fitness3)
#DummyDF <- as.data.frame(t(Header), header=T)
#colnames(DummyDF) <- t(Header)
#DummyDF[0,1:42]
DummyDF <- data.frame(matrix(ncol = 42, nrow = 0))
DummyDF <- as.data.frame(DummyDF)
colnames(DummyDF) <- c(t(Header), "Observations")

f
for (i in f) {
  #i <- "P21502_279"
  
  #Select data
  myData <- subset.data.frame(Barcoding_fitness3, grepl(paste("\\b",i,"\\b", sep = ""), Barcoding_fitness3$ID))
  #Index if below detection limit
  myData$Observations <- (ifelse(myData$Barcode_rate > -Inf, ("Yes"), myData$Barcode_rate))
  myData$Observations <- (ifelse(myData$Observations == -Inf, ("No"), myData$Observations))
  #Change -Inf to min value across all columns
  myData2 <- myData %>%
    mutate(across(
      # For every column you want...
      c(Relative_rate, Barcode_rate),
      # ...change its values where appropriate:
      ~ case_when(
        # +Inf becomes the finite max.
        . ==  Inf ~ max(.[is.finite(.)]),
        # -Inf becomes the finite min.
        . == -Inf ~ min(.[is.finite(.)]),
        # Other values stay the same.
        TRUE      ~ .
      )
    ))
  
  #Append data to empty dataframe
  DummyDF <- rbind(DummyDF, myData2)
}

#Good, now we can sort out NULL observations using Observation column (or use the min detection limit in DummyDF)
DummyDF2 <- subset.data.frame(DummyDF, grepl('Yes', DummyDF$Observations))

#Lets extract averages per population below, There are options here!!!
#omitting this step will plot all datapoints instead (use Barcoding_fitness3)
#Use DummyDF2 to remove datapoints with 0 reads observed
#use DummyDF to replace 0 reads observations with the detection limit per sample

Barcoding_Fitness_Avarages <- ddply(DummyDF, c("Population", "Treatment", "Strain"), summarise,
                                    Barcode_mean = mean(Barcode_rate), Barcode_N = n(), 
                                    Barcode_High = mean(Barcode_rate)+2*sd(Barcode_rate)/sqrt(n()), 
                                    Barcode_Low = mean(Barcode_rate)-2*sd(Barcode_rate)/sqrt(n()))

#Here we can proceed using Barcoding_fitness3 for all data points, or Barcoding_Fitness_Avarages
#Barcoding_fitness3_Cu <- subset.data.frame(Barcoding_fitness3, grepl('Copper', Barcoding_fitness3$Treatment))
#Barcoding_fitness3_C <- subset.data.frame(Barcoding_fitness3, grepl('Control', Barcoding_fitness3$Treatment))
unique(Barcoding_Fitness_Avarages$Treatment)
Barcoding_fitness3_Cu <- subset.data.frame(Barcoding_Fitness_Avarages, grepl('Copper', Barcoding_Fitness_Avarages$Treatment))
Barcoding_fitness3_C <- subset.data.frame(Barcoding_Fitness_Avarages, grepl('Control', Barcoding_Fitness_Avarages$Treatment))


Observed_Rate <- read.table(file = "Input/Strain_Rate.txt", sep = '\t', header = TRUE)
sapply(Observed_Rate, class)
Cu_Rate <- read.table(file = "Input/Cu_predicted_rate.txt", sep = '\t', header = TRUE)
sapply(Observed_Rate, class)

Barcoding_fitness3_Cu <- join(Barcoding_fitness3_Cu, Cu_Rate, by = "Strain", type = "left")
Barcoding_fitness3_C <- join(Barcoding_fitness3_C, Observed_Rate, by = "Strain", type = "left")
dim(Barcoding_fitness3_Cu)
dim(Barcoding_fitness3_C)
Barcoding_fitness4 <- rbind(Barcoding_fitness3_C, Barcoding_fitness3_Cu)



#Need to clone start values to Cu
#Clone0 <- subset.data.frame(Barcoding_Avarages2, grepl('0', Barcoding_Avarages2$Timepoint))
#Clone0$Treatment <- recode_factor(Clone0$Treatment, "Control" = "Copper")
#Barcoding_Avarages2 <- rbind(Barcoding_Avarages2, Clone0)

#Need to add GP GP2-4_45dummy label to make color match with model, and RO5AC
#DummyStrain <- subset.data.frame(Barcoding_fitness4, grepl('GP2-4_45and46', Barcoding_fitness4$Strain))
#DummyStrain$Strain <- recode_factor(DummyStrain$Strain, "GP2-4_45and46" = "GP2-4_46dummy")
#DummyR05 <- subset.data.frame(Barcoding_fitness4, grepl('GP2-4_45and46', Barcoding_fitness4$Strain))
#DummyR05$Strain <- recode_factor(DummyR05$Strain, "GP2-4_45and46" = "RO5AC")
#sapply(DummyR05, class)

#lets delete data from RO5ac
#DummyR05[, 15:41][DummyR05[, 15:41] > 0] <- 0
#DummyR05[, 15:41][DummyR05[, 15:41] < 0] <- 0
#DummyR05 %>% replace(is.na(.), 0)

#And add datasets
#Barcoding_fitness5 <- rbind(Barcoding_fitness4, DummyStrain, DummyR05)
#unique(Barcoding_fitness5$Strain)

OneStrain <- subset.data.frame(Barcoding_fitness4, grepl('GP2-4_26', Barcoding_fitness4$Strain))
plot(OneStrain$Barcode_mean~OneStrain$Growth_Rate)
sapply(OneStrain, class)

#Compute mean and sd per treatment using diffrent phenotyping approches
Barcoding_fitness5 <- Barcoding_fitness4[complete.cases(Barcoding_fitness4[,8:10]), ]
Pop_Avarages_barcode <- ddply(Barcoding_fitness5, c("Population", "Treatment"), summarise,
                                    Barcode_Mean = mean(Barcode_mean), Barcode_N = n(), Barcode_sd = sd(Barcode_mean))  
Pop_Avarages_predicted <- ddply(Barcoding_fitness5, c("Population", "Treatment"), summarise,
                      Predicted_mean = mean(Growth_Rate), Barcode_N = n(), Predicted_sd = sd(Growth_Rate))  

Fig3B <- ggplot(Barcoding_fitness4, aes(Growth_Rate, Barcode_mean)) +
  facet_grid(rows = vars(Population), cols = vars(Treatment)) +
  geom_point(mapping = aes(color = Strain), size = 3, shape=1, stroke = 0.7) + #shape = Trea
  #geom_jitter(mapping = aes(color = Strain), size = 4, shape=1, width = 0.05) + #shape = Trea
  geom_errorbar(aes(ymin = Barcode_Low, ymax = Barcode_High, width=.1, color = Strain)) +
  geom_errorbarh(aes(xmin = Low, xmax = High, width=.1, color = Strain)) +
  #geom_line(mapping = aes(colour = Metal), size = 1) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  #geom_smooth(method='lm',formula=Barcoding_fitness5$Growth_rate~Barcoding_fitness5$Barcode_rate) +
  #stat_smooth(mapping = aes(), size = 0.5, method = 'lm', formula = y ~ x, se= TRUE) + 
  #scale_color_brewer(palette = "Paired") +
  facet_grid(rows = vars(Population), cols = vars(Treatment)) +
  geom_abline(intercept = 0, slope = 1) +
  #geom_errorbar(aes(ymin = (mean-sd), ymax = (mean+sd), width=.4, color = Metal), position = position_dodge(0.1)) +
  #scale_shape_manual(values=c(15, 16, 17, 18)) +
  scale_color_manual(values=c(col_vector)) + #aesthetics = c("colour", "fill")) +
  #geom_line(mapping = aes(t_SM ~ Ratio, color = Model_SM)) +
  #scale_linetype_manual(values = 2) +
  coord_cartesian(xlim=c(-0.3, 2), ylim=c(-0.3, 2), expand = F) + #ylim=c(-10000,+10000)
  scale_x_discrete(limits=c(0 ,0.5, 1, 1.5, 2)) +
  scale_y_discrete(limits=c(0, 0.5, 1, 1.5, 2)) +
  #annotation_logticks(sides = "l") + # adds linier tickmarks
  #coord_cartesian(ylim=c(-4,4), expand = F) + #changes the y axis
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  #scale_fill_manual(values=c("#008000", "#2846FF", "#000000")) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x=expression("Modeled growth rate "~(day^{-1})), y=expression("Observed growth rate "~(day^{-1})), title = "") +
  theme(plot.title = element_text(vjust = - 14, hjust = 0.04)) +
  theme(panel.spacing = unit(0.2, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=5)) +
  theme(text=(element_text(size=12))) +
  theme(axis.text=(element_text(size=12, colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(strip.background =element_rect(fill='white')) +
  #stat_cor(aes(color = scientific_name), label.x = 3)
  theme(legend.position = "right", 
        legend.spacing.x=unit(0.1, "cm"),
        legend.spacing.y=unit(0, "cm"),
        legend.box.background = element_rect(fill='white'))
  #guides(fill=guide_legend(ncol=3), color = guide_legend(override.aes = list(size = 0.5)))


print(Fig3B)
dev.copy(pdf, "Plots/Fig3b.pdf")
dev.off()

#Lets save this data to correlate genomic features against
write.table(Barcoding_fitness4, file = "Results/Fitness.txt", sep = '\t', col.names = TRUE)

#lets check if the histogram of fitness looks different from EC50s
#But the two experiments where not reproducible in terms of selection pressure by Cu so that needs to be considered!
Fig1A <- ggplot(Barcoding_fitness4, aes(x=Growth_Rate, fill=Population, color=Population)) +
  facet_grid(rows = vars(Treatment)) +
  geom_histogram(bins = 30, alpha=0.5, position="identity", aes(y = ..density..), color="black") + #aes(y = ..density..)
  geom_density(alpha=0.7) +
  scale_color_manual(values=c("gray14", "#B85633")) +
  scale_fill_manual(values=c("gray14", "#B85633")) +
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  #coord_cartesian(xlim=c(6, 10.1), ylim=c(0,2.5), expand = F) + #ylim=c(-10000,+10000)
  #scale_x_discrete(limits=c(6,7,8,9,10)) +
  labs(x=expression("Growth rate"~(day^{-1})), y = "Strain density", title = "A) Predicted") +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 2) + # and a border around each panel
  theme(plot.title = element_text(vjust = - 8, hjust = 0.04, size = 20)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=15)) +
  theme(text=(element_text(size=20))) +
  theme(axis.text=(element_text(size=20, colour = "Black"))) +
  #theme(axis.text=(element_text(colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(strip.background =element_rect(fill='white')) +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.75, 0.40),
        legend.margin = margin(10, 10, 10, 10),
        legend.box.background = element_rect(fill='white'),
        legend.box="horizontal",
        legend.key = element_rect(fill = 'white', color = 'white'),
        legend.background = element_blank(),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig1A)
dev.copy(pdf, "Plots/Fig1A.pdf")
dev.off()

Fig1B <- ggplot(Barcoding_fitness4, aes(x=Barcode_mean, fill=Population, color=Population)) +
  facet_grid(rows = vars(Treatment)) +
  geom_histogram(bins = 30, alpha=0.5, position="identity", aes(y = ..density..), color="black") + #aes(y = ..density..)
  geom_density(alpha=0.7) +
  scale_color_manual(values=c("gray14", "#B85633")) +
  scale_fill_manual(values=c("gray14", "#B85633")) +
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  #coord_cartesian(xlim=c(6, 10.1), ylim=c(0,2.5), expand = F) + #ylim=c(-10000,+10000)
  #scale_x_discrete(limits=c(6,7,8,9,10)) +
  labs(x=expression("Growth rate"~(day^{-1})), y = "Strain density", title = "B) Observed") +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 2) + # and a border around each panel
  theme(plot.title = element_text(vjust = - 8, hjust = 0.04, size = 20)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=15)) +
  theme(text=(element_text(size=20))) +
  theme(axis.text=(element_text(size=20, colour = "Black"))) +
  #theme(axis.text=(element_text(colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(strip.background =element_rect(fill='white')) +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.75, 0.40),
        legend.margin = margin(10, 10, 10, 10),
        legend.box.background = element_rect(fill='white'),
        legend.box="horizontal",
        legend.key = element_rect(fill = 'white', color = 'white'),
        legend.background = element_blank(),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig1B)
dev.copy(pdf, "Plots/Fig1B.pdf")
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

#Plot changes in EC50 during PIPT####
PIPT_EC50 <- read.table(file = "Input/PIPT_EC50.txt", sep = '\t', header = TRUE)
sapply(PIPT_EC50, class)
PIPT_EC50 <- PIPT_EC50[complete.cases(PIPT_EC50[,9]), ]
PIPT_EC50_mean <- ddply(PIPT_EC50, c("Population", "Treatment", "Timepoint", "Experiment"), summarise,
                        Mean = mean(EC50a), N = n(), 
                        High = mean(EC50a)+2*sd(EC50a)/sqrt(n()), 
                        Low = mean(EC50a)-2*sd(EC50a)/sqrt(n())) 

Fig2B <- ggplot(PIPT_EC50_mean, aes(Timepoint, Mean)) +
  facet_grid(cols = vars(Experiment)) + # , rows = vars(Parameter),  +
  geom_point(mapping = aes(color = interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 2, shape=16, stroke = 0.7) + #shape = Trea
  geom_errorbar(aes(ymin = Low, ymax = High, width=1, color = interaction(Treatment,Population,sep="-",lex.order=TRUE))) +
  geom_line(mapping = aes(colour = interaction(Treatment,Population,sep="-",lex.order=TRUE), linetype=interaction(Treatment,Population,sep="-",lex.order=TRUE)), size = 1) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  scale_color_manual(values=c("gray14", "#008000", "#B85633", "Blue", "dodgerblue1", "Blue")) + #aesthetics = c("colour", "fill")) +
  scale_linetype_manual(values = c(1,2,1,1,2,1)) +
  geom_vline(xintercept = 42.2, color="black", linetype = "dashed") +
  coord_cartesian(xlim=c(-1, 46), ylim=c(6, 11), expand = F) + #ylim=c(-10000,+10000)
  #scale_x_discrete(limits=c(-0.5, 0 ,0.5, 1, 1.5, 2)) +
  scale_y_discrete(limits=c(6, 7, 8, 9, 10)) +
  #annotation_logticks(sides = "l") + # adds linier tickmarks
  #coord_cartesian(ylim=c(-4,4), expand = F) + #changes the y axis
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x="Time (days)", y=expression("EC50 (" *mu ~ "M Cu)"), title = "B") +
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

print(Fig2B)
dev.copy(pdf, "Plots/Fig2B.pdf")
dev.off()

##################
