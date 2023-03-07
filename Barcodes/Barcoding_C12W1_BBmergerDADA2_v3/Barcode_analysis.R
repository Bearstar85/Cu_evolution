#Houskeeping the loading packages####
#Lets clear old objects
rm(list=ls())
getwd()
#and set the work directory in case we have moved around or opened another project
setwd("~/Documents/R/Cu_evolution/Barcodes/Barcoding_C12W1_BBmergerDADA2_v3/C12W1_abundances")
dir()

#Bash commands for formation input data#####

#cd /Users/xanbjg/Documents/R/Cu_evolution/Barcodes/Barcoding_C12W1/C12W1_abundances
#rm *.Unknown_*
#for i in *.quantification*; do mv "$i" "${i/.quantification/}"; done

#Load packages and data (not adapted)
library(reshape2)
library(multtest) 
library(lattice)
library(devtools)
library(ggbiplot)
library(plyr)
library(scales)
library(ggplot2)
#library(devtools)
library(grid)
#library(gridExtra)
#library(lemon)
#library(gplots)
library(tidyverse)
#library(lubridate) # useful for working with dates
library(cowplot) # useul for combining multiple plots
library(ggthemes)
library(broom)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
#library(ggpmisc)
library(gplots)
#library(rstatix)
#install.packages("lemon")

#Data normalization####
#Make a vector list of all input file names
InputFiles <- list.files(path = "~/Documents/R/Cu_evolution/Barcodes/Barcoding_C12W1_BBmergerDADA2_v3/C12W1_abundances/Input")
InputFiles
#And read in the get Indexing file
#WARNING, certain Indexing names in Allels will not comply with differential equation function (notably x and df). 
Allels <- read.delim("../Allele_indexing.txt", sep = '\t', header=TRUE)
Indexing <- read.delim("../Indexing.txt", header=TRUE)
sapply(Allels, class)
#myData <- read.table(file = "Input/P21502_250.tsv", sep = '\t', header = TRUE)
#myData <- join(myData, Allels, by = "Barcode")

#First we setup the differential equation
solve_allele_equation <- function(data, n_obs, allele_col, equation) {
  
  # create a data.frame tol work on
  df <- data
  
  # assign the number of obsrevations to the correct allele
  for (j in 1:nrow(df)) {
    assign(x = df[[allele_col]][j], value = df[[n_obs]][j] )
  }
  
  # solve the equations
  x <- sapply(df[[equation]], function(x) eval(parse(text = x)) )
  names(x) <- NULL
  
  # add the equation solutions to the original data.frame
  df[["Total_C"]] <- x
  
  return(df)
  
}

# test the function on the test.data
#solve_allele_equation(data = Allels, 
                      #n_obs = "Total", 
                      #allele_col = "Homologs", 
                      #equation = "Differential_equation2")

#Lets loop the analysis
#First make a file to append to.
rm(All_reads)
All_reads <- data.frame(matrix(NA, ncol=26, nrow=0))[-1]

#i <-"P21502_146.tsv"
for (i in InputFiles) {
#Import data (too be looped)
input=paste("Input/",i,"", sep="")
myData <- read.table(file = input, sep = '\t', header = TRUE)

#Change Indexing of strains and allels
myData2 <- gsub('_C12W1_', ',', myData$Barcode)
myData <- cbind(myData,myData2)
myData <- myData %>% separate(myData2, c("ID.1", "allel"), sep = ",",remove = TRUE, convert = FALSE)
myData3 <- join(Allels, myData, by = "Barcode")
sapply(myData3, class)
#length(unique(myData3$Diffrential_ID))
#sort(unique(myData3$Diffrential_ID))

#Solve differential equation
myData4 <- solve_allele_equation(data = myData3, 
n_obs = "Total", 
allele_col = "Diffrential_ID", 
equation = "Differential_equation")

#The equation produces NA for 0/0, change to 0 again
myData4$Total_C[is.na(myData4$Total_C)] <- 0

#Generate sample ID for indexing experimental samples

Name <- i
ID <- gsub('.tsv', '', Name)
N <- nrow(myData4)
Sample <- rep(ID, each = N)
Sample <- as.data.frame(Sample)
colnames(Sample) <- c("ID")
Index <- join(Sample, Indexing, by = "ID")

#Fuse the sample indexing with data
myData5 <- cbind(Index, myData4)
#myData5 <- subset.data.frame(myData4, grepl("Yes", myData4$Keep_1copy))

#And normalize amplicon counts on trimmed data
Nreads <- sum(myData5$Total_C)
myData5$Relative_abundance <- (myData5$Total_C/Nreads)
Read_counts <- cbind(ID, Nreads)
sapply(myData5, class)

#Remove * on allel so it can be used for indexing, and second ID column
myData5$allel <- gsub('*', '', myData5$allel)
#myData5 <- subset(myData5, select = -c(1))

#Well save the header and export data
output <- paste("Indexed/",ID,".csv", sep="")
output2 <- paste("Read_counts/",ID,".csv", sep="")
output3 <- paste("Plots/",ID,".pdf", sep="")

header <- colnames(myData5)
write.table(header, "Indexed/1_header.tmp", sep=",",  col.names=FALSE, row.names =FALSE)
write.table(myData5, output, sep=",",  col.names=FALSE, row.names =FALSE)
write.table(Read_counts, output2, sep=",",  col.names=FALSE, row.names =FALSE)

#Append data to dataframe All_reads
All_reads <- rbind(All_reads, myData5)

#Okey lets also make some plots
Titel <- unique(myData5$Sample)

Fig1 <- ggplot(myData5, stats = "identity", aes(Strain, fill = allel)) +
  geom_col(data=myData5, aes(Strain, Relative_abundance)) +
  geom_hline(yintercept = 0.017) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = Titel) +
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
dev.copy(pdf, output3)
dev.off()

}
header <- read.csv("Indexed/1_header.tmp", header = F)
header <- t(header)
colnames(All_reads) <- (header)
write.table(All_reads, "Indexed/All_reads.csv", sep=",",  col.names=TRUE, row.names =FALSE)

#Count the occurrences of the wrong strain in the wrong population (False positives)
VGs <- subset.data.frame(All_reads, grepl('VG', All_reads$Strain))
VGs <- subset.data.frame(VGs, grepl('GP', VGs$Experiment))
VG_FP<- length(VGs$Population[VGs$Total_C > 0])

GPs <- subset.data.frame(All_reads, grepl('GP', All_reads$Strain))
GPs <- subset.data.frame(GPs, grepl('VG', GPs$Experiment))
GP_FP<- length(GPs$Population[GPs$Total_C > 0])

#The number of indisputable false positives is
VG_FP+GP_FP
#Out of
TO <- length(All_reads$Population[All_reads$Total_C > 0])
TO
#FDR is down to 0.7% of strains observations
((VG_FP+GP_FP)*2)/TO

#Graphical interpretation of data####
#Lets make a graph of only the MMs abundances.
MMsamples <- subset.data.frame(All_reads, grepl("0", All_reads$Timepoint))
MMsamples <- MMsamples[order(MMsamples$Strain),]
MMsamples$Relative_abundance[MMsamples$Relative_abundance == 0] <- NA
sapply(MMsamples, class)
AvaragesMM <- ddply(MMsamples, c("Experiment", "Strain", "allel"), summarise,
                     mean = mean(Relative_abundance), sd = sd(Relative_abundance))

#This way hard to make errorbars right.
Fig4 <- ggplot(AvaragesMM, stats = "identity", aes(Strain, fill = allel)) +
#ggplot(AvaragesMM, x = "Strain", y = "mean", add = "mean_sd", fill = aes(Strain, allel)) + #, fstats = "identity", aes(Strain, 
  geom_bar(data=AvaragesMM, aes(Strain, mean), position = "stack", stat = "identity") +
  facet_grid(rows = vars(Experiment)) +
  geom_hline(yintercept = 1/29) +
  #geom_errorbar(aes(ymin = (mean-sd), ymax = (mean+sd), width=.4, color = Metal), position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin = (mean-sd), ymax = (mean+sd), width=0.2), stat = "identity") + #position = position_dodge(0.3
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
  scale_fill_manual(values = c("forestgreen", "darkolivegreen2", "chartreuse", "blue3", "dodgerblue2", "deepskyblue1", "darkorchid1", "grey")) +
  labs (x="Strain", y="Fraction of amplicons") +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=6))) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.position= c(0.83, 0.82),
        #legend.box="horizontal",
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #scale_fill_manual(breaks=c(3)),
        legend.margin = margin(2, 2, 2, 2),
        #guides(color = FALSE), #,color = FALSE, size = FALSE,
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm")) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

print(Fig4)

dev.copy(pdf, "Plots/Fig4.pdf")
dev.off()

#This way hard to remove boarders around bars, and keep errorbbars.
Fig5 <- ggbarplot(MMsamples, x = "Strain", y = "Relative_abundance", add = "mean_sd", fill = "allel", size = 0,  facet.by = "Experiment") + #color=NA removes border but also error bars, ar
  #geom_col(data=AvaragesMM, aes(Strain, mean)) +
  #geom_bar(stat="identity",size=2) +
  facet_grid(rows = vars(Experiment)) +
  #coord_cartesian(ylim=c(0, 0.0029), expand = F) +
  geom_hline(yintercept = 1/29) +
  #geom_errorbar(aes(ymin = (mean-sd), ymax = (mean+sd), width=.4, color = Metal), position = position_dodge(0.1)) +
  #geom_errorbar(aes(ymin = (mean-sd), ymax = (mean+sd), width=0.5)) + #position = position_dodge(0.3
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  theme(plot.title = element_text(vjust = -5, hjust = 0.1)) +
  scale_fill_manual(values = c("forestgreen", "darkolivegreen2", "chartreuse", "blue3", "dodgerblue2", "deepskyblue1", "darkorchid1", "grey")) +
  labs (x="Strain", y="Fraction of amplicons") +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=6))) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6)) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=0.25) +
  theme(plot.margin=unit(c(0,1,0,0.2),"cm")) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.position= c(0.83, 0.82),
        #legend.box="horizontal",
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #scale_fill_manual(breaks=c(3)),
        legend.margin = margin(2, 2, 2, 2),
        #guides(color = FALSE), #,color = FALSE, size = FALSE,
        legend.box.background = element_rect(fill='white'),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm")) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

print(Fig5)

dev.copy(pdf, "Plots/Fig5.pdf")
dev.off()

#This code can be used to pull out any single strains alleles by name
subset.data.frame(MMsamples, grepl("GP2-4_52", MMsamples$Strain))

#Lets pick one strain at a time and plot them individualy (Then Loop this)
df_uniq <- unique(All_reads$Strain)
colourCount = length(unique(mtcars$hp))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
df_uniq

#Lets make an Dataframe to save all ratios in
#AlleleRatios <- as.data.frame("x", header('ID', "Barcode1_RA", "Barcode2_RA", "AllelRatio"))
AlleleRatios <- data.frame(ID=character(),
                           Barcode1=character(),
                           Barcode1_RA=numeric(),
                           Total1=numeric(),
                           Barcode2=character(),
                           Barcode2_RA=numeric(),
                           Total2=numeric(),
                           AllelRatio=numeric(),
                 stringsAsFactors=FALSE)

#sapply(OneStrain, class)
quartz()

#THIS LOOP KEEPS FAILING DUE TO GGPLOT ERRORWarningsS (in last FigStrains part)
#Stop loop before plots to finish rest of script


#i <- "GP2-4_44"
for (i in df_uniq) {

OneStrain <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain))
PIPT <- subset.data.frame(OneStrain, grepl(pattern = 'GP|VG', OneStrain$Experiment))
PIPT <- subset.data.frame(PIPT, grepl(pattern = 'GP|VG', PIPT$Population))
Strain2 <- subset.data.frame(OneStrain, grepl('Strain', OneStrain$Experiment))
output <- paste("Plots/Strains/",i,".pdf", sep="")
output2 <- paste("Plots/Strains/",i,"_Strains.pdf", sep="")

#Lets plot how the allel ratios vary with abundance. Should be no correlation and close to 1
OneStrain2 <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain))
Barcode1 <- subset.data.frame(OneStrain2, grepl("C12W1_1", OneStrain2$Barcode))
Barcode1 <- select(Barcode1, c('ID',  "Barcode", "Relative_abundance", "Total_C"))
Barcode2 <- subset.data.frame(OneStrain2, grepl("C12W1_2", OneStrain2$Barcode))
Barcode2 <- select(Barcode2, c('ID',  "Barcode", "Relative_abundance", "Total_C"))
OneStrain3 <- join(Barcode1, Barcode2, by = "ID")
colnames(OneStrain3) <- (c('ID', "Barcode1", "Barcode1_RA", "Total1", "Barcode2", "Barcode2_RA", "Total2"))
OneStrain3$AllelRatio <- OneStrain3$Barcode1_RA/OneStrain3$Barcode2_RA
OneStrain3$ReadSum <- OneStrain3$Total1+OneStrain3$Total2

#append results here
AlleleRatios <- rbind(AlleleRatios, OneStrain3)
#End below here to avoid crash on figures
#}

#This needs to be here for some annoying reason, loop terminates otherwise
quartz()
#sapply(PIPT, class)
PIPT$Replicate <- as.character(as.integer(PIPT$Replicate))
suppressWarnings(FigRA <- ggplot(PIPT, aes(Timepoint, Relative_abundance), message=FALSE) +
                   geom_point(mapping = aes(shape = Barcode, color = Replicate), size = 2) +
                   facet_wrap(~Experiment + Treatment) +
                   theme(strip.text.x = element_text(size=12, angle=0, face = "italic"),
                         strip.background = element_rect(colour="white", fill="white")) +
                   geom_path(mapping = aes(color = Replicate, linetype = Barcode), size = 0.5) +
                   #scale_color_brewer(palette = "Paired") +
                   #geom_errorbar(aes(xmin = (Low), xmax = (High), width=.4), position = position_dodge(0)) +
                   #geom_vline(xintercept = 0.86, color="red") +
                   #geom_vline(xintercept = 8,6, color="red", linetype = "dashed") +
                   #scale_shape_manual(values=c(19, 0, 1, 2)) +
                   scale_color_manual(values = getPalette(colourCount)) +
                   #scale_color_brewer(values = getPalette(colourCount)) +
                   #geom_line(mapping = aes(t_SM ~ Ratio, color = Model_SM)) +
                   #scale_linetype_manual(values = 2) +
                   coord_cartesian(ylim=c(0.0001, 1), xlim=c(-1, 44), expand = F) + #ylim=c(-0, 12.5)
                   #scale_y_discrete(limits=c(0,-3,-6,-9,-12,-15,-18)) +
                   scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                                 labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
                   #annotation_logticks(sides = "b", size = 0.1) + # adds linier tickmarks
                   #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
                   background_grid(major = "none", minor = "none") + # add thin horizontal lines
                   panel_border(colour = "black", size = 1) + # and a border around each panel
                   labs (x=("Time (days)"), y=("Relative abundance"), title = i) +
                   #theme(plot.title = element_text(vjust = - 13, hjust = 0.04)) +
                   theme(panel.spacing = unit(0.1, "lines")) +
                   theme(legend.title=element_text(size=12)) +
                   theme(legend.text=element_text(size=6)) +
                   theme(text=(element_text(size=12))) +
                   theme(axis.text=(element_text(size=12))) +
                   theme(panel.background = element_blank()) +
                   theme(legend.text = element_text(face = "italic")) +
                   theme(aspect.ratio=1))

suppressWarnings(print(FigRA))
#Sys.sleep(1)

suppressWarnings(dev.copy(pdf, output))
suppressWarnings(dev.off())

#End here to also make figures (sometimes crashes the loop unpredictably for me)
}

quartz()
#This figure is only relevant for plotting the occurences of alleles in the singe genotype samples (P21502_101-164)
suppressWarnings(FigStrains <- ggplot(Strain2, aes(Sample, Relative_abundance)) + #, message=FALSE
                   geom_point(mapping = aes(color = Barcode), size = 2, shape = 1) +
                   coord_cartesian(ylim=c(0.0001, 1), expand = F) + #ylim=c(-0, 12.5)
                   scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
                                 labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
                   background_grid(major = "none", minor = "none") + # add thin horizontal lines
                   panel_border(colour = "black", size = 1) + # and a border around each panel
                   labs (x=("Strain Genotype sample"), y=("Relative abundance"), title = i) +
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

#End here to also make genotype figures (crashes the loop unpredictably for me)
#}

#Allele ratios and correlations####

#sapply(AlleleRatios2, class)
#change 0 and inf to NA to fix plotting of ratios
AlleleRatios2 <- join(AlleleRatios, Indexing, by = "ID")
AlleleRatios2$Total1 <- as.numeric(as.integer(AlleleRatios2$Total1))
AlleleRatios2$Total2 <- as.numeric(as.integer(AlleleRatios2$Total2))
AlleleRatios2$Timepoint <- as.numeric(as.integer(AlleleRatios2$Timepoint))
AlleleRatios2$AllelRatio[AlleleRatios2$AllelRatio == 0] <- NA
AlleleRatios2$AllelRatio[AlleleRatios2$AllelRatio == "inf"] <- NA
AlleleRatios2$AllelRatio[sapply(AlleleRatios2$AllelRatio, is.infinite)] <- NA
AlleleRatios2 <- na.omit(AlleleRatios2)


#Lets look at the distribution of some values and parameters
min(AlleleRatios2$AllelRatio)
max(AlleleRatios2$AllelRatio)
hist(log10(AlleleRatios2$AllelRatio), breaks=80, main="AllelRatio") #xlim=c()
dev.copy(pdf, "Plots/Strains/HistRatios.pdf")
dev.off()

min(AlleleRatios2$Total1)
max(AlleleRatios2$Total1)

min(AlleleRatios2$Total2)
max(AlleleRatios2$Total2)

min(AlleleRatios2$ReadSum)
max(AlleleRatios2$ReadSum)
median(AlleleRatios2$ReadSum)
hist(log10(AlleleRatios2$ReadSum), breaks=80, main="BarcodeBoth_abundance") #xlim=c()
dev.copy(pdf, "Plots/Strains/Barcode1_RA_hist.pdf")
dev.off()

unique(AlleleRatios2$Population)

#May need to remove RO5AC samples for this
AlleleRatios2 <- subset.data.frame(AlleleRatios2, grepl(pattern = 'GP|VG', AlleleRatios2$Population))

colourCount = 6
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

FigRatio <- ggplot(AlleleRatios2, aes(Barcode1_RA, AllelRatio), message=FALSE) +
  geom_point(mapping = aes(color = interaction(Population,Timepoint,sep="-",lex.order=TRUE)), shape=1, size = 1) +
  facet_wrap(~Barcode1) +
  labs(colour="Population-Timepoint") +
  stat_smooth(mapping = aes(), size = 0.5, method = 'lm', formula = y ~ poly(x,2), se= TRUE) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  geom_hline(yintercept = 1) +
  scale_color_manual(values = c("firebrick", "orange", "gold", "blue", "cyan", "magenta")) +
  coord_cartesian(ylim=c(0.1, 10), xlim=c(0.0001, 1), expand = F) + #ylim=c(-0, 12.5)
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 2^x),
                labels = trans_format("log10", math_format(2^.x))) + #change the scale on y axis
  #annotation_logticks(sides = "b", size = 0.1) + # adds linier tickmarks
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x=("Relative abundance allel 1)"), y=("Ratio allel 1 to 2")) + #title = i
  #theme(plot.title = element_text(vjust = - 13, hjust = 0.04)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  #theme(legend.title=element_text(size=3)) +
  theme(legend.text=element_text(size=5)) +
  theme(text=(element_text(size=5))) +
  theme(axis.text=(element_text(size=5))) +
  theme(panel.background = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position="top")

print(FigRatio)

dev.copy(pdf, "Plots/Strains/FigRatios.pdf")
dev.off()

#This we can then do to remove 0+0 observations of allels for correlation (Outside loop)
AlleleRatios3 <- join(AlleleRatios, Indexing, by = "ID")
AlleleRatios3$ReadSum[AlleleRatios3$ReadSum == 0] <- NA
AlleleRatios3 <- AlleleRatios3 %>% drop_na(ReadSum)
AlleleRatios3$Total1 <- as.numeric(as.integer(AlleleRatios3$Total1))
AlleleRatios3$Total2 <- as.numeric(as.integer(AlleleRatios3$Total2))
sapply(AlleleRatios3, class)

#Remove Strain observations
AlleleRatios3 <- subset.data.frame(AlleleRatios3, grepl(pattern = 'GP|VG', AlleleRatios3$Experiment))

#Look at distribution
max(AlleleRatios3$Total1)
min(AlleleRatios3$Total2)
hist(log10(AlleleRatios3$Total1), breaks=80, main="Allel1_counts") #xlim=c()
dev.copy(pdf, "Plots/Strains/Allel1_counts.pdf")
dev.off()

hist(log10(AlleleRatios3$Total2), breaks=80, main="Allel2_counts") #xlim=c()
dev.copy(pdf, "Plots/Strains/Allel2_counts.pdf")
dev.off()

#Need to remove RO5AC samples
AlleleRatios3 <- subset.data.frame(AlleleRatios3, grepl(pattern = 'GP|VG', AlleleRatios3$Population))

#Look for a barcode
#OneStrain <- subset.data.frame(AlleleRatios3, grepl("P21502_101_C12W1_1*", AlleleRatios3$Barcode1))

StrainIndex <- subset(Allels, select = c(1,6))
colnames(StrainIndex) <- c("Barcode1", "Strain")
AlleleRatios3 <- inner_join(AlleleRatios3, StrainIndex, by = "Barcode1")
length(unique(AlleleRatios3$Strain))

FigReg <- ggplot(AlleleRatios3, aes(Total1+1, Total2+1), message=FALSE) +
  geom_point(mapping = aes(color = interaction(Population,Timepoint,sep="-",lex.order=TRUE), shape = Treatment), size = 1.2) +
  facet_wrap(~Strain) +
  labs(colour="Population-Timepoint") +
  stat_smooth(mapping = aes(), size = 0.5, method = 'lm', formula = y ~ poly(x,2), se= TRUE) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  geom_abline(intercept = 0, slope = 1) +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = c("black", "orange", "red", "blue", "cyan", "magenta")) +
  scale_shape_manual(values=c(1, 2)) +
  coord_cartesian(ylim=c(0.5, 28000), xlim=c(0.5, 28000), expand = F) + #ylim=c(-0, 12.5)
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
  #annotation_logticks(sides = "b", size = 0.1) + # adds linier tickmarks
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x=("Read count allele 1"), y=("Read count allele 2")) + #title = i
  #theme(plot.title = element_text(vjust = - 13, hjust = 0.04)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  #theme(legend.title=element_text(size=3)) +
  theme(legend.text=element_text(size=6)) +
  theme(text=(element_text(size=6))) +
  theme(axis.text=(element_text(size=6))) +
  theme(panel.background = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position="top")

print(FigReg)

dev.copy(pdf, "Plots/Strains/FigReg.pdf")
dev.off()

#Clustering of samples####

#We also need a dataframe for cbinding Matrix data
MatrixAllR <- data.frame(unique(All_reads$Sample))
colnames(MatrixAllR) <- (c('Sample'))
MatrixAllR
df_uniq

for (i in df_uniq) {
  #Lets also make a matrix format for further analysis
  #i <- "GP2-4_28"
  Reads <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain), select = c("Sample","Barcode", "Relative_abundance"))
  Barcodes <- unique(Reads$Barcode)
  Barcodes <- gsub("\\*", "", Barcodes)
  Barcodes <- gsub("\\?", "", Barcodes)
  
  Read1 <- subset.data.frame(Reads, grepl(paste(Barcodes[1]), Reads$Barcode), select = c("Sample","Relative_abundance"))
  colnames(Read1) <- (c("Sample", paste("",i,"_",paste(Barcodes[1],""), sep = "")))
  
  Read2 <- subset.data.frame(Reads, grepl(paste(Barcodes[2]), Reads$Barcode), select = c("Sample","Relative_abundance"))
  colnames(Read2) <- (c("Sample", paste("",i,"_",paste(Barcodes[2],""), sep = "")))
  Reads_matrix <- join(Read1, Read2, by = "Sample")
  
  Read3 <- subset.data.frame(Reads, grepl(paste(Barcodes[3]), Reads$Barcode), select = c("Sample","Relative_abundance"))
  colnames(Read3) <- (c("Sample", paste("",i,"_",paste(Barcodes[3],""), sep = "")))
  Reads_matrix <- join(Reads_matrix, Read3, by = "Sample")
  
  Read4 <- subset.data.frame(Reads, grepl(paste(Barcodes[4]), Reads$Barcode), select = c("Sample","Relative_abundance"))
  colnames(Read4) <- (c("Sample", paste("",i,"_",paste(Barcodes[4],""), sep = "")))
  Reads_matrix <- join(Reads_matrix, Read4, by = "Sample")
  
  #Remove NA columns since not all strains have 4 allels
  Reads_matrix <- Reads_matrix[ , apply(Reads_matrix, 2, function(x) !any(is.na(x)))]
  
  
  #Reads_matrix <- Reads_matrix %>%
  # select_if(~ !any(is.na(.)))
  #sapply(Reads_matrix, class)
  #Reads_matrix$paste("",i,"_",paste(Barcodes[4],"")) <- as.numeric(as.integer(Reads_matrix$paste("",i,"_",paste(Barcodes[4],""))))
  #AlleleRatios3$Total2 <- as.numeric(as.integer(AlleleRatios3$Total2))
  #Reads_matrix <- subset(Reads_matrix, select=c(is.na(colMeans(Reads_matrix))))
  
  #Then we just join this with growing matrix
  MatrixAllR <- join(MatrixAllR, Reads_matrix, by = "Sample")
  
}

sapply(MatrixAllR, class)

#Lets change the fucking header now
MatrixHeadR <- colnames(MatrixAllR)
MatrixHeadR <- gsub("_P21502_[0-9]{3}_C12W1", "", MatrixHeadR)
MatrixAll2R <- MatrixAllR
colnames(MatrixAll2R) <- MatrixHeadR


#Reformat to Matrix
V <- as.vector(MatrixAll2R$Sample)
MatrixAll2R <- as.matrix(MatrixAll2R)
class(MatrixAll2R) <-"numeric"
MatrixAll2R <- MatrixAll2R[,-1]
row.names(MatrixAll2R)=V
MatrixAll2R <- t(MatrixAll2R)

#I would like to blank out 0 values but i can't get the heatmap code to work then
#MatrixAll2R[MatrixAll2R == 0] <- NA

MatrixAll2R_cF <- scale(MatrixAll2R, scale = T, center = F) #center = F, scale = T)
MatrixAll2R_cT <- scale(MatrixAll2R, scale = T, center = T) #center = T, scale = T)
MatrixAll2R_log <- log10(MatrixAll2R+1)
MatrixAll2R_cFt <- t(MatrixAll2R_cF)
#MatrixAll2_cT <- t(MatrixAll2_cT)

#Make a heatmap without clustering
colRamp <- colorRampPalette(c("red", "black", "green"), space="rgb")(64)
heatmap(MatrixAll2R, col = colRamp, na.rm = TRUE, distfun = dist, cexRow=0.3, cexCol=0.5) #Rowv = NA, Colv = NA
dev.copy(pdf, "Plots/Strains/Heat1.pdf")
dev.off()

#2#Apply clustering (Eucladian distances)
#heatmap(MatrixAll2R, col = colRamp, distfun = dist, revC = T, #margins = c(8,6)
#        xlab = "Sample", ylab = "Allel", scale = "none")

#This prioritizes strain first
heatmap(MatrixAll2R_cT, col = colRamp, cexRow=0.3, cexCol=0.5, distfun = dist, revC = T, margins = c(8,6),
        xlab = "Sample", ylab = "Allel", scale = "none")
dev.copy(pdf, "Plots/Strains/Heat_cT_Eucl.pdf")
dev.off()

heatmap.2(MatrixAll2R, col = colRamp, scale="none", revC = F, margins = c(8,6),
          xlab = "Sample", ylab = "Allel", trace="none", cexRow=0.4)

breaks <- seq(min(MatrixAll2R, na.rm = T), max(MatrixAll2R, na.rm = T), length.out = 65)
heatmap.2(MatrixAll2R_cF, col = colRamp, scale="none", cexRow=0.3, cexCol=0.5, revC = F, margins = c(8,6),
          xlab = "Sample", ylab = "Allel", trace="none", na.color = "white", breaks=breaks)
dev.copy(pdf, "Plots/Strains/Heat_cT_Green_row.pdf")
dev.off()
#?heatmap.2
#symm=F, symkey=F,symbreaks=F

#Lets make a PCA and then add indexing
PCA1 <- prcomp(na.omit(MatrixAll2R_cT), center = TRUE, scale = TRUE)
Indexing2 <- as.data.frame(Barcode1$ID)
colnames(Indexing2) <- c("ID")
sapply(Indexing2, class)
Indexing3 <- join(Indexing2, Indexing, by = "ID")

PCA1
#Plot the eugenvalues of the PC components
plot(PCA1)
summary(PCA1)

#Export PCA results
PCAdata <- PCA1$rotation
class(PCAdata)
dim(PCAdata)
plot(PCAdata[,1])

#Add Population /species/Strain) vectors
PCAdata_df <- as.data.frame(cbind(PCAdata, Indexing3))
dim(PCAdata_df)
head(PCAdata_df)
#class(PCAdata_df$)
PCAdata_df$PC1 <- as.numeric(as.character(PCAdata_df$PC1))
PCAdata_df$PC2 <- as.numeric(as.character(PCAdata_df$PC2))
PCAdata_df$PC3 <- as.numeric(as.character(PCAdata_df$PC3))
unique(PCAdata_df$Population)
unique(PCAdata_df$Population)

#Make plots
plot(PCAdata_df$PC1, PCAdata_df$PC2, col = PCAdata_df$Timepoint, xlab = "PC1 (25%)", ylab = "PC2 (15%)",
     cex.lab = 1.3, cex.axis=1.3, pch = c(PCAdata_df$Population), cex =  1.5)
legend("bottomleft", cex = 1.2, legend=levels(PCAdata_df$Bottel), pch = c(1,2,3,4,5), col = c(1,2,3,4,5))
text(PCAdata, labels = ID , pos = 1, cex = 0.7)


PCA2 <- prcomp(na.omit(t(MatrixAll2R_cT)), scale = TRUE)
dim(PCA2)
head(PCA2)

ggbiplot(PCA1, ellipse=TRUE)
ggbiplot(PCA2, ellipse=TRUE)

#StrainCounts#####
#Now I will modify the analysis to make quantitative predictions about strain density using barcodes

#Lets loop the analysis
#First make a file to append to
rm(All_Strains)
All_Strains <- data.frame(matrix(NA, ncol=26, nrow=0))[-1]
header <- read.csv("Indexed/1_header.tmp", header = F)
header <- t(header)
colnames(All_reads) <- (header)

#i <-"P21502_146.tsv"
for (i in InputFiles) {
  #Import data (too be looped)
  input=paste("Input/",i,"", sep="")
  myData <- read.table(file = input, sep = '\t', header = TRUE)
  
  #Change Indexing of strains and allels
  myData2 <- gsub('_C12W1_', ',', myData$Barcode)
  myData <- cbind(myData,myData2)
  myData <- myData %>% separate(myData2, c("ID.1", "allel"), sep = ",",remove = TRUE, convert = FALSE)
  myData3 <- join(Allels, myData, by = "Barcode")
  sapply(myData3, class)
  #length(unique(myData3$Diffrential_ID))
  #sort(unique(myData3$Diffrential_ID))
  
  #Solve differential equation
  myData4 <- solve_allele_equation(data = myData3, 
                                   n_obs = "Total", 
                                   allele_col = "Diffrential_ID", 
                                   equation = "Differential_equation_strains")
  
  #The equation produces NA for 0/0, change to 0 again
  myData4$Total_C[is.na(myData4$Total_C)] <- 0
  
  #Now we can remove the second allels
  myData4 <- subset.data.frame(myData4, grepl('YES|NO', myData4$Modified_allelcounts))
  #Generate sample ID for indexing experimental samples
  
  Name <- i
  ID <- gsub('.tsv', '', Name)
  N <- nrow(myData4)
  Sample <- rep(ID, each = N)
  Sample <- as.data.frame(Sample)
  colnames(Sample) <- c("ID")
  Index <- join(Sample, Indexing, by = "ID")
  
  #Fuse the sample indexing with data
  myData5 <- cbind(Index, myData4)
  #myData5 <- subset.data.frame(myData4, grepl("Yes", myData4$Keep_1copy))
  
  #And normalize amplicon counts on trimmed data
  Nreads <- sum(myData5$Total_C)
  myData5$Relative_abundance <- (myData5$Total_C/Nreads)
  Read_counts <- cbind(ID, Nreads)
  sapply(myData5, class)
  
  #Remove * on allel so it can be used for indexing, and second ID column
  myData5$allel <- gsub('*', '', myData5$allel)
  #myData5 <- subset(myData5, select = -c(1))
  
  #Well save the header and export data
  header <- colnames(myData5)
  
  #Append data to dataframe All_reads
  All_Strains <- rbind(All_Strains, myData5)
  
  
}

#colnames(All_Strains) <- (header)
write.table(All_Strains, "Indexed/All_strains.txt", sep='\t',  col.names=TRUE, row.names =FALSE)


