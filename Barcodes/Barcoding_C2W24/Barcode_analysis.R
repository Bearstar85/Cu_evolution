#Lets clear old objects
rm(list=ls())
getwd()
#and set the work directory in case we have moved around or opened another project
setwd("~/Documents/R/Cu_evolution/Barcodes/Barcoding_C2W24/C2W24_abundancs")
dir()

#Bash commands for formation input data#####

#cd /Users/xanbjg/Documents/R/Cu_evolution/Barcodes/Barcoding_C2W24/C2W24_abundancs
#rm *.Unknown_*
#for i in *.Barcode_Quantities*; do mv "$i" "${i/.Barcode_Quantities/}"; done

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
library(RColorBrewer)
#library(ggpmisc)
library(gplots)


#Make a vector list of all input file names
InputFiles <- list.files(path = "~/Documents/R/Cu_evolution/Barcodes/Barcoding_C2W24/C2W24_abundancs/Input")
InputFiles
#And read in the get Indexing file
Allels <- read.delim("../Allele_indexing.txt", header=TRUE)
Indexing <- read.delim("../Indexing.txt", header=TRUE)
sapply(Allels, class)
#Plot how bad the overlapps are
AllelsUniqe <- subset.data.frame(Allels, grepl("Yes", Allels$Keep_1copy))
hist(AllelsUniqe$Strain_copies, xlim=c(0,7), breaks=7, main="Allele overlapp")
dev.copy(pdf, "Plots/AllelOveralpps.pdf")
dev.off()

myData <- read.table(file = "Input/P21502_139.tsv" , sep = '\t', header = TRUE)
#i <- "P21502_119.tsv"
#Lets loop the analysis

for (i in InputFiles) {
#Import data (too be looped)
input=paste("Input/",i,"", sep="")
myData <- read.table(file = input, sep = '\t', header = TRUE)
sapply(myData, class)

#Change Indexing of strains and allels
myData2 <- gsub('_C2W24_', ',', myData$Barcode)
myData <- cbind(myData,myData2)
myData <- myData %>% separate(myData2, c("ID", "allel"), sep = ",",remove = TRUE, convert = FALSE)
myData3 <- join(myData, Allels, by = "Barcode")

#Generate sample ID for indexing experimental samples

Name <- i
ID <- gsub('.tsv', '', Name)
N <- nrow(myData3)
Sample <- rep(ID, each = N)
Sample <- as.data.frame(Sample)
colnames(Sample) <- c("ID")
Index <- join(Sample, Indexing, by = "ID")

#Fuse the sample indexing with data
#as now this needs to be done seperately for VG and GP pop to get right strain ID
#Toggel Keep_VG vs Keep_GP to swich
#I am keeping as many allels as possible to detect false positives
myData4 <- cbind(Index, myData3)
myData5 <- subset.data.frame(myData4, grepl("Yes", myData4$Keep_VG))

#And normalize amplicon counts on trimmed data
Nreads <- sum(myData5$Total)
myData5$Relative_abundance <- (myData5$Total/Nreads)
Read_counts <- cbind(ID, Nreads)
sapply(myData5, class)

#Remove * on allel so it can be used for indexing
myData5$allel <- gsub('*', '', myData5$allel)

#Well save the header and export data
output <- paste("Indexed/",ID,".csv", sep="")
output2 <- paste("Read_counts/",ID,".csv", sep="")
output3 <- paste("Plots/",ID,".pdf", sep="")

header <- colnames(myData5)
write.table(header, "Indexed/1_header.csv", sep=",",  col.names=FALSE, row.names =FALSE)
write.table(myData5, output, sep=",",  col.names=FALSE, row.names =FALSE)
write.table(Read_counts, output2, sep=",",  col.names=FALSE, row.names =FALSE)

#Okey lets also make some plots
  
Fig1 <- ggplot(myData5, stats = "identity", aes(Strain, fill = interaction(allel,Strain_copies,sep="-",lex.order=TRUE))) +
  geom_col(data=myData5, aes(Strain, Relative_abundance)) +
  labs(colour="Population-Timepoint") +
  geom_hline(yintercept = 0.02) +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  #theme_classic() +
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (title = ID) +
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

#Summarize the results####

#Run this in UNIX command in Results
#cat *.csv > xxx.txt

#Indexed/GP_adapted/GP_summary.txt
#Indexed/VG_adapted/VG_summary.txt
#Read_counts/GP_adapted/GP_summary.txt
#Read_counts/VG_adapted/VG_summary.txt

#Should do something like this here instead
#temp = list.files(pattern="*.csv")
#myfiles = lapply(temp, read.delim)

GP <- read.csv("Indexed/GP_adapted/GP_summary.txt", header=FALSE)
VG <- read.csv("Indexed/VG_adapted/VG_summary.txt", header=FALSE)
Strain <- read.csv("Indexed/Strains_oneAllele/Strain.txt", header=FALSE)

All_reads <- rbind(GP, VG, Strain)
header <- read.csv("Indexed/1_header.csv", header = F)
header <- t(header)
colnames(All_reads) <- (header)
unique(All_reads$Experiment)
unique(All_reads$Treatment)
head(All_reads)
sapply(All_reads, class)
All_reads$Replicate <- as.character(as.integer(All_reads$Replicate))
unique(All_reads$Bottel)
#change 0 to NA to fix plotting
#All_reads$Relative_abundance[All_reads$Relative_abundance == 0] <- NA

#Lets pick one strain at a time (Then Loop this)
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

#THIS LOOP KEEPS FAILING DUE TO GGPLOT ERRORWarningsS (in last FigStrains part)
#Stop loop before plots to finish rest of script
df_uniq
for (i in df_uniq) {

OneStrain <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain))
PIPT <- subset.data.frame(OneStrain, grepl(pattern = 'GP|VG', OneStrain$Experiment))
PIPT <- subset.data.frame(PIPT, grepl(pattern = 'GP|VG', PIPT$Population))
Strain2 <- subset.data.frame(OneStrain, grepl('Strain', OneStrain$Experiment))
output <- paste("Plots/Strains/",i,".pdf", sep="")
output2 <- paste("Plots/Strains/",i,"_Strains.pdf", sep="")

#Lets plot how the allel ratios vary with abundance. Should be no correlation and close to 1
OneStrain2 <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain))
Barcode1 <- subset.data.frame(OneStrain2, grepl("C2W24_1", OneStrain2$Barcode))
Barcode1 <- select(Barcode1, c('ID',  "Barcode", "Relative_abundance", "Total"))
Barcode2 <- subset.data.frame(OneStrain2, grepl("C2W24_2", OneStrain2$Barcode))
Barcode2 <- select(Barcode2, c('ID',  "Barcode", "Relative_abundance", "Total"))
OneStrain3 <- join(Barcode1, Barcode2, by = "ID")
colnames(OneStrain3) <- (c('ID', "Barcode1", "Barcode1_RA", "Total1", "Barcode2", "Barcode2_RA", "Total2"))
OneStrain3$AllelRatio <- OneStrain3$Barcode1_RA/OneStrain3$Barcode2_RA
OneStrain3$ReadSum <- OneStrain3$Total1+OneStrain3$Total2

#append results here
AlleleRatios <- rbind(AlleleRatios, OneStrain3)



#This needs to be here for some annoying reason, loop terminates otherwise
quartz()

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
Sys.sleep(1)

suppressWarnings(dev.copy(pdf, output))
suppressWarnings(dev.off())
}

#This needs to be here for some annoying reason, loop terminates otherwise
quartz()

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


#Allele ratios and correlations####

sapply(AlleleRatios2, class)
#change 0 and inf to NA to fix plotting of ratios
AlleleRatios2 <- join(AlleleRatios, Indexing, by = "ID")
AlleleRatios2$Total1 <- as.numeric(as.integer(AlleleRatios2$Total1))
AlleleRatios2$Total2 <- as.numeric(as.integer(AlleleRatios2$Total2))
AlleleRatios2$Timepoint <- as.numeric(as.integer(AlleleRatios2$Timepoint))
AlleleRatios2$AllelRatio[AlleleRatios2$AllelRatio == 0] <- NA
AlleleRatios2$AllelRatio[AlleleRatios2$AllelRatio == "inf"] <- NA
AlleleRatios2$AllelRatio[sapply(AlleleRatios2$AllelRatio, is.infinite)] <- NA
AlleleRatios2 <- na.omit(AlleleRatios2)


#Lets look at the distribution of values
min(AlleleRatios2$AllelRatio)
max(AlleleRatios2$AllelRatio)
hist(log10(AlleleRatios2$AllelRatio), breaks=80, main="AllelRatio") #xlim=c()
dev.copy(pdf, "Plots/Strains/HistRatios.pdf")
dev.off()

min(AlleleRatios2$Barcode1_RA)
max(AlleleRatios2$Barcode1_RA)
hist(log10(AlleleRatios2$Barcode1_RA), breaks=80, main="Barcode1_Relative abundance") #xlim=c()
dev.copy(pdf, "Plots/Strains/Barcode1_RA_hist.pdf")
dev.off()


unique(AlleleRatios2$Timepoint)
colourCount = 7
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#Need to remove RO5AC
AlleleRatios2 <- subset.data.frame(AlleleRatios2, grepl(pattern = 'GP|VG', AlleleRatios2$Population))

FigRatio <- ggplot(AlleleRatios2, aes(Barcode1_RA, AllelRatio), message=FALSE) +
  geom_point(mapping = aes(color = interaction(Population,Timepoint,sep="-",lex.order=TRUE)), shape=1, size = 1) +
  facet_wrap(~Barcode1) +
  labs(colour="Population-Timepoint") +
  stat_smooth(mapping = aes(), size = 0.5, method = 'lm', formula = y ~ poly(x,2), se= TRUE) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  geom_hline(yintercept = 1) +
  scale_color_manual(values = c("firebrick", "orange", "gold", "blue", "cyan", "magenta", "black")) +
  coord_cartesian(ylim=c(0.01, 100000), xlim=c(0.0001, 1), expand = F) + #ylim=c(-0, 12.5)
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
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

#Look at distribution
max(AlleleRatios3$Total1)
max(AlleleRatios3$Total2)
hist(log10(AlleleRatios3$Total1), breaks=80, main="Allel1_counts") #xlim=c()
dev.copy(pdf, "Plots/Strains/Allel1_counts.pdf")
dev.off()

hist(log10(AlleleRatios3$Total2), breaks=80, main="Allel2_counts") #xlim=c()
dev.copy(pdf, "Plots/Strains/Allel2_counts.pdf")
dev.off()

#Need to remove RO5AC
AlleleRatios3 <- subset.data.frame(AlleleRatios3, grepl(pattern = 'GP|VG', AlleleRatios3$Population))

FigReg <- ggplot(AlleleRatios3, aes(Total1+1, Total2+1), message=FALSE) +
  geom_point(mapping = aes(color = interaction(Population,Timepoint,sep="-",lex.order=TRUE)), shape=1, size = 1) +
  facet_wrap(~Barcode1) +
  labs(colour="Population-Timepoint") +
  stat_smooth(mapping = aes(), size = 0.5, method = 'lm', formula = y ~ poly(x,2), se= TRUE) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  geom_abline(intercept = 0, slope = 1) +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = c("firebrick", "orange", "gold", "blue", "cyan", "magenta", "black")) +
  coord_cartesian(ylim=c(0.5, 28000), xlim=c(0.5, 28000), expand = F) + #ylim=c(-0, 12.5)
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + #change the scale on y axis
  #annotation_logticks(sides = "b", size = 0.1) + # adds linier tickmarks
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1) + # and a border around each panel
  labs (x=("Read count allele 1)"), y=("Read count allele 2")) + #title = i
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

print(FigReg)

dev.copy(pdf, "Plots/Strains/FigReg.pdf")
dev.off()

#Lets continue and make data into format for heatmaps, PCAs, and correlations matrixes
# I have code to convert if first column=allel and rest are samples

#We also need a dataframe for cbinding Matrix data
MatrixAll <- data.frame(unique(All_reads$ID))
colnames(MatrixAll) <- (c('ID'))

#df_uniq2 <- subset(df_uniq, grepl("GP",df_uniq))

for (i in df_uniq) {
  #Lets also make a matrix format for further analysis
  #i <- "GP2-4_28"
  Reads <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain), select = c("ID","Barcode", "Total"))
  Barcodes <- unique(Reads$Barcode)
  Barcodes <- gsub("\\*", "", Barcodes)
  Barcodes <- gsub("\\?", "", Barcodes)
  
  Read1 <- subset.data.frame(Reads, grepl(paste(Barcodes[1]), Reads$Barcode), select = c("ID","Total"))
  colnames(Read1) <- (c("ID", paste("",i,"_",paste(Barcodes[1],""), sep = "")))
  
  Read2 <- subset.data.frame(Reads, grepl(paste(Barcodes[2]), Reads$Barcode), select = c("ID","Total"))
  colnames(Read2) <- (c("ID", paste("",i,"_",paste(Barcodes[2],""), sep = "")))
  Reads_matrix <- join(Read1, Read2, by = "ID")
  
  #Read3 <- subset.data.frame(Reads, grepl(paste(Barcodes[3]), Reads$Barcode), select = c("ID","Total"))
  #colnames(Read3) <- (c("ID", paste("",i,"_",paste(Barcodes[3],""), sep = "")))
  #Reads_matrix <- join(Reads_matrix, Read3, by = "ID")
  
  #Read4 <- subset.data.frame(Reads, grepl(paste(Barcodes[4]), Reads$Barcode), select = c("ID","Total"))
  #colnames(Read4) <- (c("ID", paste("",i,"_",paste(Barcodes[4],""), sep = "")))
  #Reads_matrix <- join(Reads_matrix, Read4, by = "ID")
  
  #Remove NA columns since not all strains have 4 allels
  Reads_matrix <- Reads_matrix[ , apply(Reads_matrix, 2, function(x) !any(is.na(x)))]
  
  
  #Reads_matrix <- Reads_matrix %>%
  # select_if(~ !any(is.na(.)))
  #sapply(Reads_matrix, class)
  #Reads_matrix$paste("",i,"_",paste(Barcodes[4],"")) <- as.numeric(as.integer(Reads_matrix$paste("",i,"_",paste(Barcodes[4],""))))
  #AlleleRatios3$Total2 <- as.numeric(as.integer(AlleleRatios3$Total2))
  #Reads_matrix <- subset(Reads_matrix, select=c(is.na(colMeans(Reads_matrix))))
  
  #Then we just join this with growing matrix
  MatrixAll <- join(MatrixAll, Reads_matrix, by = "ID")

}

MatrixAll <- data.frame(unique(All_reads$ID))
colnames(MatrixAll) <- (c("ID"))

#df_uniq2 <- subset(df_uniq, grepl("GP",df_uniq))

for (i in df_uniq) {
  #Lets also make a matrix format for further analysis
  #i <- "GP2-4_28"
  Reads <- subset.data.frame(All_reads, grepl(paste("\\b",i,"\\b", sep = ""), All_reads$Strain), select = c("ID","Barcode", "Total"))
  Barcodes <- unique(Reads$Barcode)
  Barcodes <- gsub("\\*", "", Barcodes)
  Barcodes <- gsub("\\?", "", Barcodes)
  
  Read1 <- subset.data.frame(Reads, grepl(paste(Barcodes[1]), Reads$Barcode), select = c("ID","Total"))
  colnames(Read1) <- (c("ID", paste("",i,"_",paste(Barcodes[1],""), sep = "")))
  
  Read2 <- subset.data.frame(Reads, grepl(paste(Barcodes[2]), Reads$Barcode), select = c("ID","Total"))
  colnames(Read2) <- (c("ID", paste("",i,"_",paste(Barcodes[2],""), sep = "")))
  Reads_matrix <- join(Read1, Read2, by = "ID")
  
  Read3 <- subset.data.frame(Reads, grepl(paste(Barcodes[3]), Reads$Barcode), select = c("ID","Total"))
  colnames(Read3) <- (c("ID", paste("",i,"_",paste(Barcodes[3],""), sep = "")))
  Reads_matrix <- join(Reads_matrix, Read3, by = "ID")
  
  Read4 <- subset.data.frame(Reads, grepl(paste(Barcodes[4]), Reads$Barcode), select = c("ID","Total"))
  colnames(Read4) <- (c("ID", paste("",i,"_",paste(Barcodes[4],""), sep = "")))
  Reads_matrix <- join(Reads_matrix, Read4, by = "ID")
  
  #Remove NA columns since not all strains have 4 allels
  Reads_matrix <- Reads_matrix[ , apply(Reads_matrix, 2, function(x) !any(is.na(x)))]
  
  #Then we just join this with growing matrix
  MatrixAll <- join(MatrixAll, Reads_matrix, by = "ID")
  
}

sapply(MatrixAll, class)

#Lets change the fucking header now
MatrixHead <- colnames(MatrixAll)
MatrixHead <- gsub("_P21502_[0-9]{3}_C2W24", "", MatrixHead)
MatrixAll2 <- MatrixAll
colnames(MatrixAll2) <- MatrixHead


#Reformat to Matrix
V <- as.vector(MatrixAll2$ID)
MatrixAll2 <- as.matrix(MatrixAll2)
class(MatrixAll2) <-"numeric"
MatrixAll2 <- MatrixAll2[,-1]
row.names(MatrixAll2)=V
MatrixAll2 <- t(MatrixAll2)

MatrixAll2_cF <- scale(MatrixAll2, scale = T, center = F) #center = F, scale = T)
MatrixAll2_cT <- scale(MatrixAll2, scale = T, center = T) #center = T, scale = T)
#MatrixAll2_cF <- t(MatrixAll2_cF)
#MatrixAll2_cT <- t(MatrixAll2_cT)


#Make a heatmap without clustering
colRamp <- colorRampPalette(c("red", "black", "green"), space="rgb")(64)
heatmap(MatrixAll2, col = colRamp, Rowv = NA, Colv = NA, na.rm = TRUE, distfun = dist)

#2#Apply clustering (Eucladian distances)
heatmap.2(MatrixAll2, col = colRamp, distfun = dist, revC = T, #margins = c(8,6)
        xlab = "Sample", ylab = "Allel", scale = "none")

heatmap(MatrixAll2_cT, col = colRamp, distfun = dist, revC = T, margins = c(8,6),
        xlab = "Strain", ylab = "Metal (EC50)", scale = "none")


heatmap.2(MatrixAll2_cF, col = colRamp, scale="none", revC = T, margins = c(8,6),
          xlab = "Strain", ylab = "Metal (EC50)", trace="none")

heatmap.2(EC50_cT, col = colRamp, scale="none", revC = T, margins = c(8,6),
          xlab = "Strain", ylab = "Metal (EC50)", trace="none")
#symm=F, symkey=F,symbreaks=F


?prcomp
#Run the PCA analysis on the data (no scaling necessary in this case; change cT/cF)
PCA1 <- prcomp(EC50_cT, scale = FALSE)
summary(PCA1)
PCA1
#Plot the eugenvalues of the PC components
plot(PCA1)
summary(PCA1)

#Export PCA results
PCAdata <- PCA1$rotation
class(PCAdata)
PCAdata
plot(PCAdata[,1])

#Change to relative abundances#####

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
MatrixHeadR <- gsub("_P21502_[0-9]{3}_C2W24", "", MatrixHeadR)
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
?heatmap.2
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

