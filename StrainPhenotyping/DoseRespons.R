#Houskeeping the loading packages####
#Lets clear old objects
rm(list=ls())
getwd()
#and set the work directory in case we have moved around or opened another project
setwd("~/Users/xanbjg/Documents/R/Cu_evolution/StrainPhenotyping")
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


#first create individual .csv files based on .csv containing data of all strains
GP_Cunorm <- read.csv2("Input/GP_lowCu_norm.csv")
VG_Cunorm <- read.csv2("Input/VG_lowCu_norm.csv")#, header=FALSE
All <- rbind(GP_Cunorm, VG_Cunorm)
sapply(All, class)
All$Concentration <- as.numeric(as.character(All$Concentration))
All$Growth_rate_.day.1. <- as.numeric(as.character(All$Growth_rate_.day.1.))
All$Inhibition <- as.numeric(as.character(All$Inhibition))

#Transform Concentration from numerical to absolute Cu concentration based on ISP-MS measurments
plot(All$Concentration)
All$Concentration <- (ifelse(All$Concentration > 0, All$Concentration*0.7034-0.1421, All$Concentration))
plot(All$Concentration)

#i want to do the analysis seperate for the two pop. so that RO5 is included separately
GP <- subset.data.frame(All, grepl("GP", All$Local))
GPstrain<- unique(GP$strain)
VG <- subset.data.frame(All, grepl("VG", All$Local))
VGstrain<- unique(VG$strain)
#Lets make an Dataframe to save all ratios in
#AlleleRatios <- as.data.frame("x", header('ID', "Barcode1_RA", "Barcode2_RA", "AllelRatio"))
DRCs <- data.frame(Strain=character(),
                           Predict=numeric(),
                           Predict_low=numeric(),
                           Predict_high=numeric(),
                           EC05=numeric(),
                   EC05_SE=numeric(),
                   EC05_Low=numeric(),
                   EC05_high=numeric(),
                   EC50=numeric(),
                   EC50_SE=numeric(),
                   EC50_Low=numeric(),
                   EC50_high=numeric(),
                   EC95=numeric(),
                   EC95_SE=numeric(),
                   EC95_Low=numeric(),
                   EC95_high=numeric(),
                           stringsAsFactors=FALSE)

Header <- colnames(DRCs)
#i <- "GP2-4_26"

#Start with the GP anlysis

for (i in GPstrain) {
  #Subset strain data
  myData <- subset.data.frame(GP, grepl(paste("\\b",i,"\\b", sep = ""), GP$strain))
  Strain <- i
  
  #drm = dose-response model for error plots (change model as needed)
  fitData= subset(myData, myData$Concentration >=0)
  bestDoseresponse.WB = drm(Inhibition~Concentration, data = fitData, fct = W2.2())
  bestDoseresponse.WB #give you the intercepts
  a = coef(bestDoseresponse.WB)[1]
  b = coef(bestDoseresponse.WB)[2]
  
  #predict response at specific concentration dose(x) with 95% confidence
  Response= PR(bestDoseresponse.WB, c(8.6504))
  Response= predict(bestDoseresponse.WB, data.frame(dose=8.6504, CURVE= "1"), interval = "confidence")
  #Response= PR(bestDoseresponse.WB, c(8.6504))
  Response2 = t(Response)
  
  #Compute EC values
  EC05=ED(bestDoseresponse.WB, respLev = 0.05, type="absolute", interval="tfls") #to get ED help file put ? in front of ED and run. 
  EC50=ED(bestDoseresponse.WB, respLev = 0.5, type="absolute", interval="tfls")
  EC95=ED(bestDoseresponse.WB, respLev = 0.95, type="absolute", interval="tfls")
  
  #Add to DRCs dataframe
  myData2 <- cbind(Strain, Response2, EC05, EC50, EC95)
  rownames(myData2) <- NULL
  colnames(myData2) <- NULL
  DRCs <- rbind(DRCs, myData2)
  
  #Lets also make a plot of data, using RO5 as a Reference point
  myData3 = subset.data.frame(GP, grepl("RO5", GP$strain))
  
  #drm = dose-response model for error plots (change model as needed)
  fitData3= subset(myData3, myData3$Concentration >=0)
  bestDoseresponse_RO5.WB = drm(Inhibition~Concentration, data = fitData3, fct = W2.2())
  myData4=merge(myData,myData3, all= TRUE)
  
  #Make joint analysis for plot
  myData4=merge(myData,myData3, all= TRUE)
  fitData= subset(myData4, myData4$Concentration >=0)
  multi.m1 <- drm(Inhibition~Concentration, data = fitData, strain, fct = W2.2())

  
  #generate RGB code for Strain
  rgb(0, 128, 0, maxColorValue=255)
  
  #generate RGB code for RO5
  rgb(40, 70, 255, maxColorValue=255)
  
  #This is where the plots ends up
  output <- paste("Plots/",i,".pdf", sep="")
  
  #This is a stupid way of plotting but DRC output files are not working well
  dev.off()
  dev.new(width=2.5, height=2.5)
  par(mar=c(7,7,5,5))
  plot(multi.m1,
       type = "all",
       col=c("gray14", "#008000"), 
       broken = TRUE, bp = 3,   #this controls where the axis is cut
       bcontrol = list(factor = 1.2), #= "defult", style = "gap", width = "0.02"), #this controls how the cut is made (distance form origo)
       conName = "Control",
       cex =  2, cex.axis=2, lty=c(1, 2), lwd=1, #size of: symbols, axis numbers, lines type, line thickness,  stuff
       xlim = c(0,13), #Changes,  must be the same in the row below  
       ylim =c(-0.2,1.2),
       xt = c(0, 4, 5, 7, 9, 12), #Specify x numbers
       yt = c(-0.2, 0.2, 0.6, 1), #Specify y number
       xlab = expression("Cu concentration (" *mu ~ "M)"), line = 5, cex.lab =2,
       ylab = "Inhibition of growth rate",
       legend = FALSE)
  legend(3, 1.45, title="",
         c(Strain,"RO5"), col=c("gray14", "#008000"), lty=1:2, cex=1.75, lwd=1, pch=c(1,2), box.lty=0, x.intersp=0.25, y.intersp=0.3, 
         bg=NULL, seg.len=1)
       #legend = TRUE, legendPos = c(15, 1.4), legendText=c(Strain, "RO5"), cex.legend = 1.5, lwd.legend = 0.25, text.width = 1) #text.width=c(0.25,0.25,0.25),
       #legend = FALSE, legendPos = c(1.75, 1), legendText=c("1.5h", "3h", "5h"), cex.legend = 2, lwd.legend = 2,
       #main = title(main = "T0", adj = 0.05, line = -2.5, cex.main = 2))
  plot(bestDoseresponse_RO5.WB, col = "#008000", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
  plot(bestDoseresponse.WB, col = "gray14", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
  dev.copy(pdf, output)
  dev.off()
  
}

#VG analysis
for (i in VGstrain) {
  #Subset strain data
  myData <- subset.data.frame(VG, grepl(paste("\\b",i,"\\b", sep = ""), VG$strain))
  Strain <- i
  
  #drm = dose-response model for error plots (change model as needed)
  fitData= subset(myData, myData$Concentration >=0)
  bestDoseresponse.WB = drm(Inhibition~Concentration, data = fitData, fct = W2.2())
  bestDoseresponse.WB #give you the intercepts
  a = coef(bestDoseresponse.WB)[1]
  b = coef(bestDoseresponse.WB)[2]
  
  #predict response at specific concentration dose(x) with 95% confidence
  Response= PR(bestDoseresponse.WB, c(8.6504))
  Response= predict(bestDoseresponse.WB, data.frame(dose=8.6504, CURVE= "1"), interval = "confidence")
  #Response= PR(bestDoseresponse.WB, c(8.6504))
  Response2 = t(Response)
  
  #Compute EC values
  EC05=ED(bestDoseresponse.WB, respLev = 0.05, type="absolute", interval="tfls") #to get ED help file put ? in front of ED and run. 
  EC50=ED(bestDoseresponse.WB, respLev = 0.5, type="absolute", interval="tfls")
  EC95=ED(bestDoseresponse.WB, respLev = 0.95, type="absolute", interval="tfls")
  
  #Add to DRCs dataframe
  myData2 <- cbind(Strain, Response2, EC05, EC50, EC95)
  rownames(myData2) <- NULL
  colnames(myData2) <- NULL
  DRCs <- rbind(DRCs, myData2)
  
  #Lets also make a plot of data, using RO5 as a Reference point
  myData3 = subset.data.frame(VG, grepl("RO5", VG$strain))
  
  #drm = dose-response model for error plots (change model as needed)
  fitData3= subset(myData3, myData3$Concentration >=0)
  bestDoseresponse_RO5.WB = drm(Inhibition~Concentration, data = fitData3, fct = W2.2())
  myData4=merge(myData,myData3, all= TRUE)
  
  #Make joint analysis for plot
  myData4=merge(myData,myData3, all= TRUE)
  fitData= subset(myData4, myData4$Concentration >=0)
  multi.m1 <- drm(Inhibition~Concentration, data = fitData, strain, fct = W2.2())
  
  
  #generate RGB code for Strain
  rgb(0, 128, 0, maxColorValue=255)
  
  #generate RGB code for RO5
  rgb(40, 70, 255, maxColorValue=255)
  
  #This is where the plots ends up
  output <- paste("Plots/",i,".pdf", sep="")
  
  #This is a stupid way of plotting but DRC output files are not working well
  dev.off()
  dev.new(width=2.5, height=2.5)
  par(mar=c(7,7,5,5))
  plot(multi.m1,
       type = "all",
       col=c("#008000", "#B85633"), 
       broken = TRUE, bp = 3,   #this controls where the axis is cut
       bcontrol = list(factor = 1.2), #= "defult", style = "gap", width = "0.02"), #this controls how the cut is made (distance form origo)
       conName = "Control",
       cex =  2, cex.axis=2, lty=c(2, 1), pch=c(2,1), lwd=1, #size of: symbols, axis numbers, lines type, line thickness,  stuff
       xlim = c(0,13), #Changes,  must be the same in the row below  
       ylim =c(-0.2,1.2),
       xt = c(0, 4, 5, 7, 9, 12), #Specify x numbers
       yt = c(-0.2, 0.2, 0.6, 1), #Specify y number
       xlab = expression("Cu concentration (" *mu ~ "M)"), line = 5, cex.lab =2,
       ylab = "Inhibition of growth rate",
       legend = FALSE)
  legend(3, 1.45, title="",
         c(Strain,"RO5"), col=c("#B85633", "#008000"), lty=1:2, cex=1.75, lwd=1, pch=c(1,2), box.lty=0, x.intersp=0.25, y.intersp=0.3, 
         bg=NULL, seg.len=1)
  #legend = TRUE, legendPos = c(15, 1.4), legendText=c(Strain, "RO5"), cex.legend = 1.5, lwd.legend = 0.25, text.width = 1) #text.width=c(0.25,0.25,0.25),
  #legend = FALSE, legendPos = c(1.75, 1), legendText=c("1.5h", "3h", "5h"), cex.legend = 2, lwd.legend = 2,
  #main = title(main = "T0", adj = 0.05, line = -2.5, cex.main = 2))
  plot(bestDoseresponse_RO5.WB, col = "#008000", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
  plot(bestDoseresponse.WB, col = "#B85633", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
  dev.copy(pdf, output)
  dev.off()
  
}

colnames(DRCs) <- Header
DRCs2 <- DRCs
write.table(DRCs, file = "Results/ECvaluesStrains.txt", sep = '\t', col.names = TRUE)
#CAN START DOWNSTREAM ANALYSIS HERE
#DRCs2 <- read.table(file = "Results/ECvaluesStrains.txt", sep = '\t', header = TRUE)

DRCs2 <- DRCs2 %>% mutate_at(c(2:16), as.numeric)
sapply(DRCs2 , class)

#Graphs and modeling####
StrainMeta <- read.delim("Strain_summary.txt", sep = '\t', header=TRUE)
DRC_Meta <- join(DRCs2, StrainMeta, by = "Strain", type = "right")
sapply(DRC_Meta , class)

#Lets make a DRC curve with the median, min and max tolerant strains per pop +RO5 for References
GP_EC <- subset.data.frame(DRC_Meta, grepl("GP", DRC_Meta$Population))
#GP_EC <- GP_EC[order(EC50),]

GP_EC[which.max(GP_EC$EC50),]
GP_EC[which.min(GP_EC$EC50),]
median(GP_EC$EC50)


#We have GP2-4_56, GP2-4_28, and the median is 7.827672= between GP2-4_51/GP2-4_69
GP_x4 <- subset.data.frame(GP, grepl(pattern = 'GP2-4_56|GP2-4_28|GP2-4_51|RO5', GP$strain))
GP_56 <- subset.data.frame(GP, grepl("GP2-4_56", GP$strain))
GP_28 <- subset.data.frame(GP, grepl("GP2-4_28", GP$strain))
GP_51 <- subset.data.frame(GP, grepl("GP2-4_51", GP$strain))
RO5_GP <- subset.data.frame(GP, grepl("RO5", GP$strain))

#fitData= subset(GP_x4, GP_x4$Concentration >=0)
multi.m1 <- drm(Inhibition~Concentration, data = GP_x4, strain, fct = W2.2())
GP_56.WB = drm(Inhibition~Concentration, data = GP_56, fct = W2.2())
GP_28.WB = drm(Inhibition~Concentration, data = GP_28, fct = W2.2())
GP_51.WB = drm(Inhibition~Concentration, data = GP_51, fct = W2.2())
RO5_GP.WB = drm(Inhibition~Concentration, data = RO5_GP, fct = W2.2())

dev.off()
dev.new(width=2.5, height=2.5)
par(mar=c(7,7,5,5))
plot(multi.m1,
     type = "all",
     col=c("gray14", "gray14", "gray14", "#008000"), 
     broken = TRUE, bp = 3,   #this controls where the axis is cut
     bcontrol = list(factor = 1.2), #= "defult", style = "gap", width = "0.02"), #this controls how the cut is made (distance form origo)
     conName = "Control",
     cex =  2, cex.axis=2, lty=c(1, 1, 1, 2), pch=c(6, 1, 2, 4), lwd=2, #size of: symbols, axis numbers, lines type, line thickness,  stuff
     xlim = c(0,13), #Changes,  must be the same in the row below  
     ylim =c(-0.2,1.2),
     xt = c(0, 5, 6, 8, 10, 13), #Specify x numbers
     yt = c(-0.2, 0.2, 0.6, 1), #Specify y number
     xlab = expression("Cu concentration (" *mu ~ "M)"), line = 5, cex.lab =2,
     ylab = "Inhibition of growth rate",
     legend = FALSE)
box(lwd=2)
abline(v = median(GP_EC$EC50), col = "gray14", lwd = 2, lty = 1) #GPs median EC50
abline(v = 8.41, col = "#008000", lwd = 2, lty = 1) #RO5s EC50
legend(3, 1.2, title="", c("GP2-4_28","GP2-4_56", "GP2-4_61", "RO5AC"), col=c("gray14", "gray14", "gray14", "#008000"),
       cex=1.75, lwd=2, pch=c(6, 1, 2, 4), box.lty=0, x.intersp=0.25, y.intersp=0.4, 
       bg=NULL, seg.len=1)
plot(GP_56.WB, col = "gray14", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
plot(GP_28.WB, col = "gray14", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
plot(GP_51.WB, col = "gray14", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
plot(RO5_GP.WB, col = "#008000", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)


dev.copy(pdf, "Plots/Fig1A.pdf")
dev.off()

#Now for the VG strains
VG_EC <- subset.data.frame(DRC_Meta, grepl("VG", DRC_Meta$Population))
#VG_EC <- VG_EC[order(EC50),]

VG_EC[which.max(VG_EC$EC50),]
VG_EC[which.min(VG_EC$EC50),]
median(VG_EC$EC50)

#We have VG1-2_81, VG1-2_67, and the median is 7.954671= between VG1-2_88/VG1-2_86/
VG_x4 <- subset.data.frame(VG, grepl(pattern = 'VG1-2_81|VG1-2_67|VG1-2_86|RO5', VG$strain))
VG_81 <- subset.data.frame(VG, grepl("VG1-2_81", VG$strain))
VG_67 <- subset.data.frame(VG, grepl("VG1-2_67", VG$strain))
VG_86 <- subset.data.frame(VG, grepl("VG1-2_86", VG$strain))
RO5_VG <- subset.data.frame(VG, grepl("RO5", VG$strain))

#fitData= subset(VG_x4, VG_x4$Concentration >=0)
multi.m2 <- drm(Inhibition~Concentration, data = VG_x4, strain, fct = W2.2())
VG_81.WB = drm(Inhibition~Concentration, data = VG_81, fct = W2.2())
VG_67.WB = drm(Inhibition~Concentration, data = VG_67, fct = W2.2())
VG_86.WB = drm(Inhibition~Concentration, data = VG_86, fct = W2.2())
RO5_VG.WB = drm(Inhibition~Concentration, data = RO5_VG, fct = W2.2())

dev.off()
dev.new(width=2.5, height=2.5)
par(mar=c(7,7,5,5))
plot(multi.m2,
     type = "all",
     col=c("#B85633", "#B85633", "#B85633", "#008000"), 
     broken = TRUE, bp = 3,   #this controls where the axis is cut
     bcontrol = list(factor = 1.2), #= "defult", style = "gap", width = "0.02"), #this controls how the cut is made (distance form origo)
     conName = "Control",
     cex =  2, cex.axis=2, lty=c(1, 1, 1, 2), pch=c(6, 2, 1, 4), lwd=2, #size of: symbols, axis numbers, lines type, line thickness,  stuff
     xlim = c(0,13), #Changes,  must be the same in the row below  
     ylim =c(-0.2,1.2),
     xt = c(0, 5, 6, 8, 10, 13), #Specify x numbers
     yt = c(-0.2, 0.2, 0.6, 1), #Specify y number
     xlab = expression("Cu concentration (" *mu ~ "M)"), line = 5, cex.lab =2,
     ylab = "Inhibition of growth rate",
     legend = FALSE)
box(lwd=2)
abline(v = median(VG_EC$EC50), col = "#B85633", lwd = 2, lty = 1)
abline(v = 8.33, col = "#008000", lwd = 2, lty = 1)
legend(3, 1.2, title="", c("VG1-2_67", "VG1-2_81", "VG1-2_86", "RO5AC"), col=c("#B85633", "#B85633", "#B85633", "#008000"),
       cex=1.75, lwd=2, pch=c(6, 2, 1, 4), box.lty=0, x.intersp=0.25, y.intersp=0.4, 
       bg=NULL, seg.len=1)
plot(VG_81.WB, col = "#B85633", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
plot(VG_67.WB, col = "#B85633", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
plot(VG_86.WB, col = "#B85633", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)
plot(RO5_VG.WB, col = "#008000", type = "confidence", xlim = c(5,13), lty=c(0), add = TRUE)


dev.copy(pdf, "Plots/Fig1C.pdf")
dev.off()

#Good, now i want a histogram of both populations
Fig1B <- ggplot(DRC_Meta, aes(x=EC50, fill=Population, color=Population)) +
  geom_histogram(bins = 30, alpha=0.7, position="identity", color="black") + #aes(y = ..density..)
  #geom_density(alpha=0.7) +
  scale_color_manual(values=c("gray14", "#B85633"), labels = c("Reference inlet", "Mining inlet")) +
  scale_fill_manual(values=c("gray14", "#B85633"), labels = c("Reference inlet", "Mining inlet")) +
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  coord_cartesian(xlim=c(6, 10.1), ylim=c(0,8), expand = F) + #ylim=c(-10000,+10000)
  scale_x_discrete(limits=c(6,7,8,9,10)) +
  labs(x=expression("EC50 (" *mu ~ "M Cu)"), y = "Number of strains", title = "A") +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 2) + # and a border around each panel
  theme(plot.title = element_text(vjust = - 8, hjust = 0.04)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=30)) +
  theme(text=(element_text(size=27))) +
  theme(axis.text=(element_text(size=27, colour = "Black"))) +
  #theme(axis.text=(element_text(colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.55, 0.9),
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

#Lets do the scatter-plot version instead
#Need to index them from high to low EC50 first
DRC_Meta_GP <- subset.data.frame(DRC_Meta, grepl('GP', DRC_Meta$Population))
DRC_Meta_GP <- DRC_Meta_GP[order(DRC_Meta_GP$EC50),]
Rank <- seq(2, 29, by=1)
DRC_Meta_GP <- cbind(DRC_Meta_GP, Rank)

DRC_Meta_VG <- subset.data.frame(DRC_Meta, grepl('VG', DRC_Meta$Population))
DRC_Meta_VG <- DRC_Meta_VG[order(DRC_Meta_VG$EC50),]
Rank <- seq(1, 30, by=1)
DRC_Meta_VG <- cbind(DRC_Meta_VG, Rank)

#And join the datasets again
DRC_Meta2 <- rbind(DRC_Meta_GP, DRC_Meta_VG)
sapply(DRC_Meta2, class)
Fig1Bv2 <- ggplot(DRC_Meta2, aes(Rank, EC50)) +
  #geom_point(bins = 30, alpha=0.5, position="identity", aes(y = ..density..), color="black") + #aes(y = ..density..)
  geom_point(mapping = aes(color = Population, shape = Population), size = 5, stroke = 0.7) + # shape=c(16,17)
  geom_errorbar(aes(ymin = EC50_Low, ymax = EC50_high, width=.5, color = Population)) +
  scale_color_manual(values=c("gray14", "#B85633"), labels = c("Reference inlet", "Mining inlet")) +
  scale_fill_manual(values=c("gray14", "#B85633"), labels = c("Reference inlet", "Mining inlet")) +
  scale_shape_manual(values=c(16,17), labels = c("Reference inlet", "Mining inlet")) +
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  coord_cartesian(ylim=c(6, 10.5), xlim=c(0,31), expand = F) + #ylim=c(-10000,+10000)
  scale_x_discrete(limits=c(1,5,10,15,20,25,30)) +
  labs(x="Ranked number of strain", y =expression("EC50 (" *mu ~ "M Cu)"), title = "A") +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 2) + # and a border around each panel
  theme(plot.title = element_text(vjust = - 8, hjust = 0.04)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=30)) +
  theme(text=(element_text(size=27))) +
  theme(axis.text=(element_text(size=27, colour = "Black"))) +
  #theme(axis.text=(element_text(colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.55, 0.87),
        legend.margin = margin(10, 10, 10, 10),
        legend.box.background = element_rect(fill='white'),
        legend.box="horizontal",
        legend.key = element_rect(fill = 'white', color = 'white'),
        legend.background = element_blank(),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig1Bv2)
dev.copy(pdf, "Plots/Fig1v2B.pdf")
dev.off()

#Lets also make one with Growth_rate_.day.1.
sapply(DRC_Meta , class)
Fig1E <- ggplot(DRC_Meta, aes(x=Growth_Rate, fill=Population, color=Population)) +
  geom_histogram(bins = 30, alpha=0.5, position="identity", aes(y = ..density..), color="black") + #aes(y = ..density..)
  geom_density(alpha=0.7) +
  scale_color_manual(values=c("gray14", "#B85633"), labels = c("Reference inlet", "Mining inlet")) +
  scale_fill_manual(values=c("gray14", "#B85633"), , labels = c("Reference inlet", "Mining inlet")) +
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  coord_cartesian(xlim=c(0, 2), ylim=c(0,5), expand = F) + #ylim=c(-10000,+10000)
  scale_x_discrete(limits=c(0,0.5,1,1.5,2)) +
  labs(x=expression("Growth rate"~(day^{-1})), y = "Strain density", title = "E") +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 2) + # and a border around each panel
  theme(plot.title = element_text(vjust = - 8, hjust = 0.04)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=30)) +
  theme(text=(element_text(size=27))) +
  theme(axis.text=(element_text(size=27, colour = "Black"))) +
  #theme(axis.text=(element_text(colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.5, 0.85),
        legend.margin = margin(10, 10, 10, 10),
        legend.box.background = element_rect(fill='white'),
        legend.box="horizontal",
        legend.key = element_rect(fill = 'white', color = 'white'),
        legend.background = element_blank(),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig1E)
dev.copy(pdf, "Plots/Fig1E.pdf")
dev.off()
#Lets do the scatter-plot version instead
#Need to index them from high to low Growth_Rate first
DRC_Meta_GP <- subset.data.frame(DRC_Meta, grepl('GP', DRC_Meta$Population))
DRC_Meta_GP <- DRC_Meta_GP[order(DRC_Meta_GP$Growth_Rate),]
Rank <- seq(2, 29, by=1)
DRC_Meta_GP <- cbind(DRC_Meta_GP, Rank)

DRC_Meta_VG <- subset.data.frame(DRC_Meta, grepl('VG', DRC_Meta$Population))
DRC_Meta_VG <- DRC_Meta_VG[order(DRC_Meta_VG$Growth_Rate),]
Rank <- seq(1, 30, by=1)
DRC_Meta_VG <- cbind(DRC_Meta_VG, Rank)

#And join the datasets again
DRC_Meta2 <- rbind(DRC_Meta_GP, DRC_Meta_VG)
sapply(DRC_Meta2, class)

Fig1Cv2 <- ggplot(DRC_Meta2, aes(Rank, Growth_Rate)) +
  #geom_point(bins = 30, alpha=0.5, position="identity", aes(y = ..density..), color="black") + #aes(y = ..density..)
  geom_point(mapping = aes(color = Population, shape = Population), size = 5, stroke = 0.7) + # shape=c(16,17)
  geom_errorbar(aes(ymin = Growth_Rate-(2*SD_Grpwth_Rate./sqrt(5)), ymax = Growth_Rate+(2*SD_Grpwth_Rate./sqrt(5)), width=.5, color = Population)) +
  scale_color_manual(values=c("gray14", "#B85633"), labels = c("Reference inlet", "Mining inlet")) +
  scale_fill_manual(values=c("gray14", "#B85633"), labels = c("Reference inlet", "Mining inlet")) +
  scale_shape_manual(values=c(16,17), labels = c("Reference inlet", "Mining inlet")) +
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  coord_cartesian(ylim=c(0, 2), xlim=c(0,31), expand = F) + #ylim=c(-10000,+10000)
  scale_x_discrete(limits=c(1,5,10,15,20,25,30)) +
  labs(x="Ranked number of strain", y =expression("Growth rate"~(day^{-1})), title = "A") +
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 2) + # and a border around each panel
  theme(plot.title = element_text(vjust = - 8, hjust = 0.04)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=30)) +
  theme(text=(element_text(size=27))) +
  theme(axis.text=(element_text(size=27, colour = "Black"))) +
  #theme(axis.text=(element_text(colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.55, 0.27),
        legend.margin = margin(10, 10, 10, 10),
        legend.box.background = element_rect(fill='white'),
        legend.box="horizontal",
        legend.key = element_rect(fill = 'white', color = 'white'),
        legend.background = element_blank(),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig1Cv2)
dev.copy(pdf, "Plots/Fig1v2C.pdf")
dev.off()
#Are there statistical pop differences?
#Variance diffs EC? Yes: p-value = 0.0114
VarEC <- var.test(formula = EC50 ~ Population, data = DRC_Meta)
VarEC
#Variance diffs Rate? No: p-value = 0.472
VarRate <- var.test(formula = Growth_Rate ~ Population, data = DRC_Meta)
VarRate
#Differences in mean EC? No: p-value = 0.3473
t_EC <- t.test(formula = EC50 ~ Population, data = DRC_Meta, var.equal = FALSE)
t_EC
#Differences in mean EC? Yes: p-value = 0.02281
t_Rate <- t.test(formula = Growth_Rate ~ Population, data = DRC_Meta, var.equal = FALSE)
t_Rate

#Lets try to also plot growth curves with density plot

Density <- read.delim("Input/GrowthCurveExamples.txt", sep = '\t', header=TRUE)
sapply(Density , class)

GP_RFU <- subset.data.frame(Density, grepl("GP", Density$Population))
VG_RFU <- subset.data.frame(Density, grepl("VG", Density$Population))

  Fig1D <- ggplot(GP_RFU, aes(Time, RFU)) +
  geom_point(mapping = aes(shape = Strain, color = Strain), size = 6, stroke = 1.7) +
  geom_line(mapping = aes(colour = Strain, linetype=Strain), size = 1.5) + #Fits polynomal function to data (can be changed to lm: https://plotly.com/ggplot2/stat_smooth/) and https://stackoverflow.com/questions/31829528/specify-regression-line-intercept-r-ggplot2
  #scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin = (RFU-SD), ymax = (RFU+SD), width=1, color = Strain), position = position_dodge(0.1)) +
  scale_shape_manual(values=c(6, 2, 1, 4)) +
  scale_linetype_manual(values=c(1, 1, 1, 2)) +
  scale_color_manual(values=c("gray14", "gray14", "gray14", "#008000")) +
  #geom_line(mapping = aes(t_SM ~ Ratio, color = Model_SM)) +
  #scale_linetype_manual(values = 2) +
  coord_cartesian(xlim=c(-3.1, 4.1), ylim=c(0.001, 2), expand = F) + #ylim=c(-10000,+10000)
  scale_x_discrete(limits=c(-2,0,2,4)) +
  scale_y_continuous(trans = "log10") + #change the scale on y axis
  annotation_logticks(sides = "l") + # adds linier tickmarks
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 1.5) + # and a border around each panel
  labs (x="Time (days)", y=("Culture density (RFU)"), title = "D") +
    theme(plot.title = element_text(vjust = - 8, hjust = 0.04)) +
    theme(panel.spacing = unit(0.1, "lines")) +
    theme(legend.title=element_blank()) +
    theme(legend.text=element_text(size=25)) +
    theme(text=(element_text(size=27))) +
    theme(axis.text=(element_text(size=27, colour = "Black"))) +
    #theme(axis.text=(element_text(colour = "Black"))) +
    theme(panel.background = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(legend.text = element_text(face = "italic")) +
    theme(aspect.ratio=1) +
    theme(legend.position = c(0.65, 0.85),
          legend.margin = margin(10, 10, 10, 10),
          legend.box.background = element_rect(fill='white'),
          legend.box="horizontal",
          legend.key = element_rect(fill = 'white', color = 'white'),
          legend.background = element_blank(),
          legend.spacing.x=unit(0, "cm"),
          legend.spacing.y=unit(0, "cm"))

print(Fig1D)
dev.copy(pdf, "Plots/Fig1D.pdf")
dev.off()

Fig1F <- ggplot(VG_RFU, aes(Time, RFU)) +
  geom_point(mapping = aes(shape = Strain, color = Strain), size = 6, stroke = 1.7) +
  geom_line(mapping = aes(colour = Strain, linetype=Strain), size = 1.5) + 
  #scale_color_brewer(palette = "Paired") +
  geom_errorbar(aes(ymin = (RFU-SD), ymax = (RFU+SD), width=1, color = Strain), position = position_dodge(0.1)) +
  scale_shape_manual(values=c(4, 6, 2, 1)) +
  scale_linetype_manual(values=c(2, 1, 1, 1)) +
  scale_color_manual(values=c( "#008000","#B85633", "#B85633", "#B85633")) +
  #geom_line(mapping = aes(t_SM ~ Ratio, color = Model_SM)) +
  #scale_linetype_manual(values = 2) +
  coord_cartesian(xlim=c(-3.1, 4.1), ylim=c(0.001, 2), expand = F) + #ylim=c(-10000,+10000)
  scale_x_discrete(limits=c(-2,0,2,4)) +
  scale_y_continuous(trans = "log10") + #change the scale on y axis
  annotation_logticks(sides = "l") + # adds linier tickmarks
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #adjust text tilt and possition
  background_grid(major = "none", minor = "none") + # add thin horizontal lines
  panel_border(colour = "black", size = 2) + # and a border around each panel
  labs (x="Time (days)", y=("Culture density (RFU)"), title = "D") +
  theme(plot.title = element_text(vjust = - 8, hjust = 0.04)) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=25)) +
  theme(text=(element_text(size=27))) +
  theme(axis.text=(element_text(size=27, colour = "Black"))) +
  #theme(axis.text=(element_text(colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.60, 0.85),
        legend.margin = margin(10, 10, 10, 10),
        legend.box.background = element_rect(fill='white'),
        legend.box="horizontal",
        legend.key = element_rect(fill = 'white', color = 'white'),
        legend.background = element_blank(),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))

print(Fig1F)
dev.copy(pdf, "Plots/Fig1F.pdf")
dev.off()
  
#Model selection#####
#Lets  calculate the predicted growth rate of each strain
DRC_Meta$CuGrowth <- ((1-DRC_Meta$Predict)*(DRC_Meta$Growth_Rate))
DRC_Meta$CuGrowth_L <- ((1-DRC_Meta$Predict_low)*(DRC_Meta$Growth_Rate))
DRC_Meta$CuGrowth_H <- ((1-DRC_Meta$Predict_high)*(DRC_Meta$Growth_Rate))
sapply(DRC_Meta, class)

#Lets output this data for future use in PIPT experiment
write.table(DRC_Meta, file = "Results/DRCpredictions.txt", sep = '\t', col.names = TRUE)

#DRC_Meta$Start_density <- as.numeric(as.character(DRC_Meta$Start_density))
Strains <- unique(DRC_Meta$Strain)
Strains
P_density3 <- data.frame(Population=character(),
                    Strain=character(),
                   Treatment=character(),
                   t=numeric(),
                   Density=numeric(),
                   stringsAsFactors=FALSE)

#First we add the Cu models for the two pops
for (i in Strains) {
#i <- "GP2-4_26"
myData <- subset.data.frame(DRC_Meta, grepl(paste("\\b",i,"\\b", sep = ""), DRC_Meta$Strain))
pop <- myData$Population
  #First look at the imported data
  
  #create the model and input parameters, ie. a vector of (time=t), and set start density N0 (same as Exp ca 80000 um-2, mL-1 in surface area), and set Metal.
  t <- seq(0, 42, by=0.1)
  n <- length(t)
  u <-  myData$CuGrowth
  N0 <- myData$Start_density*100 #Specific start cell-density per strain and replicate bottle
  P_density <- N0*exp(t*u)
  P_density
  S <- rep(i, each = n)
  Treatment <- rep("Cu", each = n)
  P_density2 <- data.frame(pop, S, Treatment, t, P_density)
  P_density3 <- rbind(P_density3, P_density2)
}

#Now the control
for (i in Strains) {
  #i <- "GP2-4_26"
  myData <- subset.data.frame(DRC_Meta, grepl(paste("\\b",i,"\\b", sep = ""), DRC_Meta$Strain))
  pop <- myData$Population
  #First look at the imported data
  
  #create the model and input parameters, ie. a vector of (time=t), and set start density N0 (same as Exp ca 80000 um-2, mL-1 in surface area), and set Metal.
  t <- seq(0, 42, by=0.1)
  n <- length(t)
  u <-  myData$Growth_Rate
  N0 <- myData$Start_density*100 #Specific start cell-density per strain and replicate bottle
  P_density <- N0*exp(t*u)
  P_density
  S <- rep(i, each = n)
  Treatment <- rep("Control", each = n)
  P_density2 <- data.frame(pop, S, Treatment, t, P_density)
  P_density3 <- rbind(P_density3, P_density2)
}

sapply(P_density3, class)
unique(P_density3$pop)
unique(P_density3$Treatment)
#Lets work on modeling growth rate and changes in density of each strain


#Plot data to make sure it looks ok
ggplot(data = P_density3, aes(x = t, y = P_density)) +
  geom_point() +
  scale_y_continuous(trans = "log10") #change the scale on y axis

##Strain selection graphs#####

#Lets first pedict pop evolution in growth rate
#StrainsO <- cbind(GP_density2, VG_density2)
#longDataO <-melt(StrainsO)

library(RColorBrewer)
n <- 58
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

#Lets make it all in one graph

Fig2.C <- ggplot(data = P_density3, aes(x = t, y = P_density , fill=S)) + geom_area(position='fill') +
  #facet_wrap(~pop) +
  facet_grid(rows = vars(pop), cols = vars(Treatment)) +
  labs (x="Time (days)", y=("Relative biomass"), title = "") +
  theme(plot.title = element_text(vjust = - 8, hjust = 0.04)) +
  coord_cartesian(xlim=c(0, 42), ylim=c(0, 1), expand = F) + #ylim=c(-10000,+10000)
  scale_color_manual(values=c(col_vector), aesthetics = c("colour", "fill")) +
  scale_x_discrete(limits=c(10,20,30,40)) +
  scale_y_discrete(limits=c(0.2,0.4,0.6,0.8)) +
  panel_border(colour = "black", size = 1) +
  background_grid(major = "none", minor = "none") +# and a border around each panel
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=10)) +
  theme(text=(element_text(size=10))) +
  theme(axis.text=(element_text(size=10, colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(strip.background =element_rect(fill='white'))+
  theme(aspect.ratio=1.1) +
  theme(legend.position = "top", 
        legend.spacing.x=unit(0.2, "cm"),
        legend.spacing.y=unit(1, "cm"),
        legend.box.background = element_rect(fill='white')) +
  guides(fill=guide_legend(ncol=6), color = guide_legend(override.aes = list(size = 0.1)))

Fig2.C

dev.copy(pdf, "Plots/Fig2.C.pdf")
dev.off()

#Next step now is to do the same but for the barcode quantification of strains

################

#Predicting evolution of growth rate####

#Combine populations  biomass (cells rep)
GP_density <- subset.data.frame(P_density3, grepl("GP", P_density3$pop))
GP_density_Cu <- subset.data.frame(GP_density, grepl("Cu", GP_density$Treatment))
GP_density_C <- subset.data.frame(GP_density, grepl("Control", GP_density$Treatment))

VG_density <- subset.data.frame(P_density3, grepl("VG", P_density3$pop))
VG_density_Cu <- subset.data.frame(VG_density, grepl("Cu", VG_density$Treatment))
VG_density_C <- subset.data.frame(VG_density, grepl("Control", VG_density$Treatment))
dim(VG_density_Cu)
#Merge and summarize
GP_density_Cu2 <- t(dcast(data = GP_density_Cu,formula = S~t, fun.aggregate = sum,value.var = c("P_density"), head = FALSE))
colnames(GP_density_Cu2) <- as.character(GP_density_Cu2[1,])
GP_density_Cu3 <- GP_density_Cu2[-1, ]
class(GP_density_Cu3) <- "numeric"
dim(GP_density_Cu3)
GP_density_Cu3_v <- rowSums(GP_density_Cu3)

GP_density_C2 <- t(dcast(data = GP_density_C,formula = S~t, fun.aggregate = sum,value.var = c("P_density"), head = FALSE))
colnames(GP_density_C2) <- as.character(GP_density_C2[1,])
GP_density_C3 <- GP_density_C2[-1, ]
class(GP_density_C3) <- "numeric"
dim(GP_density_C3)
GP_density_C3_v <- rowSums(GP_density_C3)

VG_density_Cu2 <- t(dcast(data = VG_density_Cu,formula = S~t, fun.aggregate = sum,value.var = c("P_density"), head = FALSE))
colnames(VG_density_Cu2) <- as.character(VG_density_Cu2[1,])
VG_density_Cu3 <- VG_density_Cu2[-1, ]
class(VG_density_Cu3) <- "numeric"
dim(VG_density_Cu3)
VG_density_Cu3_v <- rowSums(VG_density_Cu3)

VG_density_C2 <- t(dcast(data = VG_density_C,formula = S~t, fun.aggregate = sum,value.var = c("P_density"), head = FALSE))
colnames(VG_density_C2) <- as.character(VG_density_C2[1,])
VG_density_C3 <- VG_density_C2[-1, ]
class(VG_density_C3) <- "numeric"
dim(VG_density_C3)
VG_density_C3_v <- rowSums(VG_density_C3)

#Lets output the growth rate estimates from the selection models so we can plot that for all metals later

#Sum strain Areal increase vs t here

#We have vectors for t and Area seperate in slection models so lets use that

#Drc models first, we can use dplyr's lag function in the equation
t <- seq(0, 42, by=0.1)
n <- length(t)
RateGPCu <- as.data.frame(cbind(t, rep("GP", each = n), rep("Cu", each = n), (log(GP_density_Cu3_v/lag(GP_density_Cu3_v))/(t - lag(t)))))
RateGPC <- as.data.frame(cbind(t, rep("GP", each = n), rep("Control", each = n), (log(GP_density_C3_v/lag(GP_density_C3_v))/(t - lag(t)))))
RateVGCu <- as.data.frame(cbind(t, rep("VG", each = n), rep("Cu", each = n), (log(VG_density_Cu3_v/lag(VG_density_Cu3_v))/(t - lag(t)))))
RateVGC <- as.data.frame(cbind(t, rep("VG", each = n), rep("Control", each = n), (log(VG_density_C3_v/lag(VG_density_C3_v))/(t - lag(t)))))
RateEvolution <- rbind(RateGPCu, RateGPC, RateVGCu, RateVGC)
colnames(RateEvolution) <- c("t", "Population", "Treatment", "Growth_rate")
sapply(RateEvolution, class)
RateEvolution$t <- as.numeric(as.character(RateEvolution$t))
RateEvolution$Growth_rate <- as.numeric(as.character(RateEvolution$Growth_rate))

FigSx.C <- ggplot(RateEvolution, aes(t, Growth_rate)) +
  geom_line(mapping = aes(color =Population, linetype=Treatment), size = 1.5) +
  #facet_wrap(~pop) +
  #facet_grid(rows = vars(pop), cols = vars(Treatment)) +
  labs (x="Time (days)", y=expression("Growth rate"~(day^{-1})), title = "") +
  #theme(plot.title = element_text(vjust = - 8, hjust = 0.04)) +
  coord_cartesian(xlim=c(0, 42), ylim=c(0, 2), expand = F) + #ylim=c(-10000,+10000)
  scale_color_manual(values=c("gray14", "#B85633"), labels = c("Reference inlet", "Mining inlet")) +
  scale_linetype_manual(values=c(1, 2), labels = c("Control", "Copper")) +
  scale_x_discrete(limits=c(10,20,30,40)) +
  scale_y_discrete(limits=c(0.5,1,1.5,2)) +
  panel_border(colour = "black", size = 2) +
  background_grid(major = "none", minor = "none") +# and a border around each panel
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=25)) +
  theme(text=(element_text(size=25))) +
  theme(axis.text=(element_text(size=25, colour = "Black"))) +
  theme(panel.background = element_blank()) +
  theme(legend.text = element_text(face = "italic")) +
  theme(strip.background =element_rect(fill='white'))+
  theme(aspect.ratio=1) +
  theme(legend.key.width = unit(2,"cm")) +
  theme(legend.position = c(0.50, 0.1),
        legend.margin = margin(1, 1, 1, 1),
        legend.box.background = element_rect(fill='white'),
        legend.box="horizontal",
        legend.key = element_rect(fill = 'white', color = 'white'),
        legend.background = element_blank(),
        legend.spacing.x=unit(0, "cm"),
        legend.spacing.y=unit(0, "cm"))
  #guides(linetype = guide_legend(override.aes = list(size = 10))))
  #guides(fill=guide_legend(ncol=6), color = guide_legend(override.aes = list(size = 0.1)))

FigSx.C

dev.copy(pdf, "Plots/FigSx.pdf")
dev.off()
##################

#Stats fro MAnuscript

t_EC
t_Rate
