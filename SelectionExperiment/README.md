2023-03-02
Author: Björn Andersson

## SelectionExperiment results

This analysis processes data and observations relative to the 42 day long copper tolerance selection experiment. In this experiment strains where pooled into GP (Reference inlet population) and VG (mining exposed population) origin, and subjected to selection in copper toxic (8.65 uM) and non-toxic environments (0.04 uM). In this directory, output data from **StrainPhenotyping** and **Barcodes** are combined with observational data from the selection experiment.


![Experimetal Design](https://github.com/Bearstar85/Cu_evolution/blob/master/ExperimentalDesign1.jpg)

## Input files and metadata

**All_strains.txt** [**All_strains_allreps.txt** same with contaminated samples retained]
This file is the intraspecific metabarcoding result file from the **Barcodes** analysis directory. It has the columns: *ID*; MiSeq sample indexing name, *Sample*; uniqe sample name for DNA source, *Experiment*; category of either singel genotype sample (strain), GP or VG selection experiment, *Population*; category of population that sample belongs to (GP or VG), *Treatment*; category of Control (non-toxic) or Copper (copper toxic) exposure, *Timepoint*; days since Treatment was started, *Bottel*; Uniqe experimental bottel identifier, *Replicate*; within Treatment X Population experimental replicate identifier (1-5), *Strain*; what strain the specific column pertains to quantify, *Total_C_*; raw number of amplicons (sum of both alleles) matching strain in sample, *Relative_abundance*; The relative abundance of Total_C in relation to all strains observed in the specific sample.   

*Barcode* through *Differential_ID* and *Index* through *allel*, not relevant for script.  


**Cu_predicted_rate.txt** and **Strain_Rate**
These two file contains mono-clonally estimated fitness proxies (growth rates) from the **StrainPhenotyping** analysis directory. **Strain_Rate** includes maximum growth rate (relevant for Control Treatment) wheras **Cu_predicted_rate.txt** DRC predicted growth rate estimates under 8.65 uM toxic Copper stress (relevant for the Copper treatment). It contain data on *Strain* name, *Growth_Rate* predicted under 8.65 uM Copper stress (day-1), and 95% confidence intervalls higher (*High*) and lower (*Low*) range.

**PAM_short**
This file contains observational data from a short 15 uM Copper assay employed during the the selection experiment. Strain RO5AC is included as a mono-clonal control. It indexes samples as in **All_strains.txt** based on *Population* (VG, GP or RO5AC) *Experiment* (GP or VG), Timepoint (days since *Treatment* (Control or Copper) was started), *Rep*licate number (1-5), *Exposure_h* is time since the 15 uM Copper assay was started (in hours continously measured), *SamplePoint* same as *Exposure_h* but in category format. *YII* the yield of photosystem II under 166 umol photons m-2 s-1 actinic light.

**PIPT_1value**
This file contains selected observational data pertaining to each individual *Bottel* Replicate. Col1-5 is indexing, as described above. *EC50* is the effective concentration of Copper that inhibits growth rate by 50%, after 42 days selection, with the *EC50SE* confidence (standard error of the mean). *Growth_rate9days* is the avarage growth rate of the bottel population during the first 9 days of *treatment* (same integrated timeperiod as for the first two metabarcoded strain abundance observations), *Generations_past* the populations avarage number of generations that past during 42 days (based on bulk RFU density observations), *Dominant_strain(%)*; qualitative description of the most abundant strain, and its relative abundance at day 42, *Second_strain*; other strains seen at >1% abundance, *BarcodeSampleContaminated_t42*; binary data if metabarcoded samples  aperared to be cross-contaminated by other samples (YES or NO).

**PIPT_EC50**
EC50 data across the selection experiment. This data is derived from dose-response data in the sub-directory **DRCs**. Col1-6 indexes data as described above. *EC50n* and *EC50SEn*; nominal EC50 value, *EC50a* and *EC50SEa* absolute EC50 values (uM Cu) corrected using ISP-MS quantification of Cu in media.  
 
**PIPT_wholeExp**
This file contains summary data from the selection experiment, in tall format. Col1-6 indexes as above. *Parameter* indexes which parameter is in *Value* ("Rate"=growth rated between consecutive 3 day periods (day-1), "RFU"= culture density in Relative Chl a fluorescence units, "FvFm"= maximum quantum yield of photosynthesis, and "pH"=pH)

**DRCs**
Directory containing dose-response data for the Selection experiment. This was ran iteratively and branched so that the two populations seperates first, then into diffrent timepoints, where data and scripts for the analysis are keept and can be run to process data. Note that on T30 DRC was collected for only 3 bottle replicates, but in addition to Copper, was also conducted against Cadmium, Nickel, and Lead (GP only). Metal concentrations are nominal values unless otherwise noted. The output data is summarized in the file **PIPT_EC50**.   

## Structure of script section by section

**SelectinExperiment.R**

#Houskeeping the loading packages####
load packages and set working directory to place with 3 subdirectories named Input, Plots, and Results. Put Input data in Input

#Read in and format data####
This part reads in the data for the strain quantification analysis, formats it, and makes a plot of it. 

#Phenotyping using metabarcoding####
This part of the script computes strain-specific growth rates from the metabarcoding data and the first 9 days of the selection experiments. These values are compared with the monoclonal observations. There are some important considerations regarding how to treat cases when a strain is not observed in the replicate. The script can be modified to set this to the detection limit of each sample, or ignore them (set as NA and filter out), or plot the raw data of individual replicates, rather than, computing a mean and conf. interval. At the moment it is set to ignore them, but changing -Inf ~ "NA", to "min (.[is.finite(.)]" on line 217 sets samples to detection limit instead. On line 238 avarages are computed from DummyDF, ingnoring NAs. If line 246 and 247 are activated, and line 249-250 deactivates, the Fig3B graphs below will plot individual replicates, rather than mean and 95% conf. intervals (just change *Barcode_mean* to *Barcode_rate* on line 292, turn on line 295 and off line 296-297). Optionally line 298 and 299 add regression analysis to plot. 

#PIPT Density, Fv/Fm, growth rate, RFU####
This part of the script makes graphs of control parameters during the selection experiment.

#Plot changes in EC50 during PIPT####
This part of the script plots changes in copper tolerance quantified either using DRC derived EC50s, or the Short PAM assay. 

#Strain abundances incl. contaminated samples####
This part reads in the **All_strains_allreps.txt** and makes plots of strain abundances in indvidual bottel reps and timepoints. It includes cross-contaminated samples removed from previous graphs and analyses. The code only retains the most abundant strains experiencing positive selection, and lumps the other together.