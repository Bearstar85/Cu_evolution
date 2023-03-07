2023-03-02
Author: Björn Andersson

## Strain phenotypes based mono-clonal cultures

This step of the analysis processes and plots phenotypic data in the form of dose-respons curves (DRC) to Cu, growth curve data, and microscopy observations. To execute script there needs to be a sub-directories named Input, Plots, and Results, below the home directory.

![Experimetal Design](https://github.com/Bearstar85/Cu_evolution/blob/master/ExperimentalDesign1.jpg)

## Input files and metadata

**StrainMeta.txt**
This file is metadata that is not used in R scripts. It has metadata on the strains (col *Strain*), WGS identifies (col *Seq_ID*) that is linking back to labels of DNA sample (col *DNA_ID*), and what there alias was in the experiments (col *Running_ID*). Also includer are *Site*, sediment *Core*, isotopic dating of sediment *Approx.Age*, what *Species* it is, the *Date_isolated* it was isolated, and who was the *Isolater*, what *Growth_Rate* (day-1), *Biovolume* um3, did it have when the selection experiment was seeded. *Comments* referes to qualitative observations DRC, growth or DNA extraction observations.

**Strain_summary.txt**
This file is used by scripts and contains detailed information about how it was treated, microscopy data, and quality control parameters. *Population*; (GP or VG), *Plate.ID*; an identifier of the plate where DRC was run, *Running_ID* is strains alias ID in the selection experiment. *Incubation Position* in the growth chamber, which could cause confunding effects due to diffrences in supplied *Light intensity* (umol PAR photons m-2, s-1). *Wait time (h)* between mixing copper media for DRC and adding of cells. *RO5_Growth_Rate*; the growth rate of strain RO5AC on the specific 24-well plate while *Growth rate (day-1)* refers to the mean growth rate of the 5 control treatments of the strain with *S.D.* standard deviation. *Ph* and maximum quantum yield of Photosystem II (*Fv/Fm*) of the inocoulum culure, *Start Density (RFU)* is the Relative Chl a flourescence value (mean across all strain 23 reps) at the start of the treatment (aim is 0.003, or ca. 1000 cells per mL/well), *Pre_Density*; RFU density of the pre-culture used to start DRC measurments and Selection experiment, *drc_Normalization*; Indicates if DRC data was normalized against the mean of the 5 control reps (Normal), or when a t-test showed suppression of growth in the controls, the 3 lowest DRC concentrations (Low), *Start_density*; the microscopically determined strain-specific cellular density going into the Selection experiment (cells mL-1),  *Relative_Start_Density*; same as *Start_density* but normalized to relative abundance, *Biovolume*, the avarage (N=20 observations) biovolume of a cell (um3),  *Surfacearea*, the avarage surface area of a cell (um2), *Surface_Vol_ratio*, surface to volume ratio of cell (um2 um-3)  *Diameter*, avarage diameter of cell (um),  *Diamter_Length_ratio*; avarage diameter to length ratio (um um-1).

**GP_lowCu_norm.csv and VG_lowCu_norm.csv**
This DRC data for all strains (*Concentraion* of Copper in uM versus *Inhibition* (fraction of growth rate with 1-inf set at 1). The latter data has been computed from *Growth rate* (day-1)* via *Inhibition RAW* (fraction of growth rate) which in these files have been modified to reduce the influence pf suppress growth in Control, which was assayed via t-tests comparing the 3 lowest Cu treatments growth with the Control. Whenever a strain had significantly higher growth rate in the low Cu reps, *Inhibition*  was normalized against the mean of these values, rather than the control. Other Metadata included *Local* where *strain*s were sampled from, *Date* when the DRC experiment was performed, the *Wait time (h)* between mixing copper media for DRC and adding of cells, *Incubation Position* in the growth chamber, which could cause confunding effects due to diffrences in supplied *Light intensity* (umol PAR photons m-2, s-1). *Control growth rate (day-1)* refers to the mean growth rate of the 5 control treatments of the strain with *S.D.* standard deviation, as well as the across-plate refrence control RO5AC strain which was cultured in one well on each 24-well plate. *Start Density (RFU)* is the Relative Chl a flourescence value (mean across all strain 23 reps) at the start of the treatment (aim is 0.003, or ca. 1000 cells per mL/well). *Well position* is the position of the specific sample on the plate (see DRC Plate design fig below), and *Replicate* is dilution step number from 1 (lowest) to 12 (highest), *Running_ID* is strains alias ID in the experiment.

![DRC Plate design](https://github.com/Bearstar85/Cu_evolution/blob/master/DRC_design.jpg)

**GrowthCurveExamples.csv**
This is growth curve data from the median, fastest, and slowest growing strain for each GP and VG population. Data collected daily (*Time*) over a period of 4 days. -3 corresponds to density during the pre-culture step, and first 0 sample the density of the culture on the day the selection experiment was started, and second 0, the concentration it was diluted to to start the experiment. *RFU* is mean of 5 well replicates, with Standard deviation in *SD*. 

## Structure of scripts section by section

**DoseRespons.R**

#Houskeeping the loading packages####
load packages and set working directory to place with 3 subdirectories named Input, Plots, and Results. Put Input data in Input

#Read in and format data####
This part reads in the data, formats it, translate numerica (theoretical)concentrations of Copper into absolute ones base on ISP-MS quantification using liniear regression CuA=0.7034CuN-0.1421. The loops on line 83-248 then uses the drc package to proccess the DRC data and generat output data in the form of Inhibition coefficents (EC) at certain levels (5, 50 and 95%), as well as predicting the Inhibition at the concentration of copper that was used during the selection experiment (dose=8.6504) uM, line 97 and 181. Does respons curve graphs is also created for each strain, with strain RO5AC as a refrence, and saved into Plots directory. Data is then compiled and saved into Results/ECvaluesStrains.txt.

#Graphs####
This part of the scrip does a more targeted iterative analysis of the data and makes relevant plots. First the median, most and least Cu tolerant strains of each population is identified and then plotted together. Some more population distribution plots are made then, plus some stats on that.   

#Model selection#####
This part of the script predicts the outcome of the selection experiment, based on monoclonal fitness observations from DRCs and Growth curves. It outputs Results/DRCpredictions.txt, and makes a deterministic model of the strain selection process. Plots are made.

