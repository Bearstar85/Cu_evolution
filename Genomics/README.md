2023-03-02
Author: Björn Andersson

## SelectionExperiment results

This analysis looks for correlative patterns between phenotypic (derived from **SelectionExperiment** directory) and genomics copy-number variance traits in genes involved in Copper metabolism (based on WGS data). 


![Experimetal Design](https://github.com/Bearstar85/Cu_evolution/blob/master/ExperimentalDesign1.jpg)

## Input files and metadata

**CuGeneCoverage/P14203_101_CuCoverage.tsv**, etc
These files contains avarage WGS coverage depth for selected gene models. Each file corresponds to one strains genome sequencing data. The *Gene* contains the name of the  genemodel, *Model mean*; the avarage coverage across the model, *Model median*; media, *Contig mean*/*Contig median*; the mean coverage across the contig where the gene modelel is locates. The genomic avarage and mean is on row 2.

**CuMetabolism_SM_short.txt**
This file has detailed data of the manual annotation process of the gene models involved in copper metabolism. The annotation was primarily based on literature/genome search of *Thalassiosira pseudonana* copper meatbolism, and BLASTn/p search of homologs in *S. marinoi*, but also some GO-term searches and annotations from *T. oceanica*. *Gene*; contains the name of the  genemodel, *Putative_identity*; the annotated function of the protein with *Putative_identity_INDEX* having added the contig location of each model so homologs on diffrent chromosomes can be distinguished between, *Ortholog*; if there is sufficent *support* or evidence (Nucleotide, aa, and functional domain homology) the corresponding *T. pseudonan* genemodel is shown here, with less support see *Homolog*. *Cu_function*; indicates what the putative function of the protein is, and *Cu relevance* if this function is highly relevant for tolerating copper toxicity. The *INDEL corrected?* column indicates if, and how many, manual curations of the S. marinoi reference genomes reading frame, where made across the model. Deep coverage RNAseq ILLUMINA data and NanoPore full length mRNA sequencing data was used to support such corrections, which generally fixed premature stopcodons, and presumably corrected the ORF.  

**DRCpredictions.txt**
This file is derived from the StrainPhenotyping directory analysis, and contains all strain-specific phenotypic traits that has been quantified from mono-clonal experiments (e.g. DRC, microscopy, growth curves, etc). *Strain*; the strain name, *Predict*/*Predict_low*/*Predict_high*; the predicted growth rate at 8.65 uM Copper concentration with 95% conf. interval. The following columns are repeted for *EC05* (5% inhibition of growth) *EC50* (50% inhibition of growth), *EC95* (95% inhibition of growth), as ECxx, *ECxx_SE*, low, and high, which is copper concentration in uM, Standard error, 95% low and high range. *Comments* are other observations from DRC or monoclonal experiments. *Population*; (GP or VG), *Plate.ID*; an identifier of the plate where DRC was run, *Running_ID* is strains alias ID in the selection experiment. *Incubation Position* in the growth chamber, which could cause confunding effects due to diffrences in supplied *Light intensity* (umol PAR photons m-2, s-1). *Wait time (h)* between mixing copper media for DRC and adding of cells. *RO5_Growth_Rate*; the growth rate of strain RO5AC on the specific 24-well plate while *Growth rate (day-1)* refers to the mean growth rate of the 5 control treatments of the strain with *S.D.* standard deviation. *Ph* and maximum quantum yield of Photosystem II (*Fv/Fm*) of the inocoulum culure, *Start Density (RFU)* is the Relative Chl a flourescence value (mean across all strain 23 reps) at the start of the treatment (aim is 0.003, or ca. 1000 cells per mL/well), *Pre_Density*; RFU density of the pre-culture used to start DRC measurments and Selection experiment, *drc_Normalization*; Indicates if DRC data was normalized against the mean of the 5 control reps (Normal), or when a t-test showed suppression of growth in the controls, the 3 lowest DRC concentrations (Low), *Start_density*; the microscopically determined strain-specific cellular density going into the Selection experiment (cells mL-1),  *Relative_Start_Density*; same as *Start_density* but normalized to relative abundance, *Biovolume*, the avarage (N=20 observations) biovolume of a cell (um3),  *Surfacearea*, the avarage surface area of a cell (um2), *Surface_Vol_ratio*, surface to volume ratio of cell (um2 um-3)  *Diameter*, avarage diameter of cell (um),  *Diamter_Length_ratio*; avarage diameter to length ratio (um um-1), *CuGrowth*; DRC predicted growth rate under 8.65 uM Cu,  *CuGrowth_L*; low, and *CuGrowth_H*; high, 95% conf. intervall.


**Fitness.txt** [**Fitness_DL.txt** is version where strains without barcode observations have been set to detection limit rather than omitted]
This file contains the phenotypic growth rate traits derived from metabarcoding of changes in strains abundances during the first 9 days of the Selection experiment (created in the **SelectionExperiment** analysis). The first 3 column indexes strains. *Barcode_N* specifies how many times the strain was observed in the metabarcoding data (Nmax=5 except GPxControl where Nmax=1 ). *Barcode_mean*, (*Barcode_High*, *Barcode_Low*); metabarcoded measurments of growth rates during co-cultivation for 9 days (95% conf. intervals). *Growth_Rate* (High*, *Low*); monoclonal predictions of same growth rates during 3 days based on DRC (95% conf. intervals).


**StrainMeta.txt**
This file has metadata about strains mainly regarding how they where processed. Its used by the script mainly to index WGS sample names with strains. It has metadata on the strains (col *Strain*), WGS identifies (col *Seq_ID*) that is linking back to labels of DNA sample (col *DNA_ID*), and what there alias was in the experiments (col *Running_ID*). Also includer are *Site*, sediment *Core*, isotopic dating of sediment *Approx.Age*, what *Species* it is, the *Date_isolated* it was isolated, and who was the *Isolater*, what *Growth_Rate* (day-1), *Biovolume* um3, did it have when the selection experiment was seeded. *Comments* referes to qualitative observations DRC, growth or DNA extraction observations.

## Structure of script section by section

**GenomicAnalysis.R** [*PhenotypeAnalysis.R* is a trimmed version of this that only looks at phenotypic data, no genomically derived ones]

#Houskeeping the loading packages####
load packages and set working directory to place with 3 subdirectories named Input, Plots, and Results. Put Input data in Input, with **CuGeneCoverage/P14203_101_CuCoverage.tsv** in sub-directory named CuGeneCoverage.

#Read in data and look at it####
This part reads in the data for the three diffrent analysis, formats it, merges, and makes some QC plot of it. During the gene coverage data analysis, there is an option to look at homologus genes sperately, summed across the genome, or summed across contiges. Since gene duplications, which is what the analysis is looking for, occures lokally, the script currently merges homologs located on the same contig. In most of these cases the genes are in tandem repeats in the S. marinoi reference genome (see **CuMetabolism_SM_short.txt**). This can be altered on line 194. The output file *WGScoverageMean.txt* shows coverage stats for each gene model, across all strains, with AvarageCoverage_Summed.pdf as a graphical output.    

#Changes in coverage filter####
Here the data gets filtered for genes that actually have support for tandem repeats, or deletions in strains. 

#Genotype-Phenotype correlation#####
This first plots strain variation in gene copy number (coverage) between strains against EC50 and barcoded growth rates under copper stress. We expect a positive correlation if more copies/protein product provide tolerance (e.g. chelators and efflux transporters), or negative if it makes strains more suseptable to toxicity (e.g. influx transporter). After this a more general clustering and data exploration is performed via PCA, heatmaps, and multiple regression analysis. This analysis is sensitive to missing data, which needs to be removed on line 454 to 476. The way the metabarcoding data is handled in *SelectionExperiment*, producing the **Fitness.txt** file, analysis will affect this. For this exploratory analysis i prefere to set No observation to the detection limit, rather than removing such observations which creats NAs and lots of removed strains/reduced power. For such approch change the input file to **Fitness_DL.txt* on line 40. 