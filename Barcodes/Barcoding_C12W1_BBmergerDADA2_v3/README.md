2023-02-16
Author: Björn Andersson

# Analysis of intra-specific metabarcoding using amplicon sequencing of C12W1 lucus in Skeletonema marinoi

## Script and experimental summary

This scripts **Barcode_analysis.R** is a custom made analysis pipeline of amplicon sequences of a hypervariable locus in *Skeletonema marinoi* (SM). The locus was bioinformaticaly identified based on analysis of whole genome sequences of 55 strains of SM from two Baltic Sea locations (Gåsfjärden: VG, and Gropviken: GP). It was predicted to have at least one uniqe allels for every strain used in the experiment.  In summary the script requires sample data file(s) that have already been mapped against a refrence database, i.e. numerical counts of observations of specific sequences. In our case this corresponds to counts of allels of 59 diffrente strains of SM (e.g. **C12W1_abundances/Input/P21502_101.tsv**). The **Indexing.txt** file contains metadata on each .tsv datafiles sample which comes from both control samples of sequences from individual genotypes, as well as mixed samples during a selection experiment. This experiment is described in detail in (Andersson et al. 2023, in prep for Molecular Ecology) but briefely strains from the two populations were mixed at equal cell densities, grown semi-continously in the exponental growth phase for 42 days (50-100 generations) with toxic copper stress (Cu) or without (C), and DNA samples where collected on day 0 (MaterMix, or MM), 9 (t9), or 42 (t42). The data files are then processed according to allelic information in **Allele_indexing.txt** file. This file is complex and not all information is essential to the current script function (See section *Input files and Metadata* below for more information). 

In short the script reads the input files (all files with .tsv extention in the Input folder), extracts the sample information an matches it with the correct metadata. It then parse out non-uniqe allel observations using differential equations and the abundance of the secondary, uniqe, allel of strains with shared allels. It then normalizes the allel observations to relativt abundance using the sum of all amplicons matching known allels (consequently, excluding all non-perfect matches). While doing this it generates QC figures for each sample into the **C12W1_abundances/Plots directory**, and outputs the normalized data into individual files in **C12W1_abundances/Indexed** directory. It also output the total amplicon count into **C12W1_abundances/Read_counts**. The rest of the script explores the data by generating a number of graphs that are exported into the **C12W1_abundances/Plots directory**. As a warning the script seemingly randomly crashes on my machine in some of of these loops, which is why Quartz is sometime called on and generating pop-up windows. A work around is to disable the visualization of plots.

Note that all ASV sekvenser that does not by name match any of the list of alleles **Indexing.txt** file will be ignored in the analysis. This can be modified by adding them to the indexing files.

![Experimetal Design1](https://github.com/Bearstar85/Cu_evolution/blob/master/ExperimentalDesign1.jpg)

## Input files and metadata

**C12W1_abundances/Input/P21502_101.tsv** etc..
The "Barcode" column contains the "Barcode identifier, which is what is used to join with metadata from **Indexing.txt**, and column "Total" is observations of this sequence in the sample. Note that all ASV sekvenser that does not by name match any of the list of alleles **Indexing.txt** file will be ignored in the analysis. This can be modified by adding them to the indexing files.

**Indexing.txt**
This file has metadata on the .tsv input files. Column "ID" corresponds to the sample name and is used to join with "Barcode" cin datafiles. 
"Sample" has a uniqe identifier information about the sample, e.g. if its from a singel genotype: *GP2-4_26*, or from a mixed selection sample: *PIPT_GP_Cu5_t9_C12W1*. 
"Experiment" has category information on if its a Strain/Gentype.  
"Population" has category information on what populations samples is from (GP or VG)
"Treatment" has category information on treatment (Cu or C, NA for strains)
"Timepoint" has numerical information on time point in days (NA for strains)
"Bottel" has category information on uniqe identifier of each experimental bottel (e.g. VG_C1, NA for strains)
"Replicate" has integer information on replicate number per treatment-population combination (1-5, NA for strains)

**Allele_indexing.txt**
Note that currently the **Barcode_analysis.R** script uses only information in columns: "Barcode", "Diffrential_ID", "Strain" and "Differential_equation". 

"Barcode" contains uniqe identifier of a strains allel, with asterisks signifying when an allele is shared between 2 (one asteriks), or 3 strains (two asteriks), in the database. Used to join with .tsv file so this column must be identical to that for propper matching.
"Diffrential_ID" contains a short uniqe identifier for each strains allel (eg. ac, or aa1 and aa2 for a shared allel)
"Strain" contains information about the strain that contain the allel (e.g. strain VG1-2_94, or GP2-4_45and46 which are clones, two strains VG1-2_99or65 lack genotype sample but have been putitatively linked to one out of two individual).
"Differential_equation" contains the arithmetic formula that informs the script how observations of non-uniqe allels are parsed based on abundance of other alleles in the sample. (E.g ac for the uniqe allele ac, or aa1 * (an / (an + bt)) for the shared allel aa1, where an is the second allele of GP2-4_26 and bt the second one for GP2-4_71)

Other columns are non-essential but:
"Uniqe_all" Binary indexing if allele is uniqe amongst all allels in database
"Uniqe_pop" Binary indexing if allele is uniqe amongst all allels in individual population (GP or VG)
"Homologs" alphabetical indexing if which allels have homologs or are uniqe.
"Index" conserved row number for sorting.

**Sequences/C12W1_all.withR05.v2.fasta** 
Contains the database of known allel sequences

**Sequences/ASV_seqs.fasta** 
Contains the database of all ASVs sequences of this iteration of the downstream analysis

*Directory structure*
The script needs a directory structure where the R-scripts and 

## Structure of scripts section by section

# Only files in Input 
#Barcode_analysis.R, the main script

#Houskeeping the loading packages####
Line 1-41 does houskeping, sets work directory, gives som Bash commands for prepping files with wrong name, loads packages etc. Note that Input file names need to match the Indexing files. If not, execute the Bash code on line 11-13 in the Input directory. Also the working directory need to be changed on line 6 and 44, for the script to run proper.


#Data normalization####
line 42-53 reads in the data. A WARNING, certain Indexing names in Allels (column "Diffrential_ID") will not comply with differential equation function (notably x and df cannot be used).
Line 54-70 code that setsup the differential equation command.
Line 84-198 This for-i-loop goes through each .tsv file in the Input directory and append indexing, solve the differential equation, normalizes the absolute amplicon observations to the relative abundance, outputs the indexed data to Index directory, (note that csv format is used with .txt extention to file name!), writes out a header, append the data to the dataframe "All_reads", and makes QC plots for each file and saves them as .pdfs in /Plots sub-director. Note that we have quarantined cross-contaminated samples by moving them from /Input into the /Contaminated directory

#Graphical interpretation of data####
From here on the data only gets subselected, restructured, and ploted.
Line 202-300 focus on analys of the time zero mastermix samples, and any false possitive observations in them. Most relevant results are outputed as figures into /Plots subdirectory.
Line 337-357 plots how the individual strains allels are behaving during the selection experiment, and in the individual genotype samples. This loop is complex and crashes on my machine so there are multiple places where it can be cut short to avoid this by activating or inactivating end }.
However the first lines of code in the loop must be executed for downstream parts of the script to function (Specifically the dataframe AlleleRatios needs to be filled on Line 313-352).

#Allele ratios and correlations####
This part of the script explores how the two allels within a strain covary. This script works best if all strains have 2 alleles, but it will try and plot singel alleles and should work on haploid or homologus strain data. 
Line 433-575 lots of exploration into how the allels covary and what false possitive signals we see in the data. Graphs are outputed into  Plots/Strains directory. My personal favorit graph is "FigReg.pdf". The plots are slow to be built, but usually re-runing the print command will work if it exited with an error.

#Clustering of samples####
This part of the script was used to track down deviating alleles, and cluster occurens patterns of unknown ASVs during trobleshooting. It is somewhat unfinished since it became redundant once all alleles were identified.
What is done is a restructuring tha data to a Matrix (barcode x sample) that can be used for heatmap clustering, PCAs and correlation analysis. Given how clearly the strain specific alleles covary now, and how agreeable replication samples are, this type of analysis did not yield much further insight for us, But simillar approches may be usefull for troubleshooting new datasets in an iterative fashion in the known allels ASVs, together with unknown ASVs where linked with certain strains (chimeras or common artifactual ASVs that slipped through the pipeline).

#StrainCounts####
This is a modified way of summing up the strain counts by there alleles. It can silence alleles prone to artefacts, and compensate for this through observations of the second allele. The input and output is smilar as the Data normalization step above, but uses other parts of the indexing file, e.g. the equations in the "Differential_equation_strains" column.