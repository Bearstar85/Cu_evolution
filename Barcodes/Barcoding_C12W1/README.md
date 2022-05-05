2022-04-27
Author: Björn Andersson

# Analysis of intra-specific metabarcoding using amplicon sequencing of C12W1 lucus in Skeletonema marinoi

## Script and experimental summary

This scripts **Barcode_analysis.R** is a custom made analysis pipeline of amplicon sequences of a hypervariable lucus in *Skeletonema marinoi* (SM). The lucus was bioinformaticaly identified based on analysis of whole genome sequences of 55 strains of SM from two Baltic Sea locations (Gåsfjärden: VG, and Gropviken: GP). It was predicted to have at least one uniqe allels for every strain used in the experiment.  In summary the script requires sample data file(s) that have already been mapped against a refrence database, i.e. numerical counts of observations of specific sequences. In our case this corresponds to counts of allels of 59 diffrente strains of SM (e.g. **C12W1_abundances/Input/P21502_101.tsv**). The **Indexing.txt** file contains metadata on each .tsv datafiles sample which comes from both control samples of sequences from individual genotypes, as well as mixed samples during a selection experiment. This experiment is described in detail in (REF) but briefely strains from the two populations were mixed at equal cell densities, grown semi-continously in the exponental growth phase for 42 days (50-100 generations) with toxic copper stress (Cu) or without (C), and DNA samples where collected on day 0 (MaterMix, or MM), 9 (t9), or 42 (t42). The data files are then processed according to allelic information in **Allele_indexing.txt** file. This file is complex and not all information is essential to the current script function (See section *Input files and Metadata* below for more information). 

In short the script reads the input files (all files with .tsv extention in the Input folder), extracts the sample information an matches it with the correct metadata. It then parse out non-uniqe allel observations using differential equations and the abundance of the secondary, uniqe, allel of strains with shared allels. It then normalizes the allel observations to relativt abundance using the sum of all amplicons matching known allels (consequently, excluding all non-perfect matches). While doing this it generates QC figures for each sample into the **C12W1_abundances/Plots directory**, and outputs the normalized data into individual files in **C12W1_abundances/Indexed** directory. It also output the total amplicon count into **C12W1_abundances/Read_counts**. The rest of the script explores the data by generating a number of graphs that are exported into the **C12W1_abundances/Plots directory**. As a warning the script seemingly randomly crashes on my machine in numeras of these loops, which is why Quartz is sometime called on and generating pop-up windows.

The non-perfect allel matching amplicons (unknowns) are also quantified and analysed for patterns. This is through the **unknowns.R** script uses a simillar directory and data structure as **C12W1_abundances** described above, but separerade into **C12W1_unknowns** directory and with observations of each uniqe unknown amplicon sequence, rather than matches to known allels.

![Experimetal Design1](https://github.com/Bearstar85/Cu_evolution/blob/master/ExperimentalDesign1.jpg)

## Input files and metadata

**C12W1_abundances/Input/P21502_101.tsv** etc..
The "Barcode" column contains the "Barcode identifier, which is what is used to join with metadata from **Indexing.txt**
**The Indexing.txt**. "Primer 1" and "Primer 2" has the absolute number of amplicon observations of each of two possible primer versions (primers are degenerative to account for a singel SNP site), and column "Total" is the sum of these observations. Importantly this file excludes all amplicons with imperfect matches to the satabase, and it counts non-uniqe allels multiple times. So the sum of the "Total" column is not equivalent to the sum of amplicons per sample. But the script deals with this during processing.

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

**C12W1_abundances/C12W1_unknowns/P21502_101.tsv** etc.
Here are observations of amplicans that are not know to actually exist in any strains and are presumably mostly made up of artefacts from sequencing errors or the PCR reaction. 
"Quantity" how many times this amplicon was seen in the given sample
"Unknown barcode" the sequence of the specific amplicon, where Y is used to mask the degenerative primer site.

**C12W1_all.withR05.v2.fasta** 
Contains the database of known allel sequences

*Directory structure*
The script needs a directory structure where the R-scripts and 

## Structure of scripts section by section

#Barcode_analysis.R, the main script

#Houskeeping the loading packages####
Line 1-41 does houskeping, sets work directory, gives som Bash commands for prepping files with wrong name, loads packages etc.


#Data normalization####
line 42-53 reads in the data. A WARNING, certain Indexing names in Allels (column "Diffrential_ID") will not comply with differential equation function (notably x and df cannot be used).
Line 54-70 code that setsup the differential equation command.
Line 82-196 This for-i-loop goes through each .tsv file in the Input directory and append indexing, solve the differential equation, normalizes the absolute amplicon observations to the relative abundance, outputs the indexed data to Index directory, (note that csv format is used with .txt extention to file name!), writes out a header (this file is needed before the loop is run the first time!), append the data to the dataframe "All_reads", and makes QC plots for each file and saves them as .pdfs in /Plots sub-director.

#Summarize the results####
This whole section of code should not be needed anymore.

#Graphical interpretation of data####
From here on the data only gets subselected, restructured, and ploted.
Line 230-335 focus on analys of the time zero mastermix samples, and false possitive observations in them. Most relevant results are outputed as figures into /Plots subdirectory.
Line 337-357 plots how the individual strains allels are behaving during the selection experiment. This loop is complex and crashes on my machine so there are multiple places where it can be cut short to avoid this by activating or inactivating end }.
However the first lines of code in the loop must be executed for downstream parts of the script to function (Specifically the dataframe AllelRatios needs to be filled on Line 362-385).

#Allele ratios and correlations####
This part of the script explores how the two allels within a strain covary.
Line 461-595 lots of exploration into how the allels covary and what false possitive signals we see in the data. Graphs are outputed into  Plots/Strains directory. My personal favorit graph is "FigReg.pdf".

#Clustering of samples####
This part of the script is unfinished in part since it was judged to be redundant.
What it is doing is restructuring tha data to a Matrix (barcode x sample) that can be used for heatmap clustering, PCAs and correlation analysis. Given how clearly the strain specific alleles covary, and how agreeable replication samples are, this type of analysis did not yield much further insight.

#**unknowns.R**, the artifact analysis script
This script analyses patterns of occurances in amplicons that do not match any of the known allels present in samples. These can represent a number of PCR artifacts (e.g. polymeras induced SNP, polymerase slipage, chimeric amplicons from two diffrent allele templates), misscalls from the sequencing, unknown mutations in cell-lines or sexual/asexual recombinant allels, contamination from allels of strains not sequenced. The script was originally created to also identify true allels of strains that where lacking reference whole genome sequences, or reference genotyp amplicon sequences. This was done in an iterative fashion and these putative allels where identified and added to the database of known allels, and removed as unknown amplicon sequences in this data. Note that the script also incoorperate some of the output files from **Barcode_analysis.R** but it is not essential for the function of the script that this works, only for a few selected QC plots of read counts per sample. It does require that the **Indexing.txt** file is in its proper location from the **Barcode_analysis.R**. Overall the structure of files and the script structure is very simillar to **Barcode_analysis.R** .

#Houskeeping the loading packages####
Line 1-37 does houskeping, sets work directory, gives som Bash commands for prepping files with wrong name, loads packages etc.

#Process data#####
Line 40-46 reads in data.
Line 47-92 for-i-loop that indexes and merge the input data into one tall dataframe (XXX), as well as counting the reads per sample, plotting it,  and storing output these results into Read_counts and Indexed, and Plots directory.

#Read_count_summary#####
Line 85-336, reads in read count files that needs to manual created by UNIX command cat *.csv > ReadCount.txt. and cat *.csv > Unknown_counts.txt  in the two Read_counts directory. It then contrast the number of reads that mapp to the database of allele, versus those that do not, and makes diffrent graphical interpretations of this. This part of the script can be skipped without affecting downstream code,or modified for other data set.

#Find true unknown sequences#####
This part of the script uses the MasterMix, time zero mixtures of all strains, to identify abundant unknown sequences that may be true allels (already found and removed in our data) or artefacts that are so abundant they approch concentrations of true allels in samples.
Line 338-493 reads in mastermix samples for the two populations, and merge and filtere so that only amplicon seen in all samples of a population are retained. There is an option to search for sequences here. Graphical interpretations are made of the abundances of true allels are contrasted against unknowns. There is an option to manualy set the filtering of unknown reads based on mean abundance of reads across all replicates (Line 479). Some results files are outputed into Results directory. 

#High abundant reads across all samples####
This script catalogs all unknown amplicons and explores patterns in them across all samples.
Line 494-513 This script joins all the unknown amplicons across all samples which takes some time to do (20 min on my machine). It starts from the amplicons that where identified in the Mastermixes, and gives each new amplicon a running numerical ID from there on. So the first 700 are mastermix samples. The dataframe is saved and outputed into Results/Summary_all_unknow_allels_RA.txt, or alternatively accessed file from that point to skip this step.
Line 524-579 this adds upp observations of unknown amplicons and does filtering for the most abundant once (user changable). Some steps take a few minutes to run, and then it offers filtering options. This are based on sum of realativa abundance observations across samples, which is a bit abstract. The original scriped used raw counts, which can bias towards high coverage sample, but the names of the filterstep has retains this lie 557-559 (e.g. High_RA_means100 was, and is still roughly equivelent to more than 100 amplicons observed across all sampels). The numbers in the comments are based on raw read observations, whereas the script currently uses relative abundances. This can be modified early on in script, when reading in the Input data.
Line 579-697 This focuses on one filter category (currently amplicons seen more than sum of 0.1 fractions abundance, which could mean that ot is present at 10% in one sample or 0.1% in 100 samples: High_RA_means500), INdexes samples and generates plots of where these amplicons are seen,. The for-i-loop is prone to crashing on my machine, see alternative stop points #} to mitigate this.
Line 699-758 uses the High_RA_means500 filtered data, makes a indexed matrix and finaly a heatmap that clusters samples both based on occurence patterns across all samples, and also clusters samples based on occurences of amplicons. This shows which amplicons are linked across samples and presumably are related to artifacts from the genotype sample it is formed from. At least that is my interpretation of the data. Results matrixes are outputed into Results directory.

#Adapth the plots to the abundant Mastermix samples instead####
Line760-923 repeates the analysis (from Line 557-758) but focus only on the samples seen in high abundance in the Mastermixes/t0 samples. In the first iteration this was 124 amplicons above 20 reads per sample, but we have now annoterade some as real allels, and also changed the QC filtering of data so this number changes, and can be manual changed where "unknowns2" are filtered out (I set to 10 now to get more amplicons included, 93 now) and also on Line 760 where "MM_all" is filtered out (filter need to adjusted to N rows in "unknowns2") . Note that the indexing number is used to select these strains, and that the entire script need to be modified to change this subselection.  Results matrixes are outputed into Results directory.






 

 


      

