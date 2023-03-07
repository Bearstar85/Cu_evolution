2023-03-02
Author: Björn Andersson

# Environmental and experimental evolution of copper tolerance in the coastal diatom *Skeletonema marinoi*

## Experimental overview

This project explores if, and how, contemporary populations the costal diatom *Skeletonema marinoi* have evolved in respons to mining pollution. The study system is two semi-enclosed inlet in the Baltic Sea, where one Gåsfjärden (VG) has been affected by mining pollution for ca. 400 years, while the other, Gropviken (GP), has not. Strains where isolated, and the genome sequenced for 55 individual strains, and they where phenotyped in terms of specific growth rate and dose-responses to toxic copper concentrations (6-12 uM Cu). An artifical evolution experiment was conducted by assembling 28 and 30 strains from the two lacations sepperately, and let them evolve with, and without, toxic Cu stress. 8.65 uM Copper, corresponding to the concentration that inhibits the reference *S. marinoi* strain RO5AC’s specific growth rate with 50% in acute toxic tests (Andersson et al. 2020: DOI: 10.1016/j.aquatox.2020.105551.). Strain selection models are computed according to Andersson et al. 2022: DOI: 10.1038/s41396-021-01092-9). Raw sequencing data is not processed or included in the scripts and data here, but deposited and available at NCBI under BioProject PRJNA939970. The amplicon sequencing data is analysed as outlined in https://github.com/topel-research-group/Bamboozle/wiki/Bamboozle-Part-2:-Barcode-Quantification, without Dada2 denoising and instead counting only perfect sequence matches against library of known strain allele.

The file **SessionInfo.txt** shows what R-packages and versions the scripts could be executed with. 

![Experimetal Design](https://github.com/Bearstar85/Cu_evolution/blob/master/ExperimentalDesign1.jpg)

## Directory structure

**Barcodes**
This is an analysis pipeline of amplicon sequences of a hypervariable lucus in *Skeletonema marinoi* (SM). The lucus was bioinformaticaly identified based on analysis of whole genome sequences of 55 strains of SM from two Baltic Sea locations. It was predicted to have at least one uniqe allels enabling tracking of evolution through selection on standing genetic diversity in a artificiell evolution experiment (See Fig above). Two barcode loci (C2W24 and C12W1) were sequenced in the experiment, but C12W1 had much more allelic diversity so the majority of the analysis focus on this data (see Barcodes/Barcoding_C12W1/README.md for more information). This Git repository does not contain the bioinformatic sequence analyses, but starts after raw reads have been trimmed, merged, and mapped back to the known allele sequences. Two pre-process approces are included, one based on Dada2 error-correction (Barcoding_C12W1_BBmergerDADA2_v3), and one that uses exact matches of merged amplicon sequences, with (Barcoding_C12W1). The latter is the one we used for the publication: Andersson et al. Strain-specific metabarcoding reveals rapid evolution of copper tolerance in populations of the coastal diatom Skeletonema marinoi, in prep. For the upstream bioinformatic pipelines, see:

https://github.com/topel-research-group/Bamboozle/blob/master/scripts/BarcodeQuantification.zip, with instructions at https://github.com/topel-research-group/Bamboozle/wiki/Bamboozle-Part-2:-Barcode-Quantification

Output from this analysis goes into **SelectionExperiment**

**StrainPhenotyping**
This directory contains data and code to analyse mono-clonal experiments on singel strains. This includes maximum growth rate, dose-respons curve data, cellular dimensions, and other observations 

Output from this analysis goes into **SelectionExperiment**

**SelectionExperiment**
In this directory, output data from **StrainPhenotyping** and **Barcodes** are combined with observational data from the selection experiment, and analysed.

Output from this analysis goes into **Genomics**

**Genomics**
Here the phenotypic strain observations from **SelectionExperiment** is combined with genomic data (gene copy number variance from WGS analysis of 55 strains), to look for genotype-phenotype correlations. Tries to cluster strains and explain if copy-number variance in Cu metabolism genes can explain strain differences in Copper tolerance (it can't). 
