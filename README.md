2022-05-05
Author: Björn Andersson

# Environmental and experimental evolution of copper tolerance in the coastal diatom Skeletonema marinoi

## Experimental overview

This project explores if, and how, contemporary populations the costal diatom *Skeletonema marinoi* have evolved in respons to mining pollution. The study system is two semi-enclosed inlet in the Baltic Sea, where one Gåsfjärden (VG) has been affected by mining pollution for ca. 400 years, while the other, Gropviken (GP), has not. Strains where isolated, and the genome sequenced for 55 individual strains, and they where phenotyped in terms of specific growth rate and dose-responses to toxic copper concentrations (6-12? uM Cu). An artifical evolution experiment was conducted by assembling 28 and 30 strains from the two lacations sepperately, and let them evolve with, and without, toxic Cu stress. 7.92 uM Cu, corresponding to the concentration that inhibits our reference *S. marinoi* strain RO5AC’s specific growth rate with 50% in acute toxic tests (Andersson et al. 2020). The Git repository is currently being populated with data and scripts, with only the Barcode analysis in a complete first version state.   

## Directory structure

**Barcodes**
This is an analysis pipeline of amplicon sequences of a hypervariable lucus in *Skeletonema marinoi* (SM). The lucus was bioinformaticaly identified based on analysis of whole genome sequences of 55 strains of SM from two Baltic Sea locations. It was predicted to have at least one uniqe allels enabling tracking of evolution through selection on standing genetic diversity in a artificiell evolution experiment (See Fig below). Two barcode loci were sequenced in the experiment, but C12W1 had much more allelic diversity so the majority of the analysis focus on this data (see Barcodes/Barcoding_C12W1/README.md for more information).

![Experimetal Design](https://github.com/Bearstar85/Cu_evolution/blob/master/ExperimentalDesign1.jpg)
