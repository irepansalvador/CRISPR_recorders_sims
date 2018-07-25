# CRISPR_recorders_sims

This repository will contain the scripts from the work:

"Is it possible to reconstruct an accurate cell lineage using CRISPR recorders?
Authors: Irepan Salvador-Martínez, Marco Grillo, Michalis Averof
and Maximilian J Telford
bioRxiv doi: https://doi.org/10.1101/373357

There is a short description for each script at the beginning of the file 
stating also to which figure in the manuscript is associated with.

The scripts are contained in the following folders:

### MATLAB_sims/

This folder contains the MATLAB scripts and associated files necessary to
generate all the simulations described in the manuscript (Figures 2-3 and 5-7).
The scripts were run in MATLAB version 2017a (Linux).

### seq_reads_analysis/

This folder contains the Jupyter notebooks (Python and R) that were used to
analyse the (already trimmed) sequencing reads to get parameters to be used
for the simulations (mutational outcomes, mutation rates, etc).

### make_NEXUS/

This folder contains a Perl script used for generating a NEXUS file from the
MATLAB simulations' output.
For each simulation of 65K cells, the script makes 10 random samples of 1,000
cells. 

The samples of 1000 cells are converted to a NEXUS file used as input for the 
PAUP* software, to infer the cell lineage based on the accumulated mutations.

More specifically, for each sample, a "root" taxa is added (unmutated recorder)
to the simulated targets as the necessary PAUP blocks to be read by PAUP.
In the folder the "reference tree" used for estimating the accuracy (by 
comparing it to the "inferred tree" with the CompareTree.pl software) is also
provided ("16div_32_targets_60_states_REF.nw)

### GESTALT/

This folder contains the scripts used for characterising the mutational outcome
of the GESTALT v7 recorder, obtaining the parameters to perform simulations and 
quantify the accuracy of the GESTALT project (Figure 3 and Suppl. Figure 3)

The GESTALT data was obtained from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81713

### R_accuracy_heatmap

This folder contains the data and R script for re-creating the accuracy heatmap
shown in Figure 7C.
Also in the folder are the two tiff output images from the R script, a heatmap 
and a levelplot.
