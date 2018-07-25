# CRISPR_recorders_sims

This repository will contain the scripts from the work:

"Is it possible to reconstruct an accurate cell lineage using CRISPR recorders?
Authors: Irepan Salvador-Martínez, Marco Grillo, Michalis Averof
and Maximilian J Telford
bioRxiv doi: https://doi.org/10.1101/373357

There is a short description for each script at the beginning of the file 
stating also to which figure in the manuscript is associated with.

The scripts are contained in the following folders:

## MATLAB_sims/

This folder contains the MATLAB scripts and associated files necessary to
generate all the simulations described in the manuscript.
The scripts were run in MATLAB version 2017a (Linux).

## seq_reads_analysis/

This folder contains the Jupyter notebooks (Python and R) that were used to
analyse the (already trimmed) sequencing reads to get parameters to be used
for the simulations (mutational outcomes, mutation rates, etc).

## make_NEXUS/

This folder contains a Perl script used for generating a NEXUS file from the
MATLAB simulations' output.
This NEXUS file was then used as input for the PAUP* software, to infer the 
cell lineage based on the accumulated mutations.

More specifically, the scripts adds a "Root" (the unmutated CRISPR recorder)
to the simulated targets and adds the necessary PAUP blocks to be read by PAUP.
In the folder it is the "reference tree" used for estimating the accuracy (by 
comparing it to the "inferred tree" with the CompareTree.pl software) is also
provided.

