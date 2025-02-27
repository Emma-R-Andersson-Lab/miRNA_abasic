# miRNA_abasic

This repository contains the scripts used to analyze the RNA-seq data in Kosek et al. "Mapping effective microRNA pairing beyond the seed using abasic modifications".

As input data, copy the "mRNA DGE" and "Protein DGE" sheets in Supplementary File 2 and save them as a tab-separated text files called "rna.txt" and "prot.txt" respectively in the ./data/dge/ folder. 3UTR sequences were downloaded from TargetScan (https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi) and Ensembl (https://www.ensembl.org/biomart/martview). Sequences used for the analysis can be found in the ./data/seqs/ folder. To generate the graphs in the paper, run runall.sh in the ./scripts/ folder followed by the Jupyter notebook for each figure.

There is a yml file listing the packages in the environment used to run the scripts and notebooks. To run the RNA structure predictions you need to install ViennaRNA (version 2.4.17) (https://www.tbi.univie.ac.at/RNA).
