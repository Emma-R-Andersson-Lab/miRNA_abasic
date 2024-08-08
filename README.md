# miRNA_abasic

This repository contains the scripts used to analyze the RNA-seq data in Kosek et al. "Mapping effective microRNA pairing beyond the seed using abasic modifications".

As input data, copy the "RNA-seq data" sheet in Supplementary File 2 and save it as a tab-separated text file called "dge_raw.txt" in the ./data/dge/ folder. 3UTR sequences were downloaded from TargetScan (https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi) and Ensembl (https://www.ensembl.org/biomart/martview). Sequences used for the analysis can be found in the ./data/seqs/ folder. To generate the graphs in the paper, run runall.sh in the ./scripts/ folder followed by the Jupyter notebook for each figure.

There is a yml file listing the packages in the environment used to run the scripts and notebooks. Addtional requirements are ViennaRNA (version 2.4.17) (https://www.tbi.univie.ac.at/RNA) and the UpSetPlot library (https://upsetplot.readthedocs.io/) to draw Fig. 4C.
