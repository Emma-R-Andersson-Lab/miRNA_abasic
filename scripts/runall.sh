#!/usr/bin/env bash
set -e
set -u
set -x


RAW_DGE_FILE=../data/dge/dge_raw.txt
TRANSCRIPT_FILE=../data/dge/transcripts.txt
DGE_FILE=../data/dge/dge.txt

MIR_SEQ_FILE=../data/seqs/mir34a.fa
MRNA_SEQ_FILE=../data/seqs/utr_seqs.fa

SITES_FILE=../results/site_predictions/sites.txt
STRUCT_DIR=../results/RNAsubopt/


# Remember to download external files from Ensembl and TargetScan prior
# to running these two scripts. The output files of these scripts
# are included in the depository.
python select_transcript.py
python make_fasta_files.py

# This script reads the raw differential gene expression file,
# assigns a transcript to each Gene ID, and filters out genes which
# are not protein-coding.
python process_rnaseq_data.py \
    --infile ${RAW_DGE_FILE} --outfile ${DGE_FILE} \
    --transcript_file ${TRANSCRIPT_FILE}

# This script predicts canonical seed-binding sites for miR-34a.
python predict_sites.py \
    --infile ${DGE_FILE} --outfile ${SITES_FILE} \
    --mir_seq_file ${MIR_SEQ_FILE} --mrna_seq_file ${MRNA_SEQ_FILE} \
    --site_types "8mer|7mer-m8|7mer-a1|6mer"


python predict_structures.py \
    --infile ${SITES_FILE} --outdir ${STRUCT_DIR}/ \
    --mir_seq_file ${MIR_SEQ_FILE} --mrna_seq_file ${MRNA_SEQ_FILE}


# After this, you can run the Jupyter notebooks to generate the figures.
