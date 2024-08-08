#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
This script creates a fasta file containing the UTR sequence of each
Transcript ID assigned by the select_transcript.py script. Sequences
from the TargetScan database are used for transcripts listed as the
"representative transcript" for a gene in this database, otherwise
sequences from the Ensembl database are used.
'''


import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    # Read transcript list
    transcripts = pd.read_csv('../data/dge/transcripts.txt', sep='\t')

    # Read UTR sequences from TargetScan
    ts_utr_seqs = pd.read_csv(
        '../data/TargetScan/UTR_sequences.txt', sep='\t', header=None,
        skiprows=1, usecols=[0, 4], names=['Transcript ID', 'Sequence'])  
    ts_utr_seqs['Transcript ID'] = ts_utr_seqs['Transcript ID'].apply(
        lambda x: x.split('.')[0])
    ts_utr_seqs = ts_utr_seqs[ts_utr_seqs['Transcript ID'].isin(
        transcripts[transcripts['Source'] == 'TargetScan']['Transcript ID'])]

    # Read UTR sequences from Ensembl
    ens_utr_seqs = list(SeqIO.parse(
        '../data/Ensembl/UTR_sequences.fa', 'fasta'))
    ens_utr_seqs = pd.DataFrame(
        [[record.id, str(record.seq)] for record in ens_utr_seqs])
    ens_utr_seqs.columns = ['Transcript ID', 'Sequence']
    ens_utr_seqs = ens_utr_seqs[ens_utr_seqs['Transcript ID'].isin(
        transcripts[transcripts['Source'] == 'Ensembl']['Transcript ID'])]

    # Merge and save
    utr_seqs = pd.concat([ts_utr_seqs, ens_utr_seqs])
    records = []
    for transcript_id, seq in utr_seqs.values.tolist():
        seq = seq.upper().replace('T', 'U').replace('-', '')
        record = SeqRecord(
            Seq(seq), id=transcript_id, description=transcript_id)
        records.append(record)
    SeqIO.write(records, '../data/seqs/utr_seqs.fa', 'fasta')


if __name__ == '__main__':
    main()
