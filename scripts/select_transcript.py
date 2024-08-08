#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
This script assigns an Ensembl Transcript ID to each Gene ID, for
purposes of miRNA target site prediction. The "representative
transcript" from the TargetScan database is preferentially chosen.
For genes with no transcript in this database, the Ensembl
"canonical" transcript is chosen instead.
'''


import numpy as np
import pandas as pd


def main():
    # Read gene info from Ensembl
    ens_gene = pd.read_csv(
        '../data/Ensembl/gene_info.txt', sep='\t', header=None, skiprows=1,
        names=['Gene ID', 'Transcript ID', 'Gene name',
               'Gene type', 'Canonical'])

    # Read gene info from TargetScan
    ts_gene = pd.read_csv('../data/TargetScan/gene_info.txt', sep='\t')
    # Keep the representative transcripts for each human gene
    ts_gene = ts_gene[(ts_gene['Species ID'] == 9606) &
                      (ts_gene['Representative transcript?'] == 1)]
    # Remove version numbers from IDs
    ts_gene['Gene ID'] = ts_gene['Gene ID'].apply(lambda x: x.split('.')[0])
    ts_gene['Transcript ID'] = ts_gene['Transcript ID'].apply(
        lambda x: x.split('.')[0])
    # Add an identifier for transcripts from the TargetScan database
    ts_gene['TargetScan'] = 1
    ts_gene = ts_gene[['Gene ID', 'Transcript ID', 'TargetScan']]

    # Merge Ensembl and TargetScan records
    transcripts = pd.merge(
        ens_gene, ts_gene, on=['Gene ID', 'Transcript ID'], how='left')
    transcripts_ts = transcripts.dropna(subset='TargetScan')
    transcripts_ens = transcripts[
        ~transcripts['Gene ID'].isin(transcripts_ts['Gene ID'])].copy()
    # Keep canonical Ensembl transcripts
    transcripts_ens.dropna(subset='Canonical', inplace=True)
    # Add an identifier for transcripts from the Ensembl database
    transcripts_ens['Ensembl'] = 1
    # Concatenate TargetScan/Ensembl lists
    transcripts = pd.concat([transcripts_ts, transcripts_ens])
    transcripts['Source'] = transcripts.apply(
        lambda x: 'TargetScan' if x['TargetScan'] == 1
        else 'Ensembl' if x['Ensembl'] == 1 else np.nan, axis=1)

    # Save transcript records
    transcripts.sort_values(by='Gene name', inplace=True)
    transcripts.drop(
        columns=['Canonical', 'TargetScan', 'Ensembl'], inplace=True)
    transcripts.to_csv('../data/dge/transcripts.txt', sep='\t', index=False)


if __name__ == '__main__':
    main()
