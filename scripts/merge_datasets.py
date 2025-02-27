#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
This script reads the RNA-seq and mass spec differential gene expression
files and merges the records into a single file.
'''


import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--rna_dge_file', type=str, default='../data/dge/rna.txt',
        help='Processed RNA-seq DGE data file.')
    parser.add_argument(
        '--prot_dge_file', type=str, default='../data/dge/prot.txt',
        help='Processed mass spec DGE data file.')
    parser.add_argument(
        '--transcript_file', type=str, default='../data/dge/transcripts.txt',
        help='File mapping Gene IDs to Transcript IDs.')
    parser.add_argument(
        '--outfile', type=str, default='../data/dge/merged.txt',
        help='Processed DGE data file.')
    args = parser.parse_args()

    columns = ['Gene ID', 'Transcript ID', 'Gene name', 'UniProt IDs']

    # Read RNA-seq dataset
    rna_data = pd.read_csv(args.rna_dge_file, sep='\t')
    rna_data.columns = [i+' mRNA' if not i in columns
                        else i for i in rna_data.columns]

    # Read mass spec dataset
    prot_data = pd.read_csv(args.prot_dge_file, sep='\t')
    prot_data.columns = [i+' prot' if not i in columns
                         else i for i in prot_data.columns]

    # Merge datasets
    data = pd.merge(
        rna_data, prot_data, on=['Gene ID', 'Gene name'], how='outer')

    # Get transcripts
    transcripts = pd.read_csv(args.transcript_file, sep='\t')
    data = pd.merge(data, transcripts[['Gene ID', 'Transcript ID']])
    assert data['Gene ID'].nunique() == data['Transcript ID'].nunique()

    # Save merged dataset
    data.sort_values(by='Gene name', inplace=True)
    data = data[columns+[i for i in data.columns if i not in columns]]
    data.to_csv(args.outfile, sep='\t', index=False, na_rep='NaN')


if __name__ == '__main__':
    main()
