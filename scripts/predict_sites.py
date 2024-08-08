#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
This script reads a table of Ensembl transcript IDs searches the
3'-UTRs of each transcript for canonical miR-34a seed-binding sites.

The seed matches are defined as follows:
    - 8mer: Match to positions 2-8 of the miRNA followed by A.
    - 7mer-m8: Match to positions 2-8 of the miRNA.
    - 7mer-A1: Match to positions 2-7 of the miRNA followed by A.
    - 6mer: Match to positions 2-7 of the miRNA.
    - 6mer-m8: Match to positions 3-8 of the miRNA.
    - 6mer-a1: Match to positions 2-6 of the miRNA followed by an A.
'''


import argparse
import re

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


def find_seed_sites(gene, mir_seq, mrna_dict, site_types=None):
    '''
    Search an mRNA for canonical seed-binding sites for a miRNA.
    '''
    mrna_seq = mrna_dict[gene['Transcript ID']]
    seed = mir_seq[1:8].translate(str.maketrans('ACGU', 'UGCA'))[::-1]
    match_dict = {'8mer': seed+'A',
                  '7mer-m8': '|'.join([seed+'[^A]', seed+'$']),
                  '7mer-a1': '|'.join(['[^'+seed[0]+']'+seed[1:]+'A',
                                       '^'+seed[1:]+'A']),
                  '6mer': '|'.join(['[^'+seed[0]+']'+seed[1:]+'[^A]',
                                    '[^'+seed[0]+']'+seed[1:]+'$',
                                    '^'+seed[1:]+'[^A]']),
                  '6mer-m8': '|'.join([seed[:-1]+'[^'+seed[-1]+'][^A]',
                                       seed[:-1]+'[^'+seed[-1]+']$']),
                  '6mer-a1': '|'.join(['[^'+seed[0]+'][^'+seed[1]+']'+
                                       seed[2:]+'A',
                                       '^[^'+seed[1]+']'+seed[2:]+'A'])}

    sites = []
    if not site_types:
        site_types = list(match_dict.keys())
    for site_type in site_types:
        sites += [(i.start(1), i.end(1)-1, site_type) for i in re.finditer(
            r'(?=('+match_dict[site_type]+'))', mrna_seq)]

    records = dict(list(zip(*[['Seed start', 'Seed end', 'Site type'],
                              list(zip(*sites))])))
    return pd.Series(records)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--infile', type=str, default='../data/dge/dge.txt',
        help='DGE data file with transcripts IDs')
    parser.add_argument(
        '--outfile', type=str, default='../results/site_predictions/sites.txt',
        help='Output file with predicted miRNA target sites.')
    parser.add_argument(
        '--mir_seq_file', type=str, default='../data/seqs/mir34a.fa',
        help='Fasta file containing the miRNA sequence (only the first '\
             'sequence in the file is read).')
    parser.add_argument(
        '--mrna_seq_file', type=str, default='../data/seqs/utr_seqs.fa',
        help='Fasta file containing the mRNA sequences.')
    parser.add_argument(
        '--site_types', type=str, default='8mer|7mer-m8|7mer-a1|6mer',
        help='Site types to include in the prediction, separated by "|"')
    args = parser.parse_args()

    data = pd.read_csv(args.infile, sep='\t')
    data = data[['Gene ID', 'Transcript ID', 'Gene name']]

    mir_seq = [str(mir.seq) for mir in
               list(SeqIO.parse(args.mir_seq_file, 'fasta'))][0]
    mrna_dict = {mrna.id: str(mrna.seq) for mrna in
                 list(SeqIO.parse(args.mrna_seq_file, 'fasta'))}
    site_types = args.site_types.split('|')

    tqdm.pandas()
    data[['Seed start', 'Seed end', 'Site type']] = data.progress_apply(
        find_seed_sites, args=(mir_seq, mrna_dict, site_types), axis=1)
    data.dropna(subset=['Seed start', 'Seed end', 'Site type'], inplace=True)
    data = data.explode(column=['Seed start', 'Seed end', 'Site type'])

    data.to_csv(args.outfile, sep='\t', index=False, na_rep='NaN')


if __name__ == '__main__':
    main()
