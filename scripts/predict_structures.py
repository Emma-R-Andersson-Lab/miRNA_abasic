#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Predict structures of target-miR complexes with RNAsubopt.
'''


import argparse
from subprocess import run, PIPE

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


def predict_ensemble(site, mir_seq, mrna_dict, tar_len, outdir):
    '''
    Predicts all structures of a target-miR complex in the energy
    range between the optimal structure and the seed-only structure,
    and writes them to a file.
    '''
    # Get target sequence
    tar_seq = _get_target_seq(site, mrna_dict, tar_len)

    # Get constraint
    tar_seq, constraint = _get_constraint(site, mir_seq, tar_seq)

    # Predict the MFE of the optimal structure
    opt_struct = _run_rnasubopt(tar_seq+'&'+mir_seq, 0, constraint)
    mfe_opt = opt_struct[0][1]

    # Predict the MFE of the seed-only structure
    seed_only_constraint = constraint.replace('.', 'x')
    seed_struct = _run_rnasubopt(tar_seq+'&'+mir_seq, 0, seed_only_constraint)
    mfe_seed = seed_struct[0][1]

    # Predict the full ensemble
    delta_mfe = abs(mfe_opt)-abs(mfe_seed)
    ensemble = _run_rnasubopt(tar_seq+'&'+mir_seq, delta_mfe, constraint)
    ensemble = _filter_structures(ensemble)

    # Save to file
    name = '_'.join([site['Transcript ID'], str(site['Seed start'])])
    with open(outdir+name+'.sub', 'w', encoding='utf-8') as file:
        file.write('>'+name+'\n')
        file.write(
            '\t'.join([tar_seq+'&'+mir_seq, str(mfe_opt), str(mfe_seed)])+'\n')
        for struct, mfe in ensemble:
            file.write('\t'.join([struct, str(mfe)])+'\n')


def _get_target_seq(site, mrna_dict, tar_len):
    '''
    Get the sequence of the predicted miRNA target site.
    '''
    mrna_seq = mrna_dict[site['Transcript ID']]
    tar_seq = mrna_seq[
        max([site['Seed end']-tar_len+1, 0]):site['Seed end']+1]
    return tar_seq


def _get_constraint(site, mir_seq, tar_seq):
    '''
    Generate base-pairing constraint to force pairing to the seed.
    '''
    # Enforce pairing to the seed.
    seed_len = 8 if site['Site type'] in ('8mer', '7mer-m8') else 7
    constraint = ('.'*(len(tar_seq)-seed_len)+'('*(seed_len-1)+'x&x'+
                  ')'*(seed_len-1)+'.'*(len(mir_seq)-seed_len))
    # Adjust target sequence and constraint for 7mer-m8 and 6mer sites
    # where the site end is the last base of the UTR (i.e. there is no
    # target base opposite g1).
    if (site['Seed end']-site['Seed start'] != 7
            and site['Site type'] in ('7mer-m8', '6mer')):
        tar_seq = tar_seq[1:]
        constraint = constraint.replace('x&x', '&x')
    return tar_seq, constraint


def _run_rnasubopt(seq, delta, constraint):
    '''
    Predicts all structures of an RNA sequences or a complex of two
    RNA sequences in a given energy range with RNAsubopt.
    '''
    args = ['RNAsubopt', '-e'+str(delta), '--sorted', '--en-only', '--noLP',
            '--constraint', '--enforceConstraint']
    input_seq = seq+'\n'+constraint

    rnasubopt = run(args, input=input_seq, stdout=PIPE, stderr=PIPE,
                    text=True, check=True)
    ensemble = []
    for record in rnasubopt.stdout.splitlines()[1:]:
        record = record.split()
        struct, mfe = record[0], float(record[1])
        ensemble.append((struct, mfe))
    return ensemble


def _filter_structures(ensemble):
    '''
    Removes structures with undesirable features from the ensemble.
    '''
    filtered_ensemble = []
    for struct, mfe in ensemble:
        # Remove intramolecular pairs
        if not struct.rfind('(') < struct.find('&') < struct.find(')'):
            continue
        filtered_ensemble.append((struct, mfe))
    return filtered_ensemble


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--infile', type=str, default='../results/site_predictions/sites.txt',
        help='Predicted miRNA target sites.')
    parser.add_argument(
        '--outdir', type=str, default='../results/RNAsubopt/',
        help='Directory for output files.')
    parser.add_argument(
        '--mir_seq_file', type=str, default='../data/seqs/mir34a.fa',
        help='Fasta file containing the miRNA sequence (only the first '\
             'sequence in the file is read).')
    parser.add_argument(
        '--mrna_seq_file', type=str, default='../data/seqs/utr_seqs.fa',
        help='Fasta file containing the mRNA sequences.')
    parser.add_argument(
        '--tar_len', type=int, default=32,
        help='Number of bases in the target site, including the seed.')
    args = parser.parse_args()

    sites = pd.read_csv(args.infile, sep='\t')

    mir_seq = [str(mir.seq) for mir in
               list(SeqIO.parse(args.mir_seq_file, 'fasta'))][0]
    mrna_dict = {mrna.id: str(mrna.seq) for mrna in
                 list(SeqIO.parse(args.mrna_seq_file, 'fasta'))}

    tqdm.pandas()
    sites.progress_apply(
        predict_ensemble,
        args=(mir_seq, mrna_dict, args.tar_len, args.outdir), axis=1)


if __name__ == '__main__':
    main()
