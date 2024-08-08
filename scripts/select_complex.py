#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
Select the best structure from the RNAsubopt predictions.
'''


from dataclasses import dataclass
from itertools import chain, islice

import numpy as np
import pandas as pd


@dataclass
class Constraints:
    force_paired: tuple = ({},)
    min_pairs: int = 0
    offset: tuple = (-float('inf'), float('inf'))
    no_gu: bool = False
    no_gu_bases: tuple = ()
    max_mm_count: float = float('inf')
    max_mm_size: float = float('inf')


def get_top_structure(site, subopt_dir, constraints=None, return_bool=False):
    '''
    Returns the predicted MFE structure of a target-miR complex
    which meets a specified set of constraints.
    '''
    if not constraints:
        constraints = Constraints()

    seq, ensemble = _get_struct_ensemble(
        site['Transcript ID'], site['Seed start'], subopt_dir)
    seed_len = 8 if site['Site type'] in ('8mer', '7mer-m8') else 7

    top_struct = [np.nan]*6
    for struct, mfe in ensemble:
        pairs = _get_basepairs(struct)
        non_seed_pairs = [(i, j) for i, j in pairs if j > seed_len]

        if len(non_seed_pairs) < constraints.min_pairs:
            continue

        offset = _get_offset(non_seed_pairs)
        if not constraints.offset[0] <= offset <= constraints.offset[1]:
            continue

        if not _check_forced_pairs(non_seed_pairs, constraints):
            continue

        gu_pairs = _get_gu_pairs(seq, non_seed_pairs)
        if not _check_gu_pairs(gu_pairs, constraints):
            continue

        if not _check_mismatches(non_seed_pairs, constraints):
            continue

        site_start = site['Seed end']-struct.lstrip('.').find('&')+1
        site_start = (site['Seed start'] if site_start > site['Seed start']
                      else site_start)
        tar_pairs, mir_pairs = (list(map(list, zip(*non_seed_pairs)))
                                if non_seed_pairs else [[], []])
        top_struct = [site_start, tar_pairs, mir_pairs, mfe, offset, gu_pairs]
        break

    if return_bool:
        return not pd.isnull(top_struct[0])
    records = dict(
        zip(*[['3P start', 'Target pairs', 'miRNA pairs', 'MFE', 'Offset',
               'GU pairs'], top_struct]))
    return pd.Series(records)


def _get_struct_ensemble(transcript_id, seed_start, subopt_dir):
    '''
    Read predicted structures from RNAsubopt output files.
    '''
    ensemble = []
    filepath = ''.join(
        [subopt_dir, '/', transcript_id, '_', str(seed_start), '.sub'])
    with open(filepath, 'r', encoding='utf-8') as file:
        for line in islice(file, 1, None):
            line = line.rstrip('\n').split()[:2]
            ensemble.append((line[0], float(line[1])))
    seq = ensemble[0][0]
    ensemble = ensemble[1:]
    return seq, ensemble


def _get_basepairs(struct):
    '''
    Get the positions of paired bases in a target-miRNA complex.
    Returns a list of tuples with the format (target base, miR base).
    Target bases are numbered from the 3'-end (starting with 1),
    miRNA bases from the 5'-end (starting with 1).
    '''
    tar_len = struct.find('&')
    open_pairs, pairs = [], []
    for i, base in enumerate(struct.replace('&', '')):
        if base == '(':
            open_pairs.append(i)
        elif base == ')':
            j = open_pairs.pop()
            pairs.append((abs(j-tar_len), i-tar_len+1))
    return pairs


def _get_offset(non_seed_pairs):
    '''
    Calculate the offset between the seed and 3'-pairing helices.
    '''
    try:
        offset = non_seed_pairs[0][0]-non_seed_pairs[0][1]
    except IndexError:
        offset = 0  # No pairs beyond the seed.
    return offset


def _check_forced_pairs(non_seed_pairs, constraints):
    '''
    Check if miRNA bases which must be paired according to the
    constraints are paired in the structure. 
    '''
    if not constraints.force_paired:
        return True
    paired_mir = {j for i, j in non_seed_pairs}
    for force_paired in constraints.force_paired:
        if set.union(paired_mir, force_paired) == paired_mir:
            return True
    return False


def _get_gu_pairs(seq, non_seed_pairs):
    '''
    Get miRNA bases which form GU pairs to the target.
    '''
    tar_seq, mir_seq = seq.split('&')
    gu_pairs = [j for i, j in non_seed_pairs
                if tar_seq[-i]+mir_seq[j-1] in ('GU', 'UG')]
    return gu_pairs


def _check_gu_pairs(gu_pairs, constraints):
    '''
    Check if the structure meets the selection criteria for presence
    and location of GU wobble pairs.
    '''
    if not gu_pairs:
        return True
    if constraints.no_gu:
        return False
    if constraints.no_gu_bases:
        if set.intersection(set(gu_pairs), set(constraints.no_gu_bases)):
            return False
    return True


def _check_mismatches(non_seed_pairs, constraints):
    '''
    Check if the structure meets the selection criteria for number
    and size of mismatches in the 3P helix.
    '''
    pair_dists = [(abs(i[0]-j[0]), abs(j[1]-i[1]))
                  for i, j in zip(non_seed_pairs, non_seed_pairs[1:])]
    mismatches = [(i, j) for i, j in pair_dists if (i, j) != (1, 1)]
    if not mismatches:
        return True
    if len(mismatches) > constraints.max_mm_count:
        return False
    if (max((i-1 for i in chain(*mismatches)), default=0)
            > constraints.max_mm_size):
        return False
    return True
