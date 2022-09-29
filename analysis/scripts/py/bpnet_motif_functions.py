# Functions and classes associated with generating perturbations
import sys
import os
import warnings
import pandas as pd
import numpy as np
import random
import math
import itertools
import matplotlib.pyplot as plt
from collections import Counter
from pybedtools import BedTool
from bpnet.cli.contrib import ContribFile
from bpnet.BPNet import BPNetSeqModel
from bpnet.data import NumpyDataset #creates a set of nested arrays
from kipoi.data import Dataset

warnings.filterwarnings("ignore")
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False

def myround(x, base=50):
    return base * round(x/base)
def myceiling(x, base = 10):
    return int(math.ceil(x / base)) * base
def myfloor(x, base = 10):
    return int(math.floor(x / base)) * base

def one_hot_decode_sequence(array):
    """
    Purpose: Given an array [position x 4], decode sequence to a string.
    """
    onehot_decoder = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T'
    }

    idxs = np.where(array)[1]
    return (''.join([onehot_decoder[i] for i in idxs]))

def one_hot_encode_sequence(sequence):
    """
    Kudos to Charles: /n/projects/cm2363/bpnet-nucleosomes/work/localimportance/allLocalImportances.py
    Purpose: Given a SINGLE sequence string, one-hot encode the data.
        + default control_profiles and control_logcounts is to be zeroed out
        + naively detects whether the sequence is one-hot-encoded.
    """
    onehot_mapping = {
    'A': [1,0,0,0],
    'C': [0,1,0,0],
    'G': [0,0,1,0],
    'T': [0,0,0,1],
    'a': [1,0,0,0],
    'c': [0,1,0,0],
    'g': [0,0,1,0],
    't': [0,0,0,1],
    'N': [0,0,0,0]
    }
    return np.array([onehot_mapping[x] for x in sequence])

def one_hot_encode_sequences(sequences):
    """
    Purpose: Given an array of sequences, one-hot-encode into a [region x position x 4] array.
    """
    return(np.stack([one_hot_encode_sequence(s) for s in sequences]))


"""
Before resizing, remove regions that will extend over the edges of chromsome boundaries.
Given:
    regions: pd.Df or pybedtools.BedTool of regions
    chrom_sizes: path to chrom.sizes.txt file for model organism of `regions`
    boundary_width: integer specifying how far to scan on each chromosome edge
Returns: regions pd.DataFrame object with filtered regions satisfying boundary criteria
"""

def check_chromosome_boundaries(regions, chrom_sizes, boundary_width = 1000):
    import pandas as pd
    from pybedtools import BedTool

    assert isinstance(regions, pd.DataFrame) or isinstance(regions, BedTool), 'Input needs to be a pybedtools.BedTool object or pd.DataFrame'

    chrom_sizes_df = pd.read_csv(chrom_sizes, sep = '\t', header = 0, names = ['chrom','length'])

    if isinstance(regions, BedTool):
        regions = BedTool.to_dataframe(regions)

    regions = regions[regions['start']>=boundary_width]
    regions = regions.merge(chrom_sizes_df, on = 'chrom')
    regions = regions[(regions['end'] + boundary_width)<=regions['length']]
    return(regions)

"""
Palindromic motifs will overlap on the + and - strand, but not necessarily be the exact same coordinates.
Keep the motif with a higher contribution score.
Given:
    dfi: list of motifs with CWM-scanned convention from `bpnet cwm-scan` [requires `row_idx`, `pattern_start`, `pattern_end` and `pattern_name`]
    pattern: string of pattern to filter for under `pattern_name` column of dfi
    overlap_threshold: %age of overlap that returns an unacceptable motif.
Returns: list of `row_idx` values to remove that are redundant.
"""

def remove_palindromic_motif_duplicates(dfi, pattern, overlap_threshold = .7,
                        motif_index_column = 'row_idx', motif_window_column = 'example_idx',
                        motif_start_column = 'pattern_start', motif_end_column = 'pattern_end',
                        motif_name_column = 'pattern_name', motif_length_column = 'pattern_len',
                        motif_signal_column = 'contrib_max'):
    remove_row_idx = np.array([])

    def getOverlap(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    i = dfi[dfi[motif_name_column]==pattern]

    #for each region
    for idx in i.example_idx.unique():
        i_across_r = i[i[motif_window_column]==idx]

        #Define windows to check for overlaps
        window_interv_dict = {row[motif_index_column]: (row[motif_start_column],row[motif_end_column])
                              for idx, row in i_across_r.iterrows()}

        #Define all combinations of windows to test
        overlap_combs = list(itertools.combinations(window_interv_dict.keys(),2))

        #Test each combination, mark the one with a lower contribution score for removal
        for comb in overlap_combs:
            ov_count = getOverlap(window_interv_dict[comb[0]], window_interv_dict[comb[1]])
            if ov_count/i[motif_length_column].iloc[0]: #If there is sufficient, overlap, select comb to remove.
                if i[motif_signal_column][i[motif_index_column]==comb[0]].iloc[0] > i[motif_signal_column][i[motif_index_column]==comb[1]].iloc[0]:
                    remove_row_idx = np.append(remove_row_idx, comb[1])
                else:
                    remove_row_idx = np.append(remove_row_idx, comb[0])

    return remove_row_idx

"""
Sometimes motifs will be subsets of other motifs (e.g. Dl vs a Dl-Dl). In this case, CWM-scanning will annotate both.
This can pose as a problem for some analysis, so this function will remove overlapping motifs based on a priority list.
Given:
    dfi: pd.df of motifs with CWM-scanned convention from `bpnet cwm-scan` [requires `row_idx`, `pattern_start`, `pattern_end` and `pattern_name`]
    priority_list: list of motifs ordered most->least important in terms of priority that are expected to overlap.
Returns: filtered pd.df with same format of `dfi` that will contain only the highest priority motifs in the event of an overlap
"""

def filter_overlapping_motifs_by_priority(dfi, priority_list):
    from tqdm import tqdm
    import itertools
    #Separate patterns that don't match priority
    def _grep(x, list_to_grep):
        """takes a string (x) and checks whether any string from a given
           list of strings (list_to_grep) exists in `x`"""
        for text in list_to_grep:
            if text.lower() in x.lower():
                return True
        return False

    contains_motifs = dfi.pattern_name.apply(_grep, list_to_grep=priority_list)

    #Separate based on which are motifs of interest in the "priority casacde".
    dfi_with_motifs = dfi[contains_motifs]
    dfi_without_motifs = dfi[~contains_motifs]

    row_idx_to_keep = []
    for ei in tqdm(dfi_with_motifs.example_idx.unique()):
        dfi_across_ei = dfi_with_motifs[dfi_with_motifs.example_idx == ei]

        #Sort values by priority rank and pattern center
        motif_sorter = dict(zip(priority_list,range(len(priority_list))))
        dfi_across_ei['pattern_name_rank'] = dfi_across_ei['pattern_name'].map(motif_sorter)
        dfi_across_ei.sort_values(['pattern_name_rank', 'contrib_weighted'], \
                ascending = [True, True], inplace = True)
        dfi_across_ei.drop('pattern_name_rank', 1, inplace = True)

        #Define function to compute overlaps
        def getOverlap(a, b):
            return max(0, min(a[1], b[1]) - max(a[0], b[0]))

        #For each motif (except the first one), check how it overlaps with higher priority motifs.
        #Mark if overlap occurs
        remove_bool = [False]
        for i in range(1,dfi_across_ei.shape[0]):
            ov = []
            for j in range(i):
                o = getOverlap(a = (dfi_across_ei.pattern_start.iloc[i], dfi_across_ei.pattern_end.iloc[i]),
                               b = (dfi_across_ei.pattern_start.iloc[j], dfi_across_ei.pattern_end.iloc[j]))
                ov.append(o)
            if sum(ov)>0:
                drop = True
            else:
                drop = False
            remove_bool.append(drop)

        dfi_across_ei['drop']=remove_bool
        dfi_across_ei = dfi_across_ei[~dfi_across_ei['drop']]
        row_idx_to_keep.append(np.array(dfi_across_ei['row_idx']))

    row_idx_to_keep = list(itertools.chain.from_iterable(row_idx_to_keep))
    dfi_with_motifs_filt = dfi_with_motifs[dfi_with_motifs['row_idx'].isin(row_idx_to_keep)]
    dfi_to_return = dfi_with_motifs_filt.append(dfi_without_motifs)
    return dfi_to_return
