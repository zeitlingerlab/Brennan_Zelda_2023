"""
Melanie Weilert
January 2021
Purpose: Extract PPM, CWM, and hCWMs from each TF-MoDISco pattern


Intructions:
1. Activate `bpnet` conda environment.
2. `python /n/projects/mw2098/shared_code/bpnet/bpnet_export_large_bw.py [options]`
3. For help: `python /n/projects/mw2098/shared_code/bpnet/bpnet_extract_modisco_logo_information.py -h`
"""

# Setup
import os
import json
import pandas as pd
import numpy as np
from optparse import OptionParser
from tqdm import tqdm
from bpnet.modisco.files import ModiscoFile
from bpnet.modisco.utils import trim_pssm_idx

import warnings
warnings.filterwarnings("ignore")

parser = OptionParser()
parser.add_option("-m", "--modisco_model_file",
                  help="Path to TF-MoDISco model.h5 file")
parser.add_option("-o", "--output_h5",
                  help="Path to output .h5 file")
parser.add_option("-i", "--ic_threshold", default = 0.08, type = 'float',
                  help="Information content threshold for trimming seqlet. [default = %default]")
(options, args) = parser.parse_args()

def bpnet_extract_modisco_logos(modisco_model_file,
                                output_h5,
                                ic_threshold = 0.08):
    """Extract PPM, CWM, and hCWMs from each TF-MoDISco pattern
    """

    #Function that takes in a task-separated dictionary of ModiscoFile objects and extracts all patterns
    def get_patterns(mf):
        patterns = []
        tasks = mf.tasks()
        for pattern in mf.patterns():
            pattern_name = pattern.name
            n_seqlets = mf.n_seqlets(pattern_name)
            p = mf.get_pattern(pattern_name)
            p.attrs = {'contribs': tasks, 'n_seqlets': n_seqlets}
            # for t in tasks:
            #     if t in p.contrib:
            #         continue
            #     else:
            #         p.contrib[t] = p.contrib[task]
            #         p.hyp_contrib[t] = p.hyp_contrib[task]
            # p._tasks = tasks
            patterns.append(p.copy().rename(p.name))
        return patterns

    #Import ModiscoFile
    mf = ModiscoFile(modisco_model_file)

    #Get all patterns
    all_patterns = get_patterns(mf = mf)

    pattern_info_dict = {}
    for pattern in tqdm(all_patterns):
        pattern_contribs = pattern.attrs['contribs']
        pattern_info_dict[f'{pattern.name}'] = {
            'contribs' : pattern_contribs,
            'name' : pattern.name,
            'trim_coords': {'ic_threshold': ic_threshold, 'coords': trim_pssm_idx(pattern.get_seq_ic(), frac=ic_threshold)},
            'PPM' : {'full': pattern.seq.tolist(), 'trimmed': pattern.trim_seq_ic(trim_frac=ic_threshold).seq.tolist()},
            'CWM' : {'full': {k: v.tolist() for k,v in pattern.contrib.items()},
                     'trimmed': {k: v.tolist() for k,v in pattern.trim_seq_ic(trim_frac=ic_threshold).contrib.items()}},
            'hCWM' : {'full': {k: v.tolist() for k,v in pattern.hyp_contrib.items()},
                     'trimmed': {k: v.tolist() for k,v in pattern.trim_seq_ic(trim_frac=ic_threshold).hyp_contrib.items()}}
        }

    from silx.io.dictdump import dicttoh5
    dicttoh5(pattern_info_dict, output_h5)

    return None


bpnet_extract_modisco_logos(modisco_model_file = options.modisco_model_file,
                            output_h5 = options.output_h5,
                            ic_threshold = options.ic_threshold)
