
"""
Melanie Weilert
April 2020
Purpose: Given a set of motifs and corresponding sequences/predictions, generate in silico genomic perturbations.

This will measure the predictions across the entire/a portion of the predicted window,
the pseudocounts (pc), max (relative to reference),
and sum of signal across all tasks across all motif windows for each mutation.

Intructions:
1. Activate `bpnet` or `bpnet-gpu` conda environment.
2. `python /n/projects/mw2098/shared_code/bpnet/bpnet_generate_perturbations.py [options]`
3. For help: `python /n/projects/mw2098/shared_code/bpnet/bpnet_generate_perturbations.py -h`
"""

# Setup
import sys
import os
import pyBigWig
import logging
import subprocess
import pandas as pd
import numpy as np
import tensorflow as tf
from optparse import OptionParser
from itertools import compress
from operator import itemgetter
from collections import OrderedDict
from tqdm import tqdm
from pybedtools import BedTool
from pathos.multiprocessing import ProcessingPool as Pool
from bpnet.cli.contrib import bpnet_contrib
from bpnet.cli.modisco import cwm_scan
from bpnet.utils import add_file_logging

#Get custom BPNet functions
import sys
sys.path.insert(0, f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/analysis/scripts/py')
from bpnet_perturb_functions import random_seq_onehot, generate_alt_sequences

import warnings
warnings.filterwarnings("ignore")
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#Set up options
parser = OptionParser()
parser.add_option("-d", "--dfi",
                  help="Path to the .tsv.gz file of the mapped motifs of interest. This dfi should follow a `bpnet cwm-scan` format.")
parser.add_option("-m", "--model_dir",
                  help="Path to the trained model directory (specified in `bpnet train <output_dir>`")
parser.add_option("-c", "--contrib_file",
                  help="Filepath to BPNet Contrib .h5 file that matches the `example_idx` found in dfi")
parser.add_option("-o", "--output_prefix",
                  help="Output prefix for the chromosome-separated .tsv.gz files.")
parser.add_option("-u", "--comb_max", default = None, type = "int",
                  help="Upper limit of simultaneous mutations to conduct across each region. [default = None]")
parser.add_option("-l", "--comb_min", default = None, type = "int",
                  help="Lower limit of simultaneous mutations to conduct across each region. [default = None]")
parser.add_option("-n", "--nodes", default = 6, type = "int",
                  help="Parallel workers for generating perturbation sequences.[default = 6]")
parser.add_option("-t", "--trials", default = 16, type = "int",
                  help="Number of trials for generating perturbed motifs.[default = 16]")
parser.add_option("--use_whole_window", default = False, action="store_true",
                  help="If selected, then the code will use the whole prediction window instead of the designated `--summary_window` parameter. This option overwrites the `--summary_window` approach. [default = False]")
parser.add_option("-w", "--summary_window", default = 50, type = "int",
                  help="Window around motif to compute the max/sum summary predictions. [default = 50]")
parser.add_option("-p", "--pseudo_count_quantile", default = .2, type = "int",
                  help="Threshold to compute pseudo count values across each window for each task/mutant. [default = .2]")
parser.add_option("-g", "--gpu", default = None,
                  help="which gpu to use.  [default = %default]")
parser.add_option("-x", "--memfrac-gpu", default = .45, type = "float",
                  help="what fraction of the GPU memory to use [default = %default]")
(options, args) = parser.parse_args()

def bpnet_generate_perturbations(dfi, contrib_file, model_dir, output_prefix,
                           comb_max = None, comb_min = None,
                           nodes = 6, trials = 16, use_whole_window = False,
                           summary_window = 51, pseudo_count_quantile = .2,
                           gpu = None, memfrac_gpu = .8):

    from functools import reduce
    from pathos.multiprocessing import ProcessingPool as Pool
    from bpnet.BPNet import BPNetSeqModel
    from itertools import combinations, chain
    from bpnet.cli.contrib import ContribFile
    from bpnet.utils import create_tf_session

    # Make directories
    output_dir = os.path.dirname(output_prefix)
    add_file_logging(output_dir, logger, 'bpnet-generate-perturbations')
    os.makedirs(output_dir, exist_ok=True)

    #Initialize GPU
    if gpu is not None:
        create_tf_session(gpu, per_process_gpu_memory_fraction=memfrac_gpu)

    #Get contribution scores
    contrib = ContribFile(file_path = contrib_file)

    #Get reference sequences
    ref_seqs_all = contrib.get_seq()

    #Get unique motif names
    dfi = pd.read_csv(dfi, sep = '\t')
    dfi = dfi.sort_values(['example_idx','pattern_name'])
    dfi['pattern_name_unique'] = dfi['pattern_name'] + '-' + (dfi.groupby(['example_idx','pattern_name']).cumcount()+1).astype(str)

    #Load the model for predictions
    bpnet_model = BPNetSeqModel.from_mdir(model_dir) #30s
    tasks = bpnet_model.tasks

    #For each chromosome, subset dfi, generate sequences, predict, summarize, and export.
    for chrom in dfi['example_chrom'].unique():

        logger.info('Generating perturbations from ' + chrom + '...')

        #Subset regions based on chromosome
        example_idxs = dfi[dfi['example_chrom']==chrom]['example_idx'].unique()

        #Set up lists of matching idxs, seqs, and motifs
        ref_seqs = [ref_seqs_all[i] for i in example_idxs]
        motif_sets = [dfi[dfi['example_idx'].astype(int)==i] for i in example_idxs]
        windows_n = len(example_idxs)

        #Initialize pool of nodes
        p = Pool(nodes=nodes)

        #For x trials, (1) generate sequences, (2) predict from sequences
        preds_by_trial = []

        for trial in range(trials):
            #Return list of perturbed sequences for each region...list[mut->1000x4, mut->1000x4, ...]
            perturb_seqs_list = p.map(generate_alt_sequences, ref_seqs, motif_sets,
                                      [comb_max]*windows_n, [comb_min]*windows_n)
            s = []
            pis_df = pd.DataFrame(columns = ['mut','example_idx'])
            for window in range(windows_n):
                #Collect the sequences
                perturb_seq = perturb_seqs_list[window]
                perturb_seq_array = np.array([v for v in perturb_seq.values()])
                s.append(perturb_seq_array)

                #Record the indexes
                pi = pd.DataFrame([k for k in perturb_seq.keys()], columns = ['mut'])
                pi['example_idx'] = example_idxs[window]
                pis_df = pis_df.append(pi)
            pis_df['perturb_idx'] = range(pis_df.shape[0])
            s = np.vstack(s)

            perturb_preds = bpnet_model.predict(s)
            preds_by_trial.append(perturb_preds) #list[task-> method -> region x window x channel]

        #Average trials together: task -> region x window x strand where region = each perturb_idx of pis_df
        logger.info('Collecting perturbations together' + chrom + '...')
        perturb_preds = {task: np.mean(np.array([preds_by_trial[trial][task]
                                                 for trial in range(trials)]), axis = 0)
                         for task in tasks}

        logger.info('Summarizing perturbations from ' + chrom + '...')
        summary_list = [] #append `summary_dict_by_task` for each task
        counter = 0

        for task in tasks:
            logger.info('Summarizing perturbations across ' + chrom + ', ' + task + '...')
            summary_dict_by_task = {}
            for perturb_idx in pis_df.perturb_idx:
                #Subset to the correct prediction
                perturb_pred = perturb_preds[task][perturb_idx]

                #Collect associated information
                example_idx = pis_df['example_idx'][pis_df.perturb_idx==perturb_idx]
                motif_set = dfi[dfi['example_idx']==example_idx.values[0]]
                mut = pis_df.mut[pis_df.perturb_idx==perturb_idx]

                #Get pseudocounts across the predicted region
                pc = np.percentile(perturb_pred, pseudo_count_quantile * 100)

                if use_whole_window:
                    #logger.info('--use_whole_window marked, summarizing entire window...' + chrom + '...')
                    #Summarize the whole window predictions
                    pred_sum = np.sum(np.abs(perturb_pred))
                    summary_dict_by_task[counter] = [perturb_idx, example_idx.values[0], mut.values[0], pc, pred_sum,]

                    counter += 1
                else:
                    #logger.info('Summarizing ' + summary_window +  ' bp around each motif...')
                    #Collect the associated reference prediction
                    ref_perturb_idx = pis_df[(pis_df['example_idx']==example_idx.values[0]) & (pis_df['mut']=='Reference')].perturb_idx
                    ref_perturb_pred = perturb_preds[task][ref_perturb_idx.values[0]]

                    #For each motif in a motif_set...
                    summary_across_motif_set = pd.DataFrame()

                    for idx,motif in motif_set.iterrows():

                        #Collect boundaries of current motif
                        centerpoint = int(motif['pattern_start'] + np.floor((motif['pattern_end'] - motif['pattern_start'])/2))
                        lower_bounds = int(centerpoint - np.floor(summary_window/2))
                        upper_bounds = int(centerpoint + np.floor(summary_window/2))

                        #If boundaries are exceeded beyond the window, create null entry
                        if lower_bounds<=0 or upper_bounds>=bpnet_model.input_seqlen():
                            summary_dict_by_task[counter] = [perturb_idx, example_idx.values[0], mut.values[0], motif.pattern_name_unique,
                                                                 pc, None, None]
                        else:
                            #Call ref sig across motif slice within the summary_window and find max on pos and neg strand
                            max_sig_idx = np.argmax(ref_perturb_pred[lower_bounds:upper_bounds], axis = 0)

                            #Collect values across regions of interest.
                            sig_bounded = perturb_pred[lower_bounds:upper_bounds,:] #window_len x strand
                            pred_sum = np.sum(np.abs(sig_bounded))
                            pred_max = sig_bounded[max_sig_idx,[0, 1]].sum()

                            summary_dict_by_task[counter] = [perturb_idx, example_idx.values[0], mut.values[0], motif.pattern_name_unique,
                                                                 pc, pred_sum, pred_max]
                        counter += 1

            #Collect dict objects by task and convert to df
            if use_whole_window:
                summary_df = pd.DataFrame.from_dict(summary_dict_by_task,orient='index')
                summary_df.columns = ['perturb_idx', 'example_idx', 'mut', f'{task}/pc', f'{task}/pred_sum']
            else:
                summary_df = pd.DataFrame.from_dict(summary_dict_by_task,orient='index')
                summary_df.columns = ['perturb_idx', 'example_idx','mut','motif',
                                          f'{task}/pc', f'{task}/pred_sum', f'{task}/pred_max']
            #Append to list
            summary_list.append(summary_df)

        #Collect different task columns and assign to single df
        logger.info('Collecting all perturbs into pd.df...')
        if use_whole_window:
            summary_df = reduce(lambda x,y: pd.merge(x,y, on=['perturb_idx','example_idx','mut'], how='outer'),
                                summary_list)
        else:
            summary_df = reduce(lambda x,y: pd.merge(x,y, on=['perturb_idx','example_idx','mut','motif'],
                                                     how='outer'),
                                summary_list)

        logger.info('Writing perturbs to a .tsv.gz...')
        summary_df.to_csv(f'{output_prefix}_{chrom}.tsv.gz', sep = '\t', index = False)
    return None

#Run featured function
bpnet_generate_perturbations(dfi = options.dfi,
                                contrib_file = options.contrib_file,
                                model_dir = options.model_dir,
                                output_prefix = options.output_prefix,
                                comb_max = options.comb_max,
                                comb_min = options.comb_min,
                                nodes = options.nodes,
                                trials = options.trials,
                                use_whole_window = options.use_whole_window,
                                summary_window = options.summary_window,
                                pseudo_count_quantile = options.pseudo_count_quantile,
                                gpu = options.gpu,
                                memfrac_gpu = options.memfrac_gpu)
