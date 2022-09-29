
import os
import warnings
import pandas as pd
import numpy as np
import random
import math
import itertools
import matplotlib.pyplot as plt
import tensorflow as tf
from itertools import compress
from operator import itemgetter
from collections import OrderedDict
from tqdm import tqdm
from pathos.multiprocessing import ProcessingPool as Pool
from collections import Counter
from pybedtools import BedTool
from bpnet.cli.contrib import ContribFile
from bpnet.BPNet import BPNetSeqModel
from bpnet.data import NumpyDataset #creates a set of nested arrays
from kipoi.data import Dataset

#Source custom scripts
import sys
sys.path.insert(0, f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/analysis/scripts/py')
from bpnet_plot_functions import plot_pred_motif_boxes, plot_nexus_prediction, plot_contrib_track, plot_contrib_motif_box
from bpnet_data_format_functions import pred_npd_to_tidy_df, contrib_npd_to_df

warnings.filterwarnings("ignore")
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False

def myround(x, base=50):
    return base * round(x/base)
def myceiling(x, base = 10):
    return int(math.ceil(x / base)) * base
def myfloor(x, base = 10):
    return int(math.floor(x / base)) * base

"""
Purpose: Generate a random one-hot-encoded DNA sequence
Given:
    l: integer value of how long the random sequence should be
Output:
    seq_1he = np.array with shape [l x 4]
"""

def random_seq_onehot(l):
    """Generate random sequence one-hot-encoded
    Args:
      l: sequence length
    """
    from concise.preprocessing import encodeDNA
    import random

    return encodeDNA([''.join(random.choices("ACGT", k=int(l)))])[0]

"""
Purpose: Generate all combinations of mutated sequences based on given sequence and desired mutant coordinates.
Given:
    ref_seq: np.array with shape [l x 4] of the one-hot encoded reference sequence
    motifs: pd.df of the motifs across the reference sequence
        + format based on `bpnet cwm-scan`
        + required columns: `pattern_len`, `pattern_start`, `pattern_end`,`pattern_name_unique`
    comb_max: max number of simultaneous mutations allowed
    comb_min: min number of simultaneous mutations allowed
Output:
    perturb_seqs_dict_list = trials x -> ref/mutation -> seqlen x 4
"""

def generate_alt_sequences(ref_seq, motifs, comb_max = None, comb_min = None):

    from itertools import combinations, chain

    def generate_random_sequences_1he(l):
        """Generate random sequence one-hot-encoded
        Args:
          l: sequence length
        """
        from concise.preprocessing import encodeDNA
        import random

        return encodeDNA([''.join(random.choices("ACGT", k=int(l)))])[0]

    #Insert relevant mutant sequences in each combination
    def get_alt_seq_from_combo(combo, ref_seq, motifs, mut_seqs_dict):
        alt_seq = ref_seq.copy()
        for mut in combo:
            motif = motifs[motifs['pattern_name_unique']==mut]
            alt_seq[motif.pattern_start.iloc[0]:motif.pattern_end.iloc[0]] = mut_seqs_dict[mut]
        return alt_seq

    #Mark the max depth of combinations
    c_min_depth = 1
    c_max_depth = motifs.shape[0]
    if(comb_max is not None): c_max_depth = comb_max
    if(comb_min is not None): c_min_depth = comb_min

    #Determine all combinations of unique motifs
    mut_combos = [list(combinations(motifs['pattern_name_unique'], d)) for d in range(c_min_depth, c_max_depth+1)]
    mut_combos = list(chain(*mut_combos))

    #Record mutant sequences used in this combination: motif -> motif_len x 4
    mut_seqs_dict = {motifs['pattern_name_unique'].iloc[idx]:
                          generate_random_sequences_1he(motifs['pattern_len'].iloc[idx])
                     for idx in range(len(motifs))}

    # Generate: combo -> 1000 x 4
    alt_seqs_dict = {'_'.join(list(combo)): get_alt_seq_from_combo(combo = combo, ref_seq = ref_seq,
                                                                   motifs = motifs, mut_seqs_dict = mut_seqs_dict)
                     for combo in mut_combos}
    #Add reference sequence
    perturb_seqs_dict = {**{'Reference': ref_seq}, **alt_seqs_dict}

    #Clean up results
    #motif_combo_ids = list(perturb_seqs_dict.keys())
    #perturb_seq_array = np.array([v for k,v in perturb_seqs_dict.items()])
    #return motif_combo_ids, perturb_seq_array
    return(perturb_seqs_dict)




"""
Purpose: Generate perturbation profiles across motifs
Given:
    + dfi: pd.df of the motifs across a single reference region
        + format based on `bpnet cwm-scan`
        + required columns: `pattern_len`, `pattern_start`, `pattern_end`,`pattern_name_unique`
    + model_dir = Path to the trained model directory (specified in `bpnet train <output_dir>`
    + contrib_file = BPNet Contrib .h5 file that matches the `example_idx` found in dfi
    + contrib_type = type of contribution score to extract ('profile' or 'count')
    + comb_max = Upper limit of simultaneous mutations to conduct across each region. [default = None]
    + comb_min = Lower limit of simultaneous mutations to conduct across each region. [default = None]
    + trials = Number of trials for generating perturbed motifs.[default = 16]
    + pseudo_count_quantile = Threshold to compute pseudo count values across each window for each task/mutant. [default = .2]
    + gpu = which gpu to use.  [default = %default]
    + memfrac = what fraction of the GPU memory to use [default = %default]
    + return_contrib = whether to compute contribution scores [default = %default]
Outputs:
    + perturb_preds = dictionary of profile predictions with the format: task -> [perturb x seqlen x strand]
    + perturb_contrib = dictionary of contribution scores with the format: task -> [perturb x seqlen x 4]
    + muts = list of perturbation combinations that were conducted
        + len(muts)==perturb_preds[task].shape[0]==perturb_contrib[task].shape[0]
"""

def generate_perturbs_across_window(dfi, contrib_file, model_dir,
                                    comb_max = None, comb_min = None,
                                    trials = 16, contrib_type = 'profile',
                                    pseudo_count_quantile = .2,
                                    gpu = None, memfrac_gpu = .8,
                                    return_contrib = True):

    from bpnet.BPNet import BPNetSeqModel
    from itertools import combinations, chain
    from bpnet.cli.contrib import ContribFile
    from bpnet.utils import create_tf_session

    #Import custom bpnet perturbation functions.
    import sys
    sys.path.insert(0, f'/n/projects/mw2098/shared_code/bpnet/scripts')
    from perturb_functions import random_seq_onehot, generate_alt_sequences

    #Ensure that motifs are from a single BPNet window
    assert len(dfi.example_idx.unique())==1, 'Perturb only a single region. Please subset `dfi`.'

    #Initialize GPU
    if gpu is not None:
        create_tf_session(gpu, per_process_gpu_memory_fraction=memfrac_gpu)

    #Get reference sequence
    contrib = ContribFile(file_path = contrib_file)
    ref_seq = contrib.get_seq(idx = dfi.example_idx.unique()[0])

    #Load the model for predictions
    bpnet_model = BPNetSeqModel.from_mdir(model_dir) #30s
    tasks = bpnet_model.tasks

    # For x trials, (1) generate sequences, (2) predict from sequences
    seqs_by_trial = []
    preds_by_trial = []
    contrib_by_trial = []
    for trial in range(trials):
        # Return list of perturbed sequences and convert into an array [mut x seqlen x 4]
        perturb_seqs_dict = generate_alt_sequences(ref_seq = ref_seq, motifs = dfi,
                                                   comb_max = comb_max, comb_min = comb_min)
        perturb_seqs_array = np.array([v for v in perturb_seqs_dict.values()])

        # Predict array of seq/preds/contrib
        if return_contrib:
            perturb_all = bpnet_model.predict_all(perturb_seqs_array, contrib_method = 'deeplift')

            # Seq: [perturb x window x 4]
            seq = np.array([perturb_all[i]['seq']
                            for i in range(len(perturb_seqs_dict))])
            seqs_by_trial.append(seq)

            # Preds: task->[perturb x window x strand]
            preds = {task: np.array([perturb_all[i]['pred'][task]
                                     for i in range(len(perturb_seqs_dict))])
                     for task in tasks}
            preds_by_trial.append(preds)

            # Contrib: task->[perturb x window x strand]
            contrib = {task: np.array([perturb_all[i]['contrib_score'][f'{task}/{contrib_type}']
                                     for i in range(len(perturb_seqs_dict))])
                       for task in tasks}
            contrib_by_trial.append(contrib)
        else:
            perturb_all = bpnet_model.predict(perturb_seqs_array)

            # Preds: task->[perturb x window x strand]
            preds = {task: np.array([perturb_all[task][i]
                                     for i in range(len(perturb_seqs_dict))])
                     for task in tasks}
            preds_by_trial.append(preds)

    # Collect values across trials. For contrib, multiply by 1he sequence to keep contrib scores clean.
    # Preds: task -> [perturb x seqlen x strand]
    perturb_preds = {task: np.mean(np.array([preds_by_trial[trial][task]
                                             for trial in range(trials)]), axis = 0)
                     for task in tasks}

    # Contrib: task -> [perturb x seqlen x 4]
    if return_contrib:
        perturb_contrib = {task: np.mean(np.array([contrib_by_trial[trial][task]*seqs_by_trial[trial]
                                                 for trial in range(trials)]), axis = 0)
                         for task in tasks}

        return perturb_preds, perturb_contrib, list(perturb_seqs_dict.keys())
    else:
        return perturb_preds, [], list(perturb_seqs_dict.keys())



"""
Purpose: Plot perturbation profiles across motifs
Given:
    + perturb_preds = dictionary of profile predictions with the format: task -> [perturb x seqlen x strand]
    + perturb_contrib = dictionary of contribution scores with the format: task -> [perturb x seqlen x 4]
    + muts = list of perturbation combinations that were conducted
        + len(muts)==perturb_preds[task].shape[0]==perturb_contrib[task].shape[0]
    + dfi: pd.df of the motifs across a single reference region
        + format based on `bpnet cwm-scan`
        + required columns: `pattern_len`, `pattern_start`, `pattern_end`,`pattern_name_unique`, `pattern`
    + task_color_dict: dictionary of tasks specifying which color to plot across tasks.
        + If tasks are omitted, then they will not be plotted.
    + muts_to_plot: list of mutations to subset and plot. Default = None will plot all possible mutations.
    + title: title of plot [default = None]
    + xrange: coordinate list [min, max] of window subset
    + free_ylim: Boolean indicating whether ylimits should be freely scaled [default = False]
Outputs:
    + fig: matplotlib figure
"""
def plot_perturbs_across_window(perturb_preds, perturb_contrib, muts, dfi,
                                task_color_dict, muts_to_plot = None, title = None,
                                xrange = [0, 1000], free_ylim = False):
    import itertools
    import matplotlib.pyplot as plt
    from bpnet.plot.tracks import plot_track, to_neg#, plot_seqlet_box
    from bpnet.modisco.pattern_instances import dfi2seqlets
    from plot_functions import plot_seqlet_box

    def roundup(num, dec):
        return(np.ceil(num * dec) / dec)
    def rounddown(num, dec):
        return(np.floor(num * dec) / dec)

    assert len(dfi.example_idx.unique())==1, 'Please only plot dfi across a single window.'

    if muts_to_plot==None:
        muts_to_plot=muts

    # Define figure format and dimensions
    muts_to_plot_idxs = [muts.index(mut) for mut in muts_to_plot]
    xslice = slice(xrange[0],xrange[1])
    num_cols = len(task_color_dict)
    num_conds = len(muts_to_plot)
    num_rows = (num_conds) *2 #Preds and contrib tracks
    gs_dict={'width_ratios': list(itertools.repeat(1, num_cols)),
             'height_ratios': list(np.array(list(itertools.repeat([2,1], num_conds))).flatten())
            } #2:1 ratio of pred to contrib

    fig, axs = plt.subplots(nrows = num_rows, ncols = num_cols,
                            figsize=(num_cols*4, num_conds*2.5),
                            gridspec_kw=gs_dict, constrained_layout=True) #sharex=True, figsize = (width, height)

    for t,task in enumerate(list(task_color_dict.keys())):
        task_pred_ylim = roundup(np.max(perturb_preds[task]), 200)
        task_contrib_ylim_upr = roundup(np.max(perturb_contrib[task]), 200)
        task_contrib_ylim_lwr = rounddown(np.min(perturb_contrib[task]), 200)

        for m_count,m in enumerate(muts_to_plot_idxs): #m_count is for plot indexes, m is for perturb idxs
            if len(list(task_color_dict.keys()))==1: #If you only plot 1 task, axs handles indexing differently.
                pred_idx_tuple = (2*m_count)
                contrib_idx_tuple = (2*m_count + 1)
            else:
                pred_idx_tuple = ((2*m_count),t)
                contrib_idx_tuple = ((2*m_count + 1),t)

            if free_ylim:
                pred_ylim = None
                contrib_ylim = None
            else:
                pred_ylim = (-task_pred_ylim, task_pred_ylim)
                contrib_ylim = (task_contrib_ylim_lwr,task_contrib_ylim_upr)

            #Plot prediction track
            plot_track(arr = to_neg(perturb_preds[task][m, xslice]), ax = axs[pred_idx_tuple],
                       ylim = pred_ylim, color = task_color_dict[task], track = None)
            axs[pred_idx_tuple].set_title(': '.join([task, muts[m]]))
            axs[pred_idx_tuple].margins(x=0)

            #Plot contribution track
            plot_track(arr = perturb_contrib[task][m, xslice], ax = axs[contrib_idx_tuple],
                       ylim = contrib_ylim)
            axs[contrib_idx_tuple].set_xticks(np.arange(min(xrange), max(xrange), 100), minor = False)
            axs[contrib_idx_tuple].margins(x=0)

            #Plot seqlets across window
            seqlets = dfi2seqlets(dfi, short_name=False)
            for s,seqlet in enumerate(seqlets):
                seqlet.name = dfi.pattern_name_unique.iloc[s] #Rename seqlets
                seqlet = seqlet.shift(-min(xrange)) #Shift seqlets to match reformatting window.
                plot_seqlet_box(seqlet, ax =  axs[pred_idx_tuple], add_label=False)
                plot_seqlet_box(seqlet, ax =  axs[contrib_idx_tuple], add_label=True)

                #Mark which seqlets are knocked out as red
                if seqlet.name in muts[m].split('_'):
                    plot_seqlet_box(seqlet, ax =  axs[pred_idx_tuple],
                                    box_fill = '#e50000', alpha = .25, add_label=False)
                    plot_seqlet_box(seqlet, ax =  axs[contrib_idx_tuple],
                                    box_fill = '#e50000', alpha = .25, add_label=False)
    fig.suptitle(title)
    return fig

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

def generate_random_seq(seqlen, weights = [.25, .25, .25, .25]):
    """
    Purpose: Generate a random DNA sequence of a specified length.
    """
    import random
    seq = random.choices(['A','C','G','T'], weights = weights, k=seqlen)
    return(''.join(seq))

def motif_coords(motif, position):
    """
    Kudos: I straight up jacked this from Ziga Avsec's bpnet code bpnet.bpnet.simulate#L19
    Purpose: Given motif (string) and a center position, find the motif boundaries.
    """
    start = position - len(motif) // 2
    end = start + len(motif)
    # print(start, end)
    return start, end

def insert_motif(seq, motif, position):
    """
    Kudos: I straight up jacked this from Ziga Avsec's bpnet code bpnet.bpnet.simulate#L25
    Purpose: Given a sequence, inject a motif centered on a position.
    """
    assert position < len(seq)
    start, end = motif_coords(motif, position)
    new_seq = seq[:start] + motif + seq[end:]
    assert len(new_seq) == len(seq)
    return new_seq

def predict_injected_seq(model,
                         primary_motif, primary_position,
                         secondary_motif = '', secondary_positions = [],
                         input_seqs = None,
                         input_seqlen = 1000,
                         output_seqlen = 1000,
                         trials = 512, batch_size=32, verbose = False):
    """
    Purpose:
        + Given a primary motif A, inject that motif and predict.
        + Given a secondary motif B, inject that motif with the primary motif as AB and measure synergy.
    Inputs:
        + model: BPNet model loaded by `BPNetSeqModel`
        + primary/secondary motif: string of sequence to inject
        + primary/secondary position: position within OUTPUT window to inject sequences
        + input_seq: If not None, should be a 1-hot encoded sequence array of random sequences to test across--overrides trials
        + input_seqlen: input sequence length
        + output_seqlen: output sequence length
        + trials: how many repeats to predict random sequences
    Outputs: dict of predicted values directly
    Application:
        1. Measure motif-motif syngery by assigning primary/secondary motifs as A and B, respectively.
            + `correct_for_shoulder_synergy` can be either False or True.
        2. Measure motif vs null distribution by assigning motif as secondary and null as primary.
            + `correct_for_shoulder_synergy` should be False.
        3. When you are doing analysis on factors with multiple channels per TF,
            you can use the output values to merge them after running this function.
    Math:
        1. Inject both motifs (AB) and primary motif (A) and measure AB/A
        2. Correction for shoulder syngery will take the AB component and correct for affinity of task to B.
            + AB - (B - 0), Where:
                - AB: contains both, primary and secondary motif
                - B : contains only secondary motif
                - 0 : doesn't contain any motif
    """

    from concise.preprocessing import encodeDNA

    #Compute where the input positions will be
    input_flank = (input_seqlen - output_seqlen) // 2
    input_primary_position = primary_position + input_flank
    #print("Generating sequences...")

    if input_seqs is not None:
        assert input_seqs.shape[2]==4, 'Array is not [region x position x 4]'
        assert input_seqs.shape[1]==input_seqlen, 'Array is not same as input_seqlen'
        null_seqs = input_seqs.copy()
        null_seqs = [one_hot_decode_sequence(input_seqs[s]) for s in range(input_seqs.shape[0])]
        trials = input_seqs.shape[0]
        if verbose:
            print(f'Number of trials is {trials}')
    else:
        null_seqs = [generate_random_seq(seqlen = input_seqlen) for t in range(trials)]

    A_seqs = [insert_motif(seq = t, motif = primary_motif, position = input_primary_position) for t in null_seqs]

    seqs = null_seqs + A_seqs

    for secondary_position in secondary_positions:

        input_secondary_position = secondary_position + input_flank

        #Generate repeated trials of sequence injections
        AB_seqs = [insert_motif(A, secondary_motif, input_secondary_position) for A in A_seqs]
        B_seqs = [insert_motif(N, secondary_motif, input_secondary_position) for N in null_seqs]

        #Add trials together embedded as null, A, [for each position -> B,AB]
        seqs = seqs + B_seqs + AB_seqs

    # One-hot encode sequences
    seqs_to_predict = encodeDNA(seqs)

    #predict profiles
    #print("Predicting...")
    preds_w_counts = model.predict(seqs_to_predict, batch_size = batch_size)

    #print("Collecting results...")
    #separate different injections and average across trials, which are the smallest margin
    pred_across_tasks_dict = {}
    for t in model.tasks:
        preds_w_counts_by_task = preds_w_counts[t]

        pred_dict = {}
        null_lower_range = 0 * trials
        null_upper_range = (0 + 1) * trials
        pred_dict['null'] = np.mean(preds_w_counts_by_task[null_lower_range:null_upper_range], axis = 0)

        A_lower_range = 1 * trials
        A_upper_range = (1 + 1) * trials
        pred_dict['A'] = np.mean(preds_w_counts_by_task[A_lower_range:A_upper_range], axis = 0)

        counter = 2
        labels = ['B', 'AB']
        distances_dict = {}
        for i,d in enumerate(secondary_positions):
            p_dict = {}
            for r in labels: #hard-coded because of the four states
                lower_range = counter * trials
                upper_range = (counter + 1) * trials
                p_dict[r]=np.mean(preds_w_counts_by_task[lower_range:upper_range], axis = 0)
                counter +=1
            distances_dict[d] = p_dict
        pred_dict['distances'] = distances_dict

        #Final appending across tasks
        pred_across_tasks_dict[t] = pred_dict

    return pred_across_tasks_dict


def summarize_injected_seq(preds_dict, primary_motif, primary_position,
                           secondary_motif = '', secondary_positions = [],
                           measurement_window_around_motif = 50, pseudocount_thresh = .2):
    """
    Purpose:
        + Given a pred_dict from predict_injected_seq from BPNet (extra subset of tasks), summarize values
        + TODO: Add maximum values
    Inputs:
        + model_chrombpnet: ATAC_seq coverage keras model imported by `load_chrombpnet` already
        + model_bias: Tn5 bias keras model imported by `load_chrombpnet` already
        + primary/secondary motif: string of sequence to inject
        + primary/secondary positions: position within OUTPUT window to inject sequences
        + measurement_window_around_motif: measure sum and maximum values in a window around primary motif (bp)
    Outputs: dict of predicted values directly
    Application:
        1. Measure motif-motif syngery by assigning primary/secondary motifs as A and B, respectively.
        2. Measure motif vs null distribution by assigning motif as secondary and null as primary.
        3. When you are doing analysis on factors with multiple channels per TF,
            you can use the output values to merge them after running this function.
    Math (not in function):
        1. Inject both motifs (AB) and primary motif (A) and measure AB/A
        2. Correction for shoulder syngery will take the AB component and correct for affinity of task to B.
            + AB - (B - 0), Where:
                - AB: contains both, primary and secondary motif
                - B : contains only secondary motif
                - 0 : doesn't contain any motif
    """
    perturb_df = pd.DataFrame()

    for t,pred_dict in preds_dict.items():
        windowA_min = primary_position - (measurement_window_around_motif // 2)
        windowA_max = primary_position + (measurement_window_around_motif // 2)

        results_dict = {}
        results_dict['task'] = t
        results_dict['pc_A'] = np.quantile(pred_dict['A'], pseudocount_thresh)
        results_dict['pc_null'] = np.quantile(pred_dict['null'], pseudocount_thresh)
        results_dict['all_sum_A'] = pred_dict['A'].sum()
        results_dict['all_sum_null'] = pred_dict['null'].sum()
        results_dict['windowA_sum_A'] = pred_dict['A'][windowA_min:windowA_max].sum()
        results_dict['windowA_sum_null'] = pred_dict['null'][windowA_min:windowA_max].sum()

        # all_motifA_max_pos = pred_dict['A'].argmax()
    #     windowA_motifA_max_pos = pred_dict['A'][windowA_min:windowA_max].argmax()
    #     results_dict['windowA_max_A'] = pred_dict['A'][windowA_min:windowA_max].max()
    #     results_dict['windowA_max_null'] = pred_dict['null'][windowA_min:windowA_max][windowA_motifA_max_pos].max()

        for d in secondary_positions:
            results_across_dist_dict = results_dict.copy()

            window_distance = d
            windowB_min = window_distance - (measurement_window_around_motif // 2)
            windowB_max = window_distance + (measurement_window_around_motif // 2)

            #Compute results across the whole window
            results_across_dist_dict['all_sum_AB'] = pred_dict['distances'][d]['AB'].sum()
            results_across_dist_dict['all_sum_B'] = pred_dict['distances'][d]['B'].sum()
            results_across_dist_dict['pc_B'] = np.quantile(pred_dict['distances'][d]['B'], pseudocount_thresh)

            #Compute results across partial window at motifA
            results_across_dist_dict['windowA_sum_AB'] = pred_dict['distances'][d]['AB'][windowA_min:windowA_max].sum()
            results_across_dist_dict['windowA_sum_B'] = pred_dict['distances'][d]['B'][windowA_min:windowA_max].sum()
    #         results_dict['windowA_max_AB'] = pred_dict['distances'][d]['AB'][windowA_min:windowA_max][windowA_motifA_max_pos]
    #         results_dict['windowA_max_B'] = pred_dict['distances'][d]['B'][windowA_min:windowA_max][windowA_motifA_max_pos]

            #Compute results across partial window at motifB
            results_across_dist_dict['windowB_sum_AB'] = pred_dict['distances'][d]['AB'][windowB_min:windowB_max].sum()
            results_across_dist_dict['windowB_sum_B'] = pred_dict['distances'][d]['B'][windowB_min:windowB_max].sum()
            results_across_dist_dict['windowB_sum_A'] = pred_dict['A'][windowB_min:windowB_max].sum()
            results_across_dist_dict['windowB_sum_null'] = pred_dict['null'][windowB_min:windowB_max].sum()
    #         windowB_motifB_max_pos = pred_dict['distances'][d]['B'][windowB_min:windowB_max].argmax()
    #         results_dict['windowB_max_AB'] = pred_dict['distances'][d]['AB'][windowB_min:windowB_max][windowB_motifB_max_pos].max()
    #         results_dict['windowB_max_A'] = pred_dict['distances'][d]['B'][windowB_min:windowB_max][windowB_motifB_max_pos].max()
    #         results_dict['windowB_max_B'] = pred_dict['A'][windowB_min:windowB_max].max()
    #         results_dict['windowB_max_null'] = pred_dict['null'][windowB_min:windowB_max][windowB_motifB_max_pos]

            #Record metadata
            results_across_dist_dict['primary_motif'] = primary_motif
            results_across_dist_dict['primary_position'] = primary_position
            results_across_dist_dict['secondary_motif'] = secondary_motif
            results_across_dist_dict['secondary_position'] = window_distance

            r_df = pd.DataFrame(results_across_dist_dict, index=[0])
            perturb_df = perturb_df.append(r_df)
    return(perturb_df)



















"""
########################################################################################
########################################################################################
########################################################################################

All functions below this point are deprecated for only single-mutation perturbation plots.
Use functions above for any-combination perturbation plots.

########################################################################################
########################################################################################
########################################################################################
"""


"""
Purpose: Collect perturbation sequences across trials
Given:
    region_idx: int mapping to the `example_idx` from a CWM-scanned region set
    perturb_seq_dataset: PerturbSeqDataset class object
    region_idx: region associated with perturb_seq_dataset
    trials: number of trials to generate sequences across
Output:
    perturb_seqs_dict_list = trials x -> ref/mutation -> seqlen x 4
"""

def collect_perturb_sequences(perturb_seq_dataset, region_idx, trials = 64):
    #Collect in trials
    perturb_seqs_dataset_list = [perturb_seq_dataset[region_idx] for trial in range(trials)]

    #Collect sequences across many trials
    perturb_seqs_dataset_list = [perturb_seq_dataset[region_idx] for trial in range(trials)]
    perturb_seqs_dict_list = [perturb_seqs_dataset_list[trial][0] for trial in range(trials)]

    return perturb_seqs_dict_list

"""
Purpose: Generate perturb predictions over many trials using bpnet.predict
Given:
    region_idx: int mapping to the `example_idx` from a CWM-scanned region set
    perturb_seqs: list of perturbation sequences: trials x -> ref/mutation -> seqlen x 4 (comes from collect_perturb_sequences)
    model: BPNetSeqModel object
Output:
    summary_dict = nested dict of task -> motif -> ref/mutation -> summary value
"""

def predict_perturbations(region_idx, perturb_seqs_dict_list, model, tasks):

    trials = len(perturb_seqs_dict_list)

    #Predict the signals across many trials (here, this is profile * exp(counts))
    #trials x task -> muts x seqlen x 2
    perturb_preds_list = [model.predict(np.array(list(perturb_seqs_dict_list[trial].values()))) for trial in range(trials)]

    #Collect average results from trials
    perturb_preds_array = {task: np.mean(np.array([perturb_preds_list[trial][task]
                                               for trial in range(trials)]), axis = 0)
                           for task in tasks}
    return perturb_preds_array

"""
Purpose: Collect perturbations from `predict_perturbations` and assign nested dict to label mutations.
Comments: Function goes from task -> muts x seqlen x 2     ==>     task -> ref/mutation -> seqlen x strands
Given:
    perturb_combinations: list of combinations of perturbations. MUST match the perturb_preds_array axis 0 of each task dict array.
    perturb_preds_array: dict with format: task -> combinations x seqlen x strand
Output:
    signals_dict = nested dict of task -> ref/mutation -> seqlen x strands
"""
def collect_perturb_preds(perturb_combinations, perturb_preds_array):

    #Store signals within relevant windows: task -> regions -> mutations -> window_len x strands
    signal_dict = {}
    for task in perturb_preds_array.keys():
        signal_dict[task]={}
        for mut_idx, mut_id in enumerate(perturb_combinations):
            signal = perturb_preds_array[task][mut_idx,:,:] #window_len x strand
            signal_dict[task][mut_id] = signal

    return signal_dict

"""
Purpose: For a signals_dict object summarize profile predictions within a window (sum, mean, max) and optionally apply pseudocounts.
Note: For max, will use the "Reference" profile as the baseline idx by which to measure the max values.
Given:
    signals_dict: nested dict of task -> ref/mutation -> seqlen x strands (comes from collect_perturb_preds)
    perturb_insts_df: pd.df of motif instances from PerturbComboSeqDataset's insts output
    summary_window: window width considered when computing summary statistics at the center of each motif
    summary_function: allows 'mean', 'sum', 'max'. If 'max', then will use reference max as anchored point.
    pseudo_count_quantile: If not None, compute quantile value from the ENTIRE reference profile, then add counts to both reference.
Output:
    summary_dict = nested dict of task -> motif -> ref/mutation -> summary value
"""
def collect_perturb_summary_preds(signal_dict, perturb_insts_df, summary_window = 70, summary_function = 'sum', pseudo_count_quantile = None):
    pc = 0

    #Get summary scores from the signal_dict
    summary_dict = {}
    for task in signal_dict.keys():
        summary_dict[task]={}

        #Add pseudo counts (pc) to all the signals
        if (pseudo_count_quantile is not None) and (pseudo_count_quantile != 0) :
            pc = np.percentile(signal_dict[task]['Reference'], pseudo_count_quantile * 100)

        for inst_idx, inst in perturb_insts_df.iterrows():
            motif=inst['pattern_name_unique']
            summary_dict[task][motif]={}

            #Collect boundaries of current motif
            centerpoint = int(inst['pattern_start'] + np.floor((inst['pattern_end'] - inst['pattern_start'])/2))
            lower_bounds = int(centerpoint - np.floor(summary_window/2))
            upper_bounds = int(centerpoint + np.floor(summary_window/2))

            #Call ref sig across motif slice within the summary_window and find max on pos and neg strand
            ref_sig_bounded = signal_dict[task]['Reference'][lower_bounds:upper_bounds,:]
            max_idx = np.argmax(ref_sig_bounded, axis = 0)

            for mut_id in signal_dict[task].keys():
                #Collect values across regions of interest.
                sig_bounded = signal_dict[task][mut_id][lower_bounds:upper_bounds,:] #window_len x strand

                if summary_function=='mean': summary_signal = np.mean(np.abs(sig_bounded)) + pc
                elif summary_function=='sum': summary_signal = np.sum(np.abs(sig_bounded)) + pc

                #Note: Find max value on pos and neg strand and add together
                elif summary_function == 'max': summary_signal = sig_bounded[max_idx,[0, 1]].sum() + pc
                else: summary_signal="Invalid summary function."
                #Summarize based on parameters
                summary_dict[task][motif][mut_id] = summary_signal
    return summary_dict

def random_seq_onehot(l):
    """Generate random sequence one-hot-encoded
    Args:
      l: sequence length
    """
    from concise.preprocessing import encodeDNA
    return encodeDNA([''.join(random.choices("ACGT", k=int(l)))])[0]

def decode_onehot_seq(onehot_seq):
    s = pd.DataFrame(onehot_seq)
    s.columns = ['A','C','G','T']
    seq_str = s.idxmax(1).str.cat()
    return seq_str

class PerturbComboSeqDataset(Dataset):
    """
    Purpose: Generate all combinatorial alternate sequences for a single region.
    Given:
        + dfi = (CWM-scanned pd.df of motifs) [requires `example_idx`, `pattern_name`, 'pattern_len`]
        + seqs = reference sequences in format [example_idx x seqlen x 4]
        + comb_max = max number of motifs to consider in combinations
        + comb_min = min number of motifs to consider in combinations
    Output:
        + perturb_seqs_dict: dictionary
    """

    def __init__(self, dfi, seqs, comb_max = None, comb_min = None):
#         assert len(dfi.example_idx.unique())==1, 'Only 1 region per class call.'
        self.dfi = dfi
        self.seqs = seqs
        self.comb_max = comb_max
        self.comb_min = comb_min

    def __len__(self):
        return len(self.dfi)

    def __getitem__(self, region_idx): #region_idx is equivalent to `example_idx`
        #Note: If there are empty regions that do not correspond to `dfi`, then it will return an dempty dict.

        insts = self.dfi[self.dfi['example_idx']==region_idx]
        assert all(insts.example_idx == region_idx), 'Not parsing `region_idx` correctly.'

        #Get reference sequence
        ref_seq = self.seqs[region_idx] #1000 x 4

        #Label unique motifs
        insts['pattern_name_unique'] = insts['pattern_name'] + '-' + (insts.groupby(['pattern_name']).cumcount()+1).astype(str)

        #Mark the max depth of combinations
        combination_min_depth = 1
        combination_max_depth = insts.shape[0]
        if(self.comb_max is not None): combination_max_depth = self.comb_max
        if(self.comb_min is not None): combination_min_depth = self.comb_min

        #Determine all combinations of unique motifs
        mut_combos = list(itertools.chain(*[list(itertools.combinations(insts['pattern_name_unique'], d))
                                                  for d in range(combination_min_depth,combination_max_depth+1)]))

        #Record mutant sequences used in this combinatorial trial: motif -> motif_len x 4
        mut_seqs_dict = {insts['pattern_name_unique'].iloc[idx]:
                          random_seq_onehot(insts['pattern_len'].iloc[idx]) for idx in range(len(insts))}
#         insts['mutation_seqs'] = [decode_onehot_seq(mut_seq) for mut_seq in mut_seqs_dict.values()]

        #Insert relevant mutant sequences in each combination
        def get_alt_seq_from_combo(combo, ref_seq, insts, mut_seqs_dict):
            alt_seq = ref_seq.copy()
            for motif in combo:
                inst = insts[insts['pattern_name_unique']==motif]
                alt_seq[inst.pattern_start.iloc[0]:inst.pattern_end.iloc[0]] = mut_seqs_dict[motif]
            return alt_seq

        # For each combination, generate: combo -> 1000 x 4
        alt_seqs_dict = {'_'.join(list(combo)):
            get_alt_seq_from_combo(combo = combo, ref_seq = ref_seq, insts = insts, mut_seqs_dict = mut_seqs_dict)
        for combo in mut_combos}

        perturb_seqs_dict = {**{'Reference': ref_seq}, **alt_seqs_dict}

        return perturb_seqs_dict, insts

"""
Purpose: Get motif coordinates as pd.df with a unique mutation name to loop through.
Inputs: region_idx = index of region to plot, dfi = pd.df of cwm-scanned regions
Outputs: mutant_coords_df = tidy pd.df of the mutant coordinates and their unique ids assigned
"""

def collect_mutant_coords(region_idx, dfi):
    #Collect motifs for
    mutant_coords_df = dfi[dfi.example_idx == region_idx]
    mutant_coords_df['mutation_idx'] = mutant_coords_df.row_idx

    #Get dictionary of motif occurrence counts
    patterns = list(mutant_coords_df['pattern_name'])
    pattern_dict = {a:list(range(1, b+1)) if b>1 else '' for a,b in Counter(patterns).items()}

    #Create list from dictionary with duplicated motifs being uniquely named
    mutant_coords_df['mutation_name'] = [i + '_mut' + str(pattern_dict[i].pop(0)) if len(pattern_dict[i]) else i + '_mut' for i in patterns]
    return(mutant_coords_df)

from bpnet.cli.contrib import ContribFile
from bpnet.BPNet import BPNetSeqModel
from bpnet.data import NumpyDataset #creates a set of nested arrays
from kipoi.data import Dataset

"""
Purpose: Given a contribution file, get all predictions from this file.

Required inputs:
model_dir = directory containing BPNet model
contrib_filepath = directory containing contrib.h5 from `bpnet contrib` on your regions of interest
tasks = tasks model was trained on

Outputs:
ref_npd = NumpyDataset of actual predictions (should be same number of regions as the contrib file)
ref_seqs = [region x seqlen x 4] array of 1-hot encoded sequences
region_coords = Coordinates by which contrib was run
"""

def collect_ref_preds(model_dir, contrib_filepath, tasks):

    print('Loading contributions...')

    #load contrib file and the coordinate ranges
    region_contribs = ContribFile(contrib_filepath)
    region_coords = region_contribs.get_ranges()

    #load the model
    bpnet_model = BPNetSeqModel.from_mdir(model_dir)

    print('Loading reference predictions...')
    # get 1-hot encoded sequence. Contrib idx == `example_idx` from dfi
    # will return a [number of regions originally considered x seq_length x 4]
    ref_seqs = region_contribs.get_seq()

    # Predict actual ChIP-nexus values
    ref_preds = bpnet_model.predict(ref_seqs)

    # Predict actual contribution scores
    ref_contrib = bpnet_model.contrib_score_all(ref_seqs, method='deeplift', aggregate_strand=True)
    ref_contrib_dict = {k: v * ref_seqs for k, v in ref_contrib.items()}

    #For each task, record the reference predictions and contrib scores
    ref_npd = NumpyDataset({t: {
        "pred": ref_preds[t],
        "contrib": {
            "profile": ref_contrib_dict[f"{t}/profile"],
            "count": ref_contrib_dict[f"{t}/count"],
        }} for t in tasks}, attrs={'index': 'example_idx'})

    return ref_npd, ref_seqs, region_coords


"""
Purpose: Given a region, plot BPNet profile predictions of each mutation event.
Required inputs:
contrib_filepath = directory containing contrib.h5 from `bpnet contrib` on your regions of interest
model_dir = directory containing BPNet model
dfi = motif regions (or regions of interest) following `cwm-scan` convention. This allows us to create mutations.
trials = number of times to generate random mutations that we will average predictions/contribs across
tasks = tasks model was trained on

Outputs:
alt_npd = NumpyDataset of mutated predictions  (should be same number of regions as the dfi)
"""

def collect_alt_preds(contrib_filepath, model_dir, dfi, tasks, trials = 24):

    print('Loading contributions...')

    #load contrib file and the coordinate ranges
    region_contribs = ContribFile(contrib_filepath)

    #load the model
    bpnet_model = BPNetSeqModel.from_mdir(model_dir)

    print('Loading reference predictions...')
    # get 1-hot encoded sequence. Contrib idx == `example_idx` from dfi
    # will return a [number of regions originally considered x seq_length x 4]
    ref_seqs = region_contribs.get_seq()

    #Get numpy ndarray with dimensions [trials x dfi_len x seq_width x 4]
    alt_seqs = np.array([PerturbSeqDataset(dfi = dfi, seqs = ref_seqs).load_all() for trial in range(trials)])

    #Structure is: trials -> task -> regions x seqlength x strands
    alt_preds = np.array([bpnet_model.predict(seq = alt_seqs[trial, :, :, :]) for trial in range(trials)])

    #Structure is: contrib_type -> regions x seqlength x 4
    alt_contribs = [bpnet_model.contrib_score_all(seq = alt_seqs[trial, :, :, :], method = 'deeplift', aggregate_strand = True) for trial in range(trials)]

    #Embedded function to loop across tasks
    def collect_trial_preds_and_contribs_across_task(alt_preds, alt_contribs, task):

        preds_across_trials = [alt_preds[trial][task] for trial in range(trials)]
        pred_mean_across_trials = sum(preds_across_trials)/trials #region x seqlen x 2

        profile_contrib_across_trials = [alt_contribs[trial][f'{task}/profile'] for trial in range(trials)]
        profile_contrib_mean_across_trials = sum(profile_contrib_across_trials)/trials #region x seqlen x 4

        counts_contrib_across_trials = [alt_contribs[trial][f'{task}/count'] for trial in range(trials)]
        counts_contrib_mean_across_trials = sum(counts_contrib_across_trials)/trials #region x seqlen x 4

        summary_dict = {"pred": pred_mean_across_trials,
                            "contrib": {
                                "profile": profile_contrib_mean_across_trials,
                                "count": counts_contrib_mean_across_trials}}
        return summary_dict

    #Collect data across trials
    alt_npd = NumpyDataset({t: collect_trial_preds_and_contribs_across_task(alt_preds = alt_preds,
                                                            alt_contribs = alt_contribs,
                                                            task = t) for t in tasks}, attrs={'index': 'example_idx'})
    return alt_npd


"""
NOTE: THIS FUNCTION IS DEPRECATED. USE collect_ref_preds, and collect_alt_preds instead.


Purpose: Given a region, plot BPNet profile predictions of each mutation event.
Required inputs:
model_dir = directory containing BPNet model
regions_df = regions of interest (only use this to assign the `name` column to keep track of regions)
contrib_filepath = directory containing contrib.h5 from `bpnet contrib` on your regions of interest

Outputs:
ref_npd = NumpyDataset of actual predictions (should be same number of regions as the contrib file)
alt_npd = NumpyDataset of mutated predictions  (should be same number of regions as the dfi)
region_coords = Coordinates by which contrib was run
"""

def collect_preds_and_contribs(model_dir, contrib_filepath, tasks, dfi, regions_df):

    print('Loading contributions...')
    #load contrib file and the coordinate ranges
    region_contribs = ContribFile(contrib_filepath)
    region_coords = region_contribs.get_ranges()
    region_coords['region_name'] = regions_df.name

    #load the model
    bpnet_model = BPNetSeqModel.from_mdir(model_dir)

    print('Loading reference predictions...')
    # get 1-hot encoded sequence. Contrib idx == `example_idx` from dfi
    # will return a [number of regions originally considered x seq_length x 4]
    ref_seqs = region_contribs.get_seq()

    # Predict actual ChIP-nexus values
    ref_preds = bpnet_model.predict(ref_seqs)

    # Predict actual contribution scores
    ref_contrib = bpnet_model.contrib_score_all(ref_seqs, method='deeplift', aggregate_strand=True)
    ref_contrib_dict = {k: v * ref_seqs for k, v in ref_contrib.items()}

    #For each task, record the reference predictions and contrib scores
    ref_npd = NumpyDataset({t: {
        "pred": ref_preds[t],
        "contrib": {
            "profile": ref_contrib_dict[f"{t}/profile"],
            "count": ref_contrib_dict[f"{t}/count"],
        }} for t in tasks}, attrs={'index': 'example_idx'})

    print('Loading alternate predictions...')

    #Generate alternate sequence for each motif
    # will return a [number of motifs found x seq_length x 4]
    alt_seqs = PerturbSeqDataset(dfi, ref_seqs).load_all()
    alt_preds = bpnet_model.predict(seq = alt_seqs)

    alt_contrib = bpnet_model.contrib_score_all(seq = alt_seqs, method = 'deeplift', aggregate_strand = True)
    alt_contrib_dict = {k: v * alt_seqs for k, v in alt_contrib.items()}

    #For each task, record the mutated predictions and contrib scores
    alt_npd = NumpyDataset({t: {
        "pred": alt_preds[t],
        "contrib": {
            "profile": alt_contrib_dict[f"{t}/profile"],
            "count": alt_contrib_dict[f"{t}/count"],
        }} for t in tasks}, attrs={'index': 'example_idx'})

    return ref_npd, alt_npd, region_coords

"""
Purpose: For a region, parse NumpyDatasets to generate tidy pd.dfs of preds
Inputs:
region_idx = index of region to plot
ref_npd = NumpyDataset of actual predictions (should be same number of regions as the contrib file)
alt_npd = NumpyDataset of mutated predictions  (should be same number of regions as the dfi)
dfi = pd.df of cwm-scanned regions

Outputs:
preds_df = tidy pd.df of the prediction profiles
mutant_coords_df = tidy pd.df of the mutant coordinates and their unique ids assigned
"""

def collect_dfs_across_ref_and_mutants(region_idx, ref_npd, alt_npd, dfi):

    #Collect mutant coordinates with unique names
    mutant_coords_df = collect_mutant_coords(region_idx, dfi)

    # Process reference predictions
    ref_preds_tidy_df = pred_npd_to_tidy_df(npd = ref_npd, npd_idx = region_idx) # Obtain reference profile pd.df predictions
    ref_preds_tidy_df['mutation_idx'] = 'Reference' # based on row_idx from the `dfi`, hard set.
    ref_preds_tidy_df['mutation_name'] = 'Reference'

    # Process alternate predictions (allow for multiple same motif mutations)
    mutation_preds_tidy_df = pd.DataFrame()
    for idx in list(dfi.row_idx[dfi.example_idx == region_idx]): #from each `row_idx` associated with `example_idx`
        alt_preds_tidy_df = pred_npd_to_tidy_df(npd = alt_npd, npd_idx = idx)
        alt_preds_tidy_df['mutation_idx'] = idx
        alt_preds_tidy_df['mutation_name'] = mutant_coords_df['mutation_name'][mutant_coords_df['row_idx']==idx].iloc[0]
        mutation_preds_tidy_df = mutation_preds_tidy_df.append(alt_preds_tidy_df)

    # Combine mutations and predictions
    preds_df = ref_preds_tidy_df.append(mutation_preds_tidy_df)

    return(preds_df, mutant_coords_df)


"""
Purpose: For a region, plot the cwm-scanned motifs
Inputs:
    preds_perturb_df = tidy pd.df of the perturbation predictions (REQUIRES 'plot_position', 'preds_profile', 'strand', 'mutation_idx' and 'task' columns)
    ax = plot object to add to
    task = task to filter `preds_perturb_df` by
    mutation_idx = mutation to filter `preds_perturb_df` by
    mutant_coords_df = tidy pd.df of the cwm-scanned coordinates with unique `mutation_idx` and `mutation_name` columns and 'plot_start', 'plot_end'
    region_coord = Coordinate by which contrib was run (has `region_name` attribute)

Outputs:
    ax = plot object to save as
"""
def plot_nexus_prediction_perturb(preds_perturb_df, ax, task, mutation_idx, mutant_coords_df, region_coord, xlimits = None, ylimits = None):

    #Get motifs without mutant specification
    motif_coords_df = mutant_coords_df[['pattern_start', 'pattern_end',
                                                            'pattern_start_abs', 'pattern_end_abs',
                                                            'row_idx','pattern_name', 'example_chrom',
                                                            'plot_start', 'plot_end']]

    # Subset regions to match prediction
    pred_subset_df = preds_perturb_df[(preds_perturb_df['task'] == task) & (preds_perturb_df['mutation_idx']==mutation_idx)]
    mutant_coord_df = mutant_coords_df[(mutant_coords_df['mutation_idx']==mutation_idx)]

    #Format mutation names for the subplot labels
    if mutation_idx == 'Reference':
        label_name = f'ChIP: {task}, Condition: Reference'
    else:
        mutant_name = ''.join(map(str, mutant_coords_df['mutation_name'][(mutant_coords_df['mutation_idx']==mutation_idx)]))
        label_name = f'ChIP: {task}, Condition: {mutant_name}'

    #Plot motif boxes
    plot_pred_motif_boxes(motif_coords_df = motif_coords_df, ax = ax)

    #Plot mutant box
    if mutant_coord_df.shape[0] > 0:
        plot_pred_motif_boxes(motif_coords_df = mutant_coord_df, ax = ax, color = '#B20000')

    # Plot profile predictions
    plot_nexus_prediction(preds_df = pred_subset_df, ax = ax, xlimits = xlimits, ylimits = ylimits)

    #Label aesthetics
    ax.set_title(label = label_name, loc = 'left', fontsize='x-large')
    ax.set_xlabel(f'Genomic position (bp, across {region_coord.chrom})')
    ax.ticklabel_format(axis = 'both', style = 'plain', useOffset = False)
    return(ax)


"""
Purpose: Given a region, plot the single-mutation predictions and contributions
Inputs:
    region_idx: `example_idx` of desired region
    ref_npd: NumpyDataset of actual predictions (should be same number of regions as the contrib file)
    alt_npd: NumpyDataset of mutated predictions  (should be same number of regions as the dfi)
    region_coord: from `contrib.get_ranges().iloc[region_idx]`
    dfi: pd.df of cwm-scanned regions containing desired region
    tasks: list of tasks to plot and find regions across
    output_prefix: prefix of output .pdf and .png
    xlimits: within-window coordinates to plot across
    plot_genomic_coords: bool deciding whether to plot across genomic coordinates. If True, xlimits will be automatically updated.
"""

def plot_perturb_single_mutation_metaplot(region_idx, ref_npd, alt_npd,
                                          region_coord, dfi, tasks, output_prefix,
                                          xlimits = None, plot_genomic_coords = True):

    # Get predictions and mutant scores
    preds_df, mutant_coords_df = collect_dfs_across_ref_and_mutants(region_idx = region_idx, ref_npd = ref_npd, alt_npd = alt_npd, dfi = dfi)

    # Set reference to be first column
    condition_order = list(preds_df['mutation_name'].unique())
    mutant_coords_df['mutation_name_cat'] = pd.Categorical(mutant_coords_df['mutation_name'], categories = condition_order)
    preds_df['mutation_name_cat'] = pd.Categorical(preds_df['mutation_name'], categories = condition_order)

    # Add in genomic position of data
    preds_df['genomic_position']=preds_df['position'] + region_coord['start']

    #Choose whether to plot genomic coordinates or window coordinates
    if plot_genomic_coords:
        preds_df['plot_position']=preds_df['genomic_position']
        mutant_coords_df['plot_start'] = mutant_coords_df['pattern_start_abs']
        mutant_coords_df['plot_end'] = mutant_coords_df['pattern_end_abs']
        if xlimits is not None: xlimits = xlimits + region_coord['start']
    else:
        preds_df['plot_position']=preds_df['position']
        mutant_coords_df['plot_start'] = mutant_coords_df['pattern_start']
        mutant_coords_df['plot_end'] = mutant_coords_df['pattern_end']


    #Define figure parameters
    num_cols = len(tasks)
    num_conds = len(dfi[dfi['example_idx']==region_idx]) + 1
    num_rows = (num_conds) *2 #Preds and contrib
    gs_dict={'width_ratios': list(itertools.repeat(1, num_cols)),
             'height_ratios': list(np.array(list(itertools.repeat([2,1], num_conds))).flatten())
            } #2:1 ratio of pred to contrib

    mutations_idx = np.append(['Reference'], mutant_coords_df['mutation_idx'])

    fig, axs = plt.subplots(nrows = num_rows, ncols = num_cols,
                            figsize=(num_cols*8, num_conds*5),
                            gridspec_kw=gs_dict, constrained_layout=True) #sharex=True, figsize = (width, height)
    fig.suptitle(f'Predictions of mutations across {region_coord.region_name}', fontsize='xx-large', fontweight = 'bold')

    for task_idx in range(len(tasks)):
        task = tasks[task_idx]

        #Set the same ylimits across task
        ymax_pred_across_task = myceiling(preds_df['preds_profile'][(preds_df['task'] == task)].abs().max(), base = 0.1)
        ylimits_pred_across_task = [-ymax_pred_across_task, ymax_pred_across_task]
        ylimits_contrib_across_task = [myfloor((np.array(contrib_npd_to_df(npd = ref_npd, task = task, region_idx = region_idx)).flatten().min()), base = .01),
                                       myceiling((np.array(contrib_npd_to_df(npd = ref_npd, task = task, region_idx = region_idx)).flatten().max()), base = .01)]

        for contrib_idx in range(len(mutations_idx)):
            mut_idx = mutations_idx[contrib_idx]

            if mutations_idx[contrib_idx]=='Reference':
                #Plot prediction
                mutant_coord_df = mutant_coords_df[mutant_coords_df['mutation_idx']==mutations_idx[contrib_idx]]
                plot_nexus_prediction_perturb(preds_perturb_df = preds_df, ax = axs[(contrib_idx*2),task_idx], task = task,
                                              mutation_idx = mut_idx, mutant_coords_df = mutant_coords_df,
                                              region_coord = region_coord, xlimits = xlimits, ylimits = ylimits_pred_across_task)

                #Plot contribution
                contrib_ref_df = contrib_npd_to_df(npd = ref_npd, task = task, region_idx = region_idx)
                if plot_genomic_coords: contrib_ref_df = contrib_ref_df.set_index(contrib_ref_df.index.tolist() + region_coord['start'])

                ref_logo = plot_contrib_track(contrib_ref_df, ax = axs[(contrib_idx*2)+1, task_idx],
                                              xlimits = xlimits, ylimits = ylimits_contrib_across_task)
                for mut_box_idx in mutant_coords_df['mutation_idx']:
                    plot_contrib_motif_box(ref_logo, motif_df = mutant_coords_df[mutant_coords_df['mutation_idx']==mut_box_idx])

            else:
                mutant_coord_df = mutant_coords_df[mutant_coords_df['mutation_idx']==int(mutations_idx[contrib_idx])]
                #Plot prediction
                plot_nexus_prediction_perturb(preds_perturb_df = preds_df, ax = axs[(contrib_idx*2),task_idx], task = task,
                                              mutation_idx = int(mut_idx), mutant_coords_df = mutant_coords_df,
                                              region_coord = region_coord, xlimits = xlimits, ylimits = ylimits_pred_across_task)

                #Plot contribution
                contrib_alt_df = contrib_npd_to_df(npd = alt_npd, task = task, region_idx = int(mut_idx))
                if plot_genomic_coords: contrib_alt_df = contrib_alt_df.set_index(contrib_alt_df.index.tolist() + region_coord['start'])

                alt_logo = plot_contrib_track(contrib_alt_df, ax = axs[(contrib_idx*2)+1, task_idx],
                                              xlimits = xlimits, ylimits = ylimits_contrib_across_task)
                for mut_box_idx in mutant_coords_df['mutation_idx']:
                    if mut_box_idx == int(mut_idx):
                        plot_contrib_motif_box(alt_logo, motif_df = mutant_coords_df[mutant_coords_df['mutation_idx']==mut_box_idx])
                        plot_contrib_motif_box(alt_logo, motif_color ="red", motif_df = mutant_coords_df[mutant_coords_df['mutation_idx']==mut_box_idx])
                    else:
                        plot_contrib_motif_box(alt_logo, motif_df = mutant_coords_df[mutant_coords_df['mutation_idx']==mut_box_idx])

    fig.savefig(f'{output_prefix}.pdf')
    fig.savefig(f'{output_prefix}.png')

class PerturbSeqDataset(Dataset):
    #Kudos to Ziga
    #Given a CWM-scanned pd.df of motifs, generate an alternative sequence by randomizing sequences at the motif site.
    #Based on the `row_idx`, parse through the seqs

    def __init__(self, dfi, seqs):
        self.dfi = dfi
        self.seqs = seqs

    def __len__(self):
        return len(self.dfi)

    def __getitem__(self, idx):
        inst = self.dfi.iloc[idx]
        assert inst.row_idx == idx, 'Your row_idx does not match actual rows. Please reassign if using load_all.' #if returning error, does row_idx exist in your dfi? Need to add it manually.
        ref_seq = self.seqs[inst.example_idx]

        # generate the alternative sequence
        alt_seq = ref_seq.copy()
        alt_seq[int(inst.pattern_start):int(inst.pattern_end)] = random_seq_onehot(inst.pattern_end - inst.pattern_start)

        return alt_seq
