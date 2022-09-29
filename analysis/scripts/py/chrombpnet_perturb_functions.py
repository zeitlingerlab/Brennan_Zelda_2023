"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions to applying perturbations
"""

import os
import sys
import json
import pandas as pd
import numpy as np
import re

#Custom functions
sys.path.insert(0, f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/analysis/scripts/py')
from chrombpnet_predict_functions import predict_chrombpnet, one_hot_decode_sequence

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

def generate_injected_seq(primary_motif, primary_position,
                          secondary_motif='', secondary_distances=[], seqlen=1000,
                          weights = [.25, .25, .25, .25]):
    """
    Kudos: I straight up jacked this from Ziga Avsec's bpnet code bpnet.bpnet.simulate#L37
    Purpose: Given a sequence, inject a motif centered on a position.
        + If you want to inject 2 motifs to compare distances, use:
            + side_motif to define the second motif sequence
            + side_distances to define the second motif center
    """
    random_seq = generate_random_seq(seqlen = seqlen, weights = weights)
    injected_seq = insert_motif(seq = random_seq, motif = primary_motif, position = primary_position)
    if len(secondary_distances)>1:
        print('Warning! You have entered multiple side distances.',
              'This will inject multiple secondary motifs into the SAME sequence.',
              'If you do not want this, loop this function through single side distances.')
    for d in secondary_distances:
        injected_seq = insert_motif(injected_seq, secondary_motif, d)
    return injected_seq

def predict_injected_seq(model_chrombpnet, model_bias,
                         primary_motif, primary_position,
                         secondary_motif = '', secondary_positions = [],
                         input_seqs = None,
                         input_seqlen = 2114,
                         output_seqlen = 1000,
                         trials = 512, batch_size=32, verbose = False):
    """
    Purpose:
        + Given a primary motif A, inject that motif and predict.
        + Given a secondary motif B, inject that motif with the primary motif as AB and measure synergy.
    Inputs:
        + model_chrombpnet: ATAC_seq coverage keras model imported by `load_chrombpnet` already
        + model_bias: Tn5 bias keras model imported by `load_chrombpnet` already
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

    #predict profiles
    #print("Predicting...")
    _, preds_w_counts = predict_chrombpnet(model_chrombpnet = model_chrombpnet, model_bias = model_bias,
                                           seqs = seqs, batch_size=batch_size, no_bias = True)

    #print("Collecting results...")
    #separate different injections and average across trials, which are the smallest margin
    pred_dict = {}
    null_lower_range = 0 * trials
    null_upper_range = (0 + 1) * trials
    pred_dict['null'] = np.mean(preds_w_counts[null_lower_range:null_upper_range], axis = 0)

    A_lower_range = 1 * trials
    A_upper_range = (1 + 1) * trials
    pred_dict['A'] = np.mean(preds_w_counts[A_lower_range:A_upper_range], axis = 0)

    counter = 2
    labels = ['B', 'AB']
    distances_dict = {}
    for i,d in enumerate(secondary_positions):
        p_dict = {}
        for r in labels: #hard-coded because of the four states
            lower_range = counter * trials
            upper_range = (counter + 1) * trials
            p_dict[r]=np.mean(preds_w_counts[lower_range:upper_range], axis = 0)
            counter +=1
        distances_dict[d] = p_dict
    pred_dict['distances'] = distances_dict

    return pred_dict

def summarize_injected_seq(pred_dict, primary_motif, primary_position,
                           secondary_motif = '', secondary_positions = [],
                           measurement_window_around_motif = 50, pseudocount_thresh = .2):
    """
    Purpose:
        + Given a pred_dict from predict_injected_seq, summarize values
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

    windowA_min = primary_position - (measurement_window_around_motif // 2)
    windowA_max = primary_position + (measurement_window_around_motif // 2)

    results_dict = {}
    results_dict['pc_A'] = np.quantile(pred_dict['A'], pseudocount_thresh)
    results_dict['pc_null'] = np.quantile(pred_dict['null'], pseudocount_thresh)
    results_dict['all_sum_A'] = pred_dict['A'].sum()
    results_dict['all_sum_null'] = pred_dict['null'].sum()
    results_dict['windowA_sum_A'] = pred_dict['A'][windowA_min:windowA_max].sum()
    results_dict['windowA_sum_null'] = pred_dict['null'][windowA_min:windowA_max].sum()

    # all_motifA_max_pos = pred_dict['A'].argmax()
    windowA_motifA_max_pos = pred_dict['A'][windowA_min:windowA_max].argmax()
    results_dict['windowA_max_A'] = pred_dict['A'][windowA_min:windowA_max].max()
    results_dict['windowA_max_null'] = pred_dict['null'][windowA_min:windowA_max][windowA_motifA_max_pos]

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
        results_dict['windowA_max_AB'] = pred_dict['distances'][d]['AB'][windowA_min:windowA_max][windowA_motifA_max_pos]
        results_dict['windowA_max_B'] = pred_dict['distances'][d]['B'][windowA_min:windowA_max][windowA_motifA_max_pos]

        #Compute results across partial window at motifB
        results_across_dist_dict['windowB_sum_AB'] = pred_dict['distances'][d]['AB'][windowB_min:windowB_max].sum()
        results_across_dist_dict['windowB_sum_B'] = pred_dict['distances'][d]['B'][windowB_min:windowB_max].sum()
        results_across_dist_dict['windowB_sum_A'] = pred_dict['A'][windowB_min:windowB_max].sum()
        results_across_dist_dict['windowB_sum_null'] = pred_dict['null'][windowB_min:windowB_max].sum()
        windowB_motifB_max_pos = pred_dict['distances'][d]['B'][windowB_min:windowB_max].argmax()
        results_dict['windowB_max_AB'] = pred_dict['distances'][d]['AB'][windowB_min:windowB_max][windowB_motifB_max_pos]
        results_dict['windowB_max_A'] = pred_dict['distances'][d]['B'][windowB_min:windowB_max][windowB_motifB_max_pos]
        results_dict['windowB_max_B'] = pred_dict['A'][windowB_min:windowB_max].max()
        results_dict['windowB_max_null'] = pred_dict['null'][windowB_min:windowB_max][windowB_motifB_max_pos]

        #Record metadata
        results_across_dist_dict['primary_motif'] = primary_motif
        results_across_dist_dict['primary_position'] = primary_position
        results_across_dist_dict['secondary_motif'] = secondary_motif
        results_across_dist_dict['secondary_position'] = window_distance

        r_df = pd.DataFrame(results_across_dist_dict, index=[0])
        perturb_df = perturb_df.append(r_df)

    return(perturb_df)

def generate_alt_sequences(ref_seq, motifs, motif_unique_col, sequence_window_start_col, sequence_window_end_col,
                           comb_max = None, comb_min = None):
    """
    Purpose: Generate all combinations of mutated sequences based on given sequence and desired mutant coordinates.
    Given:
        ref_seq: np.array with shape [l x 4] of the one-hot encoded reference sequence
        motifs: pd.df of the motifs across the reference sequence
        motif_unique_col: column name designating the unique motif label
        sequence_window_start_col: column name designating where in the input sequence the motifs start
        sequence_window_end_col: column name designating where in the input sequence the motifs end
        comb_max: max number of simultaneous mutations allowed
        comb_min: min number of simultaneous mutations allowed
    Output:
        perturb_seqs_dict_list = trials x -> ref/mutation -> seqlen x 4
    """
    import sys
    sys.path.insert(0, f'/n/projects/mw2098/shared_code/chrombpnet/functions')
    from perturb import generate_random_seq
    from predict import one_hot_encode_sequence
    from itertools import combinations, chain


    #Insert relevant mutant sequences in each combination
    def get_alt_seq_from_combo(combo, ref_seq, motifs, motif_unique_col, sequence_window_start_col, sequence_window_end_col, mut_seqs_dict):
        alt_seq = ref_seq.copy()
        for mut in combo:
            motif = motifs[motifs[motif_unique_col]==mut]
            alt_seq[motif[sequence_window_start_col].iloc[0]:motif[sequence_window_end_col].iloc[0]] = mut_seqs_dict[mut]
        return alt_seq

    #Mark the max depth of combinations
    c_min_depth = 1
    c_max_depth = motifs.shape[0]
    if(comb_max is not None): c_max_depth = comb_max
    if(comb_min is not None): c_min_depth = comb_min

    #Ensure these coordinate positions are integers
    motifs[sequence_window_start_col] = motifs[sequence_window_start_col].astype(int)
    motifs[sequence_window_end_col] = motifs[sequence_window_end_col].astype(int)

    #Determine all combinations of unique motifs
    mut_combos = [list(combinations(motifs[motif_unique_col], d)) for d in range(c_min_depth, c_max_depth+1)]
    mut_combos = list(chain(*mut_combos))

    #Record mutant sequences used in this combination: motif -> motif_len x 4
    mut_seqs_dict = {motifs[motif_unique_col].iloc[idx]:
                          one_hot_encode_sequence(generate_random_seq(motifs[sequence_window_end_col].iloc[idx] - motifs[sequence_window_start_col].iloc[idx]))
                     for idx in range(len(motifs))}

    # Generate: combo -> 1000 x 4
    alt_seqs_dict = {'_'.join(list(combo)): get_alt_seq_from_combo(combo = combo,
                                                                   ref_seq = ref_seq,
                                                                   motifs = motifs,
                                                                   motif_unique_col = motif_unique_col,
                                                                   sequence_window_start_col = sequence_window_start_col,
                                                                   sequence_window_end_col = sequence_window_end_col,
                                                                   mut_seqs_dict = mut_seqs_dict)
                     for combo in mut_combos}

    #Add reference sequence
    perturb_seqs_dict = {**{'Reference': ref_seq}, **alt_seqs_dict}

    #Clean up results
    return(perturb_seqs_dict)

def perturb_single_enhancer(
    enhancer_chrom, enhancer_start, enhancer_end, enhancer_id,
    motifs_df_path, model_file, bias_file, fasta_file,
    input_seqlen = 2114,
    output_seqlen = 1000,
    region_index_column = 'region_id',
    motif_name_column = 'pattern_name',
    motif_start_column = 'start',
    motif_end_column = 'end',
    motif_chrom_column = 'chrom',
    comb_min = 1,
    comb_max = 1,
    trials = 4):
    """
    Purpose: Generate 2 pd_dfs (preds and motifs) of mutated motif combinations across a SINGLE genomic region.
    Given:
        + enhancer_chrom: string designating which chromosome the enhancer is across
        + enhancer_start: string designating which start the enhancer is across
        + enhancer_end: string designating which end the enhancer is across
        + enhancer_id: string designating which name the enhancer is
        + motifs_df_path: pd.df of the motif set you are using (0-based)
        + model_file: Path to the trained model (specified in `chrombpnet`
        + bias_file: Path to the trained bias model (specified in `chrombpnet`
        + fasta_file: Path to the .fasta file to extract sequences from.
        + input_seqlen: Input sequence length of model.
        + output_seqlen: Output profile length of model.
        + comb_max: Upper limit of simultaneous mutations to conduct across each region.
        + comb_min: Lower limit of simultaneous mutations to conduct across each region.
        + region_index_column: Column name that motifs_df and motif_groups_df have in common to keep region combinations consistent.
        + motif_start_column: Column name that motifs_df contains to denote GENOMIC coordinates of motifs.
        + motif_end_column: Column name that motifs_df contains to denote GENOMIC coordinates of motifs.
        + motif_name_column: Column name that motifs_df contains to denote pattern/motif names for unique labeling.
        + motif_chrom_column: Column name that motifs_df contains to denote chromosome names.
        + trials: How many motif mutations should each combination average across?
    Output:
        + preds_df: tidy pd.df of profile predictions
        + motifs_df: coordinates and motif information for motifs that belong to the enhancer
    """

    #Packages
    sys.path.insert(0, f'/n/projects/mw2098/shared_code/basepairmodels/functions')
    from motifs import extract_seqs_from_df, resize_coordinates

    sys.path.insert(0, f'/n/projects/mw2098/shared_code/chrombpnet/functions')
    from predict import one_hot_encode_sequence, load_chrombpnet, predict_chrombpnet
    from perturb import generate_alt_sequences


    #Extract sequence from enhancer of interest
    regions_df = pd.DataFrame([enhancer_chrom, enhancer_start, enhancer_end, '*']).transpose()
    regions_df.columns =['chrom', 'start', 'end', 'strand']
    regions_input_df = resize_coordinates(regions_df, input_seqlen, 'center')
    regions_output_df = resize_coordinates(regions_df, output_seqlen, 'center')
    seqs = extract_seqs_from_df(coords_df = regions_input_df, fasta_path = fasta_file, chrom_column = 'chrom', start_column = 'start', end_column = 'end')

    assert len(seqs)==1, 'This code is designed for one region only.'

    seqs_1he = one_hot_encode_sequence(seqs[0])

    #Find set of mapped motifs across this region
    motifs_df = pd.read_csv(motifs_df_path, sep = '\t')
    motifs_df = motifs_df.sort_values([region_index_column, motif_name_column])
    motifs_df['pattern_name_unique'] = motifs_df[motif_name_column] + '-' + (motifs_df.groupby([region_index_column, motif_name_column]).cumcount()+1).astype(str)
    motif_set_df = motifs_df[(motifs_df[motif_chrom_column].values == regions_df.chrom.values) & \
                            (motifs_df[motif_start_column].values >= regions_df.start.values) &
                            (motifs_df[motif_end_column].values <= regions_df.end.values)]

    #Mark the coordinates of the motifs relative to the windows
    motif_set_df['motif_window_start_input'] = motif_set_df[motif_start_column] - regions_input_df.start.values
    motif_set_df['motif_window_end_input'] = motif_set_df[motif_end_column] - regions_input_df.start.values
    motif_set_df['motif_window_start_output'] = motif_set_df[motif_start_column] - regions_output_df.start.values
    motif_set_df['motif_window_end_output'] = motif_set_df[motif_end_column] - regions_output_df.start.values

    #Load model of interest
    # print('Loading model...')
    chrombpnet_model, bias_model = load_chrombpnet(model_file, bias_file)

    preds_by_trial = []
    #For each trial
    for trial in range(trials):
        s = []
        pis_df = pd.DataFrame(columns = ['mut', region_index_column])
        #Generate all mutated sequences,
        perturb_seq = generate_alt_sequences(ref_seq = seqs_1he, motifs = motif_set_df, motif_unique_col = 'pattern_name_unique',
                                  sequence_window_start_col = 'motif_window_start_input',
                                  sequence_window_end_col = 'motif_window_end_input',
                                  comb_max = comb_max, comb_min = comb_min)

        perturb_seq_array = np.array([v for v in perturb_seq.values()])
        s.append(perturb_seq_array)

        #Record the indexes
        pi = pd.DataFrame([k for k in perturb_seq.keys()], columns = ['mut'])
        pi[region_index_column] = enhancer_id
        pis_df = pis_df.append(pi)
        pis_df['perturb_idx'] = range(pis_df.shape[0])
        s = np.vstack(s)

        #Predict the profile
        _, preds_w_counts = predict_chrombpnet(model_chrombpnet = chrombpnet_model, model_bias = bias_model, seqs = s, no_bias = True)
        preds_by_trial.append(preds_w_counts) #list[region x window x channel]

    #Average trials together: task -> region x window x strand where region = each perturb_idx of pis_df
    # print('Collecting perturbations together ...')
    perturb_preds = {'atac': np.mean(np.array([preds_by_trial[trial] for trial in range(trials)]), axis = 0)}

    #Save predictions as a tidy.df for plotting.
    preds_df = pd.DataFrame(perturb_preds['atac'], columns = [str(i) for i in list(range(output_seqlen))])
    preds_df[pis_df.columns] = pis_df[pis_df.columns]
    preds_df = preds_df.melt(id_vars = list(pis_df.columns), var_name = 'position', value_name = 'pred')
    preds_df.position = preds_df.position.astype(int)

    #Return both predictions and motif information
    return preds_df, motif_set_df
