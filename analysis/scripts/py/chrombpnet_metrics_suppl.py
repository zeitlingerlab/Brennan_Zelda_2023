"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions to computing metrics not covered by `chrombpnet.CLI`
"""

import os
import json
import pandas as pd
import numpy as np
import logging

##################################################################
#Compute auPRC
##################################################################

def compute_auprc(regions_path, profileA, profileB,
                  valid_chroms = None, test_chroms = None,
                  pos_min_threshold=0.05,
                  neg_max_threshold=0.01,
                  required_min_pos_counts=2.5,
                  binsizes=[1, 2, 4, 10, 20, 50],
                  seqlenth = 1000, exclude_zero_profiles = False):
    """
    Kudos to Zahoor, code modified from his basepairmodels.basepairmodels.cli.metrics.py
    Description:
        + collect count information across a set of regions,
        + these need to be exponentiated, not logits.
    Inputs: 2 bigwigs and a set of .narrowPeak regions
    Outputs: data.frame of collected counts and 2 correlation scalars
    """
    import pyBigWig
    from tqdm import tqdm
    from scipy.ndimage import gaussian_filter1d
    from scipy.spatial.distance import jensenshannon
    from scipy.stats import pearsonr, spearmanr, multinomial
    from scipy.special import logsumexp

    # Read in information
    peaks_df = pd.read_csv(regions_path, sep='\t', header=None,
                           names=['chrom', 'st', 'end', 'name', 'score',
                                  'strand', 'signal', 'p', 'q', 'summit'])

    # create new column for peak pos
    peaks_df['summit_pos'] = peaks_df['st'] + peaks_df['summit']

    # create new column for start pos
    peaks_df['start_pos'] = peaks_df['summit_pos'] - seqlenth // 2

    # create new column for end pos
    peaks_df['end_pos'] = peaks_df['summit_pos'] + seqlenth // 2

    # select only the chrom & summit positon columns
    allPositions = peaks_df[['chrom', 'start_pos', 'end_pos']]
    allPositions = allPositions.reset_index(drop=True)

    # open the two bigWig files
    try:
        bigWigProfileA = pyBigWig.open(profileA)
        bigWigProfileB = pyBigWig.open(profileB)

    except Exception as e:
        logging.error("Problems occurred when opening one of the input files: "
                      "{}".format(str(e)))

    coverageA = []
    coverageB = []
    for idx, row in tqdm(allPositions.iterrows(), total=allPositions.shape[0]):

        #Collect region information
        chrom = row['chrom']
        start = row['start_pos']
        end = row['end_pos']

        #Extract bigwig profiles
        try:
            profileA = np.nan_to_num(np.array(
                bigWigProfileA.values(chrom, start, end)))
            profileB = np.nan_to_num(np.array(
                bigWigProfileB.values(chrom, start, end)))
        except Exception as e:
            raise NoTracebackException(
                "Error retrieving values {}, {}, {}".format(chrom, start, end))

        coverageA.append(profileA)
        coverageB.append(profileB)

    # Extract information into an array
    matA = np.stack(coverageA)
    matB = np.stack(coverageB)

    # Compute auPRC metrics
    if valid_chroms or test_chroms:
        valid_idxs = allPositions.index[allPositions['chrom'].isin(valid_chroms)].tolist()
        test_idxs = allPositions.index[allPositions['chrom'].isin(test_chroms)].tolist()
        both_idxs = allPositions.index[allPositions['chrom'].isin(valid_chroms + test_chroms)].tolist()
        train_idxs = np.ones(allPositions.shape[0], dtype=bool)
        train_idxs[both_idxs] = False

        test_auprc_df = eval_profile(yt = matA[test_idxs], yp = matB[test_idxs],
                                pos_min_threshold = pos_min_threshold, neg_max_threshold = neg_max_threshold,
                                required_min_pos_counts = required_min_pos_counts, binsizes = binsizes)
        test_auprc_df['chrom_cat'] = 'test'
        valid_auprc_df = eval_profile(yt = matA[valid_idxs], yp = matB[valid_idxs],
                                pos_min_threshold = pos_min_threshold, neg_max_threshold = neg_max_threshold,
                                required_min_pos_counts = required_min_pos_counts, binsizes = binsizes)
        valid_auprc_df['chrom_cat'] = 'valid'
        train_auprc_df = eval_profile(yt = matA[train_idxs], yp = matB[train_idxs],
                                pos_min_threshold = pos_min_threshold, neg_max_threshold = neg_max_threshold,
                                required_min_pos_counts = required_min_pos_counts, binsizes = binsizes)
        train_auprc_df['chrom_cat'] = 'train'

        auprc_df = pd.concat([test_auprc_df, valid_auprc_df, train_auprc_df])

    else:
        auprc_df = eval_profile(yt = matA, yp = matB,
                                pos_min_threshold = pos_min_threshold, neg_max_threshold = neg_max_threshold,
                                required_min_pos_counts = required_min_pos_counts, binsizes = binsizes)
    return(auprc_df)

def eval_profile(yt, yp,
                 pos_min_threshold=0.05,
                 neg_max_threshold=0.01,
                 required_min_pos_counts=2.5,
                 binsizes=[1, 2, 4, 10, 20, 50]):
    """
    Evaluate the profile in terms of auPRC
    Args:
      yt: matrix of true profile (counts) [regions x position]
      yp: matrix of predicted profile (fractions) [regions x position]
      pos_min_threshold: fraction threshold above which the position is considered to be a positive
      neg_max_threshold: fraction threshold bellow which the position isconsidered to be a negative
      required_min_pos_counts: smallest number of reads the peak should be supported by.
          + All regions where `pos_min_threshold` of the total reads would be
          less than `required_min_pos_counts` are excluded
    """

    # The filtering criterion assures that each position in the positive class is
    # supported by at least required_min_pos_counts  of reads.
    do_eval = yt.sum(axis=1) > required_min_pos_counts / pos_min_threshold

    # make sure everything sums to one
    yp = yp / yp.sum(axis=1, keepdims=True)
    fracs = yt / yt.sum(axis=1, keepdims=True)

    yp_random = permute_array(permute_array(yp[do_eval], axis=1), axis=0)
    out = []
    for binsize in binsizes:
        is_peak = (fracs >= pos_min_threshold).astype(float)
        ambigous = (fracs < pos_min_threshold) & (fracs >= neg_max_threshold)
        is_peak[ambigous] = -1
        y_true = np.ravel(bin_counts_amb(is_peak[do_eval], binsize))

        imbalance = np.sum(y_true == 1) / np.sum(y_true >= 0)
        n_positives = np.sum(y_true == 1)
        n_ambigous = np.sum(y_true == -1)
        frac_ambigous = n_ambigous / y_true.size

        # Define bins by the maximum summit value within them.
        try:
            res = auprc(y_true,
                        np.ravel(bin_counts_max(yp[do_eval], binsize)))
            res_random = auprc(y_true,
                               np.ravel(bin_counts_max(yp_random, binsize)))
        except Exception:
            res = np.nan
            res_random = np.nan

        out.append({"binsize": binsize,
                    "auprc": res,
                    "random_auprc": res_random,
                    "n_positives": n_positives,
                    "frac_ambigous": frac_ambigous,
                    "imbalance": imbalance
                    })

    return pd.DataFrame.from_dict(out)

def auprc(y_true, y_pred):
    """
    Kudos to Ziga Avsec: bpnet.bpnet.metrics.py#L480
    Area under the precision-recall curve
    Remove numbers classified as -1 (ambiguious regions) and NaN
    """
    import sklearn.metrics as skm

    y_true, y_pred = _mask_value_nan(y_true, y_pred)
    return skm.average_precision_score(y_true, y_pred)

def _mask_nan(y_true, y_pred):
    """
    Kudos to Ziga Avsec: bpnet.bpnet.metrics.py#L435
    Remove NaNs
    """
    mask_array = ~np.isnan(y_true)
    if np.any(np.isnan(y_pred)):
        print("WARNING: y_pred contains {0}/{1} np.nan values. removing them...".
              format(np.sum(np.isnan(y_pred)), y_pred.size))
        mask_array = np.logical_and(mask_array, ~np.isnan(y_pred))
    return y_true[mask_array], y_pred[mask_array]


def _mask_value(y_true, y_pred, mask=-1):
    """
    Kudos to Ziga Avsec: bpnet.bpnet.metrics.py#L421
    Remove desired mask values
    """
    mask_array = y_true != mask
    return y_true[mask_array], y_pred[mask_array]


def _mask_value_nan(y_true, y_pred, mask=-1):
    """
    Kudos to Ziga Avsec: bpnet.bpnet.metrics.py#L430
    Remove desired mask values and NaNs
    """
    y_true, y_pred = _mask_nan(y_true, y_pred)
    return _mask_value(y_true, y_pred, mask)

def bin_counts_max(x, binsize=2):
    """
    Kudos to Ziga Avsec: bpnet.bpnet.metrics.py#L37
    Bin the counts where x is an array of [regions x profile]
    """
    if binsize == 1:
        return x
    outlen = x.shape[1] // binsize
    xout = np.zeros((x.shape[0], outlen))
    for i in range(outlen):
        xout[:, i] = x[:, (binsize * i):(binsize * (i + 1))].max(1)
    return xout

def bin_counts_amb(x, binsize=2):
    """
    Kudos to Ziga Avsec: bpnet.bpnet.metrics.py#L50
    Bin the counts (ambiguous -1 values included) where x is an array of [regions x profile]
    """
    if binsize == 1:
        return x
    outlen = x.shape[1] // binsize
    xout = np.zeros((x.shape[0], outlen)).astype(float)
    for i in range(outlen):
        iterval = x[:, (binsize * i):(binsize * (i + 1))]
        has_amb = np.any(iterval == -1, axis=1)
        has_peak = np.any(iterval == 1, axis=1)
        # if no peak and has_amb -> -1
        # if no peak and no has_amb -> 0
        # if peak -> 1
        xout[:, i] = (has_peak - (1 - has_peak) * has_amb).astype(float)
    return xout


def permute_array(arr, axis=0):
    """
    Kudos to Ziga Avsec: bpnet.bpnet.stats.py#L149
    Permute array along a certain axis
    Args:
      arr: numpy array
      axis: axis along which to permute the array
    """
    if axis == 0:
        return np.random.permutation(arr)
    else:
        return np.random.permutation(arr.swapaxes(0, axis)).swapaxes(0, axis)


##################################################################
#Collect counts for plotting together
##################################################################

def collect_counts(regions_path, profileA, profileB,
                   profileA_name, profileB_name,
                   valid_chroms = None, test_chroms = None,
                   seqlenth = 1000, exclude_zero_profiles = False):
    """
    Kudos to Zahoor, code modified from his basepairmodels.basepairmodels.cli.metrics.py
    Description:
        + collect count information across a set of regions,
        + these need to be normal profiles, not logits, and not exponentiated scalars over the region
    Inputs: 2 bigwigs and a set of .narrowPeak regions
    Outputs: data.frame of collected counts
    """
    import pyBigWig
    from tqdm import tqdm
    from scipy.ndimage import gaussian_filter1d
    from scipy.spatial.distance import jensenshannon
    from scipy.stats import pearsonr, spearmanr, multinomial
    from scipy.special import logsumexp

    # Read in information
    peaks_df = pd.read_csv(regions_path, sep='\t', header=None,
                           names=['chrom', 'st', 'end', 'name', 'score',
                                  'strand', 'signal', 'p', 'q', 'summit'])

    # create new column for peak pos
    peaks_df['summit_pos'] = peaks_df['st'] + peaks_df['summit']

    # create new column for start pos
    peaks_df['start_pos'] = peaks_df['summit_pos'] - seqlenth // 2

    # create new column for end pos
    peaks_df['end_pos'] = peaks_df['summit_pos'] + seqlenth // 2

    # select only the chrom & summit positon columns
    allPositions = peaks_df[['chrom', 'start_pos', 'end_pos']]
    allPositions = allPositions.reset_index(drop=True)

    # open the two bigWig files
    try:
        bigWigProfileA = pyBigWig.open(profileA)
        bigWigProfileB = pyBigWig.open(profileB)

    except Exception as e:
        logging.error("Problems occurred when opening one of the input files: "
                      "{}".format(str(e)))

    # for pearson on counts
    countsA = []
    countsB = []

    # initialize arrays to hold metrics values
    array_len = len(allPositions.index)

    for idx, row in tqdm(allPositions.iterrows(), total=allPositions.shape[0]):

        #Collect region information
        chrom = row['chrom']
        start = row['start_pos']
        end = row['end_pos']

        #Extract bigwig profiles
        try:
            profileA = np.nan_to_num(np.array(
                bigWigProfileA.values(chrom, start, end)))
            profileB = np.nan_to_num(np.array(
                bigWigProfileB.values(chrom, start, end)))
        except Exception as e:
            raise NoTracebackException(
                "Error retrieving values {}, {}, {}".format(chrom, start, end))

        #Compute counts information
        valsCountsA = np.sum(profileA)
        valsCountsB = np.sum(profileB)

        # add to the counts list
        countsA.append(np.sum(valsCountsA))
        countsB.append(np.sum(valsCountsB))

    #Collect information into a pd.df
    counts_df = pd.DataFrame([countsA, countsB]).transpose()
    counts_df.columns = [profileA_name,profileB_name]
    counts_df['log10_' + profileA_name] = np.log(counts_df[profileA_name])
    counts_df['log10_' + profileB_name] = np.log(counts_df[profileB_name])
    counts_df['chrom'] = peaks_df.chrom

    #Mark chromosome categories
    if valid_chroms or test_chroms:
        counts_df['chrom_cat'] = 'train'
        counts_df.loc[counts_df['chrom'].isin(valid_chroms), 'chrom_cat'] = 'valid'
        counts_df.loc[counts_df['chrom'].isin(test_chroms), 'chrom_cat'] = 'test'

    return(counts_df)
