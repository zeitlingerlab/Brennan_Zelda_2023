#General []_to_[] functions associated with BPNet
import sys
import os
import math
import warnings
import pandas as pd
import numpy as np
import random
from collections import Counter
from pybedtools import BedTool
from bpnet.cli.contrib import ContribFile
from bpnet.BPNet import BPNetSeqModel
from bpnet.data import NumpyDataset #creates a set of nested arrays
from kipoi.data import Dataset

warnings.filterwarnings("ignore")
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False

#Round to nearest base
def myround(x, base=50):
    return base * round(x/base)
def myceiling(x, base = 10):
    return int(math.ceil(x / base)) * base
def myfloor(x, base = 10):
    return int(math.floor(x / base)) * base

"""
Convert 1-hot encoded sequence to a sequence string
"""
def one_hot_decode(arr): #arr=[seqlen x 4]
    _, seq_idx, = np.where(arr==1) #Get value
    decoder = {0:'A', 1:'C', 2:'G', 3:'T'}
    seq = ''.join([decoder[i] for i in seq_idx])
    return seq

"""
Convert pd.df with coordinate headers to a set of intervals.
bed3 = True, Assuming df has a `chrom`, `start` and `end`
bed3 = False, Assuming df has a `strand`, `name`, and `score`
"""
def df_to_intervals(df, bed3 = True):
    assert ('chrom' in df) and ('start' in df) and ('end' in df), (
        "Input pd.df does not contain either chrom, start or end columns.")
    if bed3:
        coord_str_list = ['\t'.join([df.chrom[coord], df.start[coord], df.end[coord]]) for coord in range(df.shape[0])]
        coord_str = '\n'.join(coord_str_list)
        coord_interv = BedTool(coord_str, from_string = True)
    else:
        assert ('strand' in df) and ('name' in df) and ('score' in df), (
            "Input pd.df does not contain either strand, name or score columns. Consider changing bed3 = True.")
        coord_str_list = ['\t'.join([df.chrom[coord], df.start[coord], df.end[coord],
                                     df.name[coord], df.score[coord], df.strand[coord]]) for coord in range(df.shape[0])]
        coord_str = '\n'.join(coord_str_list)
        coord_interv = BedTool(coord_str, from_string = True)
    return(coord_interv)

#Given a NumpyDataset object for a single region, extract predictions as tidy pd.df
def pred_npd_to_tidy_df(npd, npd_idx, pred_type = 'pred'):
    pred_across_tasks = pd.DataFrame()
    for task in list(npd[npd_idx].keys()): #For each task...
        assert npd[npd_idx][task][pred_type].shape[1]==2, "There are not two strands of profiles. Is this ChIP-seq data?"
        pred_df = pd.DataFrame(npd[npd_idx][task][pred_type]).reset_index()
        pred_df.columns = ['position', 'pos', 'neg']
        pred_df.neg = pred_df.neg * -1 #flip strand on negative profiles
        pred_df['task'] = task
        pred_across_tasks = pred_across_tasks.append(pred_df)
    pred_across_tasks = pd.melt(pred_across_tasks, id_vars = ['position', 'task'],
                                var_name = 'strand', value_name = 'preds_profile')
    return(pred_across_tasks)

#Given a NumpyDataset object for a single region, extract contribs based on contrib_type as tidy pd.df
def contrib_npd_to_tidy_df(npd, npd_idx, contrib_type = 'profile'):
    contrib_across_tasks = pd.DataFrame()
    for task in list(npd[npd_idx].keys()): #For each task...
        assert npd[npd_idx][task]['contrib'][contrib_type].shape[1]==4, "There are not 4 tracks. Is this 1-hot encoded?"
        #flatten 1-hot encoded results
        contrib_df = pd.DataFrame(npd[npd_idx][task]['contrib'][contrib_type].sum(axis = 1)).reset_index()
        contrib_df.columns = ['position', 'contrib' + '_' + contrib_type]
        contrib_df['strand'] = 'pos'
        contrib_df['task'] = task
        contrib_across_tasks = contrib_across_tasks.append(contrib_df)
    return(contrib_across_tasks)

# Given a npd of contribution scores and a task to parse, get a [position x ACGT(4)] pd.df
# `genomic_region_start` is a scalar value of the start site of the region of interest.
# If selected, axis will use that.
def contrib_npd_to_df(npd, region_idx, task, contrib_type = 'profile'):
    contrib_df = pd.DataFrame(npd[region_idx][task]['contrib'][contrib_type],
                              columns = ['A','C','G','T'])
    return(contrib_df)

#Convert predictions to a pandas df for custom visualization
#Get a tidy pd.df from the output of the function `BPNetSeqModel.predict_all(seqs)`
#predictions = BPNetSeqModel.predict_all(seqs)[region]['pred']
def tidy_bpnet_predictions_nexus(predictions, tasks, flip_strands = False):
    df_all = pd.DataFrame()
    for task in tasks:
        pred = predictions[task]
        if flip_strands:
            pos_pred = np.flip(pred[:, 1])
            neg_pred = np.flip(pred[:, 0] * -1)

        else:
            pos_pred = pred[:, 0]
            neg_pred = pred[:, 1] * -1

        pos_df = pd.DataFrame(pos_pred)
        pos_df.columns = ['prediction']
        pos_df['position'] = pos_df.axes[0] #get region
        pos_df['task'] = task
        pos_df['strand'] = '+'
        neg_df = pd.DataFrame(neg_pred)
        neg_df.columns = ['prediction']
        neg_df['position'] = neg_df.axes[0] #get region
        neg_df['task'] = task
        neg_df['strand'] = '-'
        #Combine
        df = pos_df.append(neg_df)
        df_all = df_all.append(df)
    return(df_all)

#Function to save contribution information with 4 columns and a task (does not use hypoth. contrib)
def tidy_bpnet_contributions(seq, contrib, tasks, contrib_type = 'profile', flip_strands = False):
    contrib_df = pd.DataFrame()
    for task in tasks:
        c = contrib[f'{task}/{contrib_type}'] * seq
        if flip_strands:
            df = pd.DataFrame(np.flip(c, axis = 0), columns = ['T','G','C','A'])
            df = df[['A','C','G','T']]
        else:
            df = pd.DataFrame(c, columns = ['A','C','G','T'])
        df['task'] = task
        df['position'] = range(seq.shape[0])
        contrib_df=contrib_df.append(df)
    return(contrib_df)
