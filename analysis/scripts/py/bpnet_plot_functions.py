#General plotting functions associated with BPNet predictions
import sys
import os
import math
import warnings
import pandas as pd
import numpy as np
import random
import logomaker
import itertools
import matplotlib.pyplot as plt
from collections import Counter
from pybedtools import BedTool
from bpnet.cli.contrib import ContribFile
from bpnet.BPNet import BPNetSeqModel
from bpnet.data import NumpyDataset #creates a set of nested arrays
from kipoi.data import Dataset
from matplotlib import gridspec

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

# Iteration of bpnet.plot.tracks.plot_seqlet_box but with customizable box colors!
def plot_seqlet_box(seqlet, ax, box_fill = '#939393', alpha = .5, add_label=False):
    """
    Args:
      seqlet: object with start, end, name, strand attribues
      if Seqname is available, then we can plot it to the right position
    """
    from bpnet.plot.utils import (simple_yaxis_format, strip_axis, spaced_xticks,
                                  draw_box, spine_subset, draw_hline)
    xlim = ax.get_xlim()

    xmin = seqlet.start + 0.5
    xmax = seqlet.end + 0.5
    if xmax < 0 or xmin > xlim[1] - xlim[0]:
        return

    draw_box(xmin, xmax, ax, col = box_fill, alpha = alpha)  # trimmed pattern location

    if add_label:
        y = ax.get_ylim()[1] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.15

        ax.text(xmin + 0.5, y,
                s=seqlet.strand + str(seqlet.name),
                fontsize=6)




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














#Given a [position x ACGT(4)] pd.df, plot contribution scores
def plot_contrib_track(contrib_df, ax, xlimits = None, ylimits = None, ytitle = 'Contrib (profile)'):
    contrib_logo = logomaker.Logo(contrib_df, ax = ax)
    contrib_logo.ax.set_ylabel(ytitle, labelpad=-1)
    contrib_logo.ax.set_xlim(xlimits)
    contrib_logo.ax.set_ylim(ylimits)
    contrib_logo.ax.ticklabel_format(axis = 'both', style = 'plain', useOffset = False)
    return(contrib_logo)

# Given a Logo object, loop through motif coordinates to plot the correct motif.
def plot_contrib_motif_box(contrib_logo, motif_df, motif_color = "gray"):
    contrib_logo_setting = contrib_logo.highlight_position_range(pmin = np.asscalar(np.array(motif_df['plot_start'])),
                                                                 pmax = np.asscalar(np.array(motif_df['plot_end'])),
                                                                 color=motif_color, alpha=.2)
    return(contrib_logo_setting)

"""
Purpose: For a region, plot the predictions
Inputs:
    preds_df = tidy pd.df of the prediction profiles (REQUIRES 'plot_position', 'preds_profile', 'strand' column)
    ax = plot object to add to
Outputs:
    ax = plot object to save as
"""
def plot_nexus_prediction(preds_df, ax, xlimits = None, ylimits = None):

    #Subset regions to appease matplotlib
    preds_pos_df = preds_df[preds_df['strand'] == 'pos']
    preds_neg_df = preds_df[preds_df['strand'] == 'neg']

    #Plot profiles
    ax.plot(preds_pos_df['plot_position'], preds_pos_df['preds_profile'], color = 'black')
    ax.plot(preds_neg_df['plot_position'], preds_neg_df['preds_profile'], color = '#383838')

    if xlimits is not None: ax.set_xlim(xlimits[0], xlimits[1])
    if ylimits is not None: ax.set_ylim(ylimits[0], ylimits[1])

    # Add aesthetics
    ax.set_ylabel('BPNet predictions')
    return(ax)

"""
Purpose: For a region, plot the cwm-scanned motifs
Inputs:
    motif_coords_df = tidy pd.df of thecwm-scanned coordinates (REQUIRES 'plot_start', 'plot_end' column)
    ax = plot object to add to
Outputs:
    ax = plot object to save as
"""

def plot_pred_motif_boxes(motif_coords_df, ax, color = '#B0B0B0'):
    #Plot motif boxes
    for motif_idx, motif in motif_coords_df.iterrows():
        ax.axvspan(motif['plot_start'], motif['plot_end'], facecolor=color, alpha = .2)
    return(ax)
