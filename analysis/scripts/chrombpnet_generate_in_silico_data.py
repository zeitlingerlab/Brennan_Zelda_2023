#! /home/mw2098/anaconda3/envs/chrombpnet/bin/python

##################################################################
# Computational setup
##################################################################
#Modules
import os
import json
import sys
import pickle
import random
import pandas as pd
import numpy as np
import tensorflow as tf
from tqdm import tqdm

#Custom functions
sys.path.insert(0, f'scripts/py')
from chrombpnet_perturb_functions import predict_injected_seq
from chrombpnet_predict_functions import load_chrombpnet, predict_chrombpnet, one_hot_encode_sequences

#Independent Variables
BASE_DIR='/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/analysis/'
os.chdir(BASE_DIR)
INPUT_SEQLEN=2114
OUTPUT_SEQLEN=1000
trials = 512
primary_position = 300
step  = 1

#Dependent variables
MODEL_DIR=BASE_DIR + 'chrombpnet/models/model_fold1'
os.environ['LD_LIBRARY_PATH']='/usr/lib64:/usr/local/cuda-11.5/lib64:/usr/local/cuda-11.5/extras/CUPTI/lib64:/lib/nccl/cuda-11.5:/usr/local/cuda-11.0/lib64:/usr/local/cuda-11.0/extras/CUPTI/lib64:/lib/nccl/cuda-11.0:/n/apps/CentOS7/install/rbenv-0.3.0/rbenv/libexec:/n/apps/CentOS7/lib64::/n/apps/CentOS7/install/guppy-2.3.1/lib:/n/apps/CentOS7/lib:/n/apps/CentOS7/proteomics/lib'

##################################################################
# Generate random seqs
##################################################################

def generate_random_seq(seqlen, weights = [.25, .25, .25, .25]):
    """
    Purpose: Generate a random DNA sequence of a specified length.
    """
    import random
    seq = random.choices(['A','C','G','T'], weights = weights, k=seqlen)
    return(''.join(seq))

random.seed(10)
random_seqs = [generate_random_seq(seqlen = INPUT_SEQLEN) for i in range(trials)]
random_seqs_1he = one_hot_encode_sequences(random_seqs)
random_seqs_1he.shape

##################################################################
# Identify motifs
##################################################################

motifs_dict = {
    'Zld': 'CAGGTAG',
    'Zld-lowaff': 'CAGGTA',
    'Zld-lowestaff': 'TAGGTAG',
    'Gaf': 'GAGAGAGAGAGAGAGAG',
    'Cad': 'TTTTATGGCC',
    'Dl': 'GGGAAAACCC',
    'Twi': 'AACACATGTT',
    'Bcd': 'TTAATCC'
}
timepoints = ['1to15','15to2','2to25','25to3']

#Identify other parameters
distances_dict = {
             '0to100': list(range(primary_position + 20, primary_position + 100, step)),
             '100to200': list(range(primary_position + 100, primary_position + 200, step)),
             '200to300': list(range(primary_position + 200, primary_position + 300, step)),
             '300to400': list(range(primary_position + 300, primary_position + 400, step))
             }

##################################################################
# Generate data
##################################################################

for t in timepoints:
    print(t)
    for motifA,seqA in tqdm(motifs_dict.items()):
        for motifB,seqB in motifs_dict.items():
            for d_name, d_set in distances_dict.items():

                if os.path.exists(f'tsv/perturbs/accessibility/insilico/in_silico_injections_paired_profiles_timepoint_{t}_motifA_{motifA}_motifB_{motifB}_distance_{d_name}.pkl') or \
                    os.path.exists(f'tsv/perturbs/accessibility/insilico/in_silico_injections_paired_profiles_timepoint_{t}_motifA_{motifB}_motifB_{motifA}_distance_{d_name}.pkl'):
                    print(f'Combination of motifA {motifA} and motifB {motifB} has already been tested')
                else:
                    ##################################################################
                    # Load model
                    ##################################################################
                    tf.keras.backend.clear_session()
                    model_path = f'{MODEL_DIR}_{t}.h5'
                    bias_model_path = f'{MODEL_DIR}_{t}.adjusted_bias_model.h5'
                    model_chrombpnet, model_bias = load_chrombpnet(model_path = model_path, bias_model_path = bias_model_path)

                    ##################################################################
                    # Predict data
                    ##################################################################
                    pred_dict = predict_injected_seq(model_chrombpnet = model_chrombpnet, model_bias = model_bias,
                                                    primary_motif = seqA,
                                                    primary_position = primary_position,
                                                    secondary_motif = seqB,
                                                    secondary_positions = d_set,
                                                    input_seqs = random_seqs_1he,
                                                    input_seqlen = INPUT_SEQLEN,
                                                    trials = trials, batch_size = 16)

                    with open(f'tsv/perturbs/accessibility/insilico/in_silico_injections_paired_profiles_timepoint_{t}_motifA_{motifA}_motifB_{motifB}_distance_{d_name}.pkl', 'wb') as f:
                        pickle.dump(pred_dict, f, pickle.HIGHEST_PROTOCOL)
