"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions to predicting using chrombpnet model
"""
import sys
import os
import json
import pandas as pd
import numpy as np
import tensorflow as tf

sys.path.insert(0, f'/home/mw2098/bin/chrombpnet-lite/src/')
from metrics import softmax
from utils.loss import multinomial_nll
from utils import one_hot

def load_chrombpnet(model_path, bias_model_path):
    """
    Purpose: Load the keras model.h5 from `chrombpnet`. Requires CustomObject.
    """
    physical_devices = tf.config.list_physical_devices('GPU')
    try:
        tf.config.experimental.set_memory_growth(physical_devices[0], True)
    except:
        # Invalid device or cannot modify virtual devices once initialized.
        pass

    with tf.keras.utils.CustomObjectScope({'multinomial_nll': multinomial_nll, 'tf': tf}):
        model_bias = tf.keras.models.load_model(bias_model_path)
        model_chrombpnet = tf.keras.models.load_model(model_path)
    return model_chrombpnet, model_bias

def predict_chrombpnet(model_chrombpnet, model_bias, seqs, batch_size=128, verbose = False, no_bias = False):
    """
    Purpose: Load the keras model.h5 from `chrombpnet`. Requires CustomObject.
    Inputs:
        + model_chrombpnet: ChromBPnet model with true ATAC-seq predictions
        + model_bias: Bias model attached to ChromBPnet
        + seqs : [region x position x 4] array of one-hot encoded sequences
        + no_bias: if True, do not predict bias and incorporate it into output predictions
    Note: will automatically detect if sequences are one-hot encoded
    """
    sys.path.insert(0, f'/home/mw2098/bin/chrombpnet-lite/src/')
    from metrics import softmax

    #If the sequence is not already one-hot-encoded, then assign./
    seqs = np.array(seqs)
    if (seqs.shape[-1] != 4) and len(seqs.shape)==1:
        seqs = one_hot_encode_sequences(seqs)
    # print(seqs.shape)

    if not no_bias:
        #predict bias on peaks and nonpeaks
        pred_bias_logits, pred_bias_logcts = \
                model_bias.predict(seqs,
                                   batch_size = batch_size,
                                   verbose=verbose)

        # predict chrombpnet on peaks and nonpeaks
        pred_logits, pred_logcts = \
                model_chrombpnet.predict([seqs,
                                          pred_bias_logits,
                                          pred_bias_logcts],
                                        batch_size=batch_size,
                                        verbose=verbose)

    # pred_logits_wo_bias, pred_logcts_wo_bias = \
    #         model_chrombpnet.predict([seqs,
    #                                   np.zeros_like(pred_bias_logits),
    #                                   np.zeros_like(pred_bias_logcts)],
    #                                   batch_size = batch_size,
    #                                   verbose=verbose)


    pred_logits_wo_bias, pred_logcts_wo_bias = \
            model_chrombpnet.predict([seqs,
                                      np.zeros((seqs.shape[0], model_bias.outputs[0].get_shape().as_list()[1])),
                                      np.zeros((seqs.shape[0], model_bias.outputs[1].get_shape().as_list()[1]))],
                                      batch_size = batch_size,
                                      verbose=verbose)

    #turn into profile
    if not no_bias:
        profile_w_bias = softmax(pred_logits) * (np.exp(pred_logcts)-1)

    profile_wo_bias = softmax(pred_logits_wo_bias) * (np.exp(pred_logcts_wo_bias)-1)

    if not no_bias:
        return(profile_w_bias, profile_wo_bias)
    else:
        return(None, profile_wo_bias)

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


# def softmax(array):
#     """
#     Kudos: https://machinelearningmastery.com/softmax-activation-function-with-python/
#     """
#     assert len(array.shape)==1, 'Input array is not a 1D vector.'
#     e = np.exp(array)
#     return(e / np.sum(e))

def pred_logits_to_softmax(preds_logits):
    """
    Purpose: Convert [regions x position x channel] logits output of `predict_basepairmodel` to softmax
    """
    assert len(preds_logits.shape)==3, 'Input array is not a 3D tensor [regions x position x task].'
    preds_softmax = np.empty(preds_logits.shape)
    for i in range(preds_logits.shape[0]):
        for j in range(preds_logits.shape[2]):
            preds_softmax[i, :, j] = softmax(preds_logits[i, :, j])
    return(preds_softmax)

def predict_from_seq_list_to_df(seqs_list, seqs_names, model_file, bias_file):
    """
    Purpose: Predict from list of sequences to a tidy pd.df
    Inputs:
        + seqs_list: list of sequences
        + seqs_names: list of names associated with sequences for unique identification
        + model_file: filepath to model .h5
        + bias_file: filepath to bias model .h5
    Output: pd.df of predictions in tidy format
    """
    #One-hot encode sequences
    seqs_1he = one_hot_encode_sequences(seqs_list)

    #Load model of interest
    print('Loading model...')
    chrombpnet_model, bias_model = load_chrombpnet(model_file, bias_file)

    #Predict with counts integrated
    _, preds_w_counts = predict_chrombpnet(model_chrombpnet = chrombpnet_model,
                                           model_bias = bias_model,
                                           seqs = seqs_1he,
                                           no_bias = True)

    preds_df = pd.DataFrame(preds_w_counts)
    preds_df['seq_name'] = seqs_names
    preds_df = preds_df.melt(id_vars = ['seq_name'], var_name = 'position', value_name = 'pred')
    preds_df = preds_df.sort_values(['seq_name', 'position'])
    return(preds_df)
