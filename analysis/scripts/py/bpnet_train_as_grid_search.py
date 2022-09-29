
"""
Melanie Weilert
January 2021
Purpose: This code is meant to train BPNet models based on input grid search.

Intructions:
1. Activate `bpnet` conda environment.
2. `python /n/projects/mw2098/shared_code/bpnet/bpnet_train_as_grid_search.py [options]`
3. For help: `python /n/projects/mw2098/shared_code/bpnet/bpnet_train_as_grid_search.py -h`
"""

# Setup
import os
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
from bpnet.utils import add_file_logging, read_pkl, create_tf_session
from argh.decorators import named, arg

import warnings
warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

parser = OptionParser()
parser.add_option("-d", "--dataspec",
                  help="Path to the dataspec.yml file specifying BPNet input data.")
parser.add_option("-g", "--config",
                  help="Path to the config.gin file specifying BPNet default parameters.")
parser.add_option("-o", "--output-directory",
                  help="Directory that trained model subdirectories will be saved under.")
parser.add_option("-m", "--manually-run", default = False, action="store_true",
                help="Boolean indicating whether to run the .sh script generated or to leave it. [default = %default]")
parser.add_option("-x", "--memfrac-gpu", default = .6, type = "float",
                  help="what fraction of the GPU memory to use [default = %default]")
parser.add_option("-p", "--tmp-dir", default = 'tmp', type = "str",
                  help="Temporary directory to write the script to.")
parser.add_option("-f", "--filters", default = '64', type = "str",
                  help="Number of filters used in convolutions. Multiple values can be given in a comma-separated format. [default = %default]")
parser.add_option("-c", "--conv-kernel-size", default = '25', type = "str",
                  help="Kernel length of the first convolutional layer. Multiple values can be given in a comma-separated format. [default = %default]")
parser.add_option("-t", "--tconv-kernel-size", default = '25', type = "str",
                  help="Kernel length of the last deconvolutional layer when generating profile predictions. Multiple values can be given in a comma-separated format. [default = %default]")
parser.add_option("-n", "--n-dil-layers", default = '9', type = "str",
                  help="Number of dilated convolutional `kmer` layers. Multiple values can be given in a comma-separated format. [default = %default]")
parser.add_option("-l", "--loss-weight", default = '10', type = "str",
                  help="Also known as lambda, the weighted scalar multiplied to counts loss. Multiple values can be given in a comma-separated format. [default = %default]")
parser.add_option("-r", "--learning-rate", default = '0.004', type = "str",
                  help="Learning rate of the Adam optimizer during training. Multiple values can be given in a comma-separated format. [default = %default]")
parser.add_option("-w", "--seq-width", default = '1000', type = "str",
                  help="Width of the sequence window evaluated by the model. Multiple values can be given in a comma-separated format. [default = %default]")
parser.add_option("-a", "--num-workers", default = 16, type = int,
                  help="Number of works to load the data with while training. If training fails, try c=1. [default = %default]")

(options, args) = parser.parse_args()

def bpnet_train_as_grid_search(dataspec,
                          config,
                          output_directory,
                          memfrac_gpu=0.6,
                          filters='64',
                          conv_kernel_size='25',
                          tconv_kernel_size='25',
                          n_dil_layers='9',
                          loss_weight='10',
                          learning_rate='0.004',
                          seq_width='1000',
                          manually_run=False,
                          tmp_dir = 'tmp',
                          num_workers = 16):
    """
    Train multiple BPNet models in combinations as per the given hyperparameters.
    """
    import os
    import sys
    import subprocess
    import numpy as np
    import pandas as pd
    from datetime import datetime
    from itertools import product, compress

    # Make directories
    output_directory = os.path.dirname(output_directory)
    add_file_logging(output_directory, logger, 'bpnet-train-grid-search')
    os.makedirs(output_directory, exist_ok=True)

    logger.info(f'Collecting hyperparameter settings...')
    hp_dict = {
            "seq_width": seq_width.replace(' ','').split(','),
            "lr": learning_rate.replace(' ','').split(','),
            "lambda": loss_weight.replace(' ','').split(','),
            "n_dil_layers": n_dil_layers.replace(' ','').split(','),
            "conv_kernel_size": conv_kernel_size.replace(' ','').split(','),
            "tconv_kernel_size": tconv_kernel_size.replace(' ','').split(','),
            "filters": filters.replace(' ','').split(',')
    }

    # Find all combinations of different hyperparameters
    def expand_grid(dictionary):
        #Kudos: https://stackoverflow.com/questions/12130883/r-expand-grid-function-in-python
        return pd.DataFrame([row for row in product(*dictionary.values())], columns=dictionary.keys())
    trials_df = expand_grid(hp_dict)

    #Collect hyperparameter settings
    hps = []
    for idx, row in trials_df.iterrows():
        keep = [r!='' for r in row.tolist()]
        hp = list(zip(trials_df.columns[keep], list(compress(row.tolist(), keep))))
        hp = ';'.join(['='.join(tups) for tups in hp])
        hps.append(hp)

    #Collect hyperparameter names
    hp_names = [hp.replace('=','').replace(';','-') for hp in hps]

    tmp_dir = 'tmp'
    date = datetime.now().strftime("%Y%m%d%H%M%S")

    cmds = [f"bpnet train --num-workers {num_workers} --vmtouch --config={config} --override='{i[0]}' --run-id {i[1]} --memfrac-gpu {memfrac_gpu} {dataspec} {output_directory}" for i in zip(hps, hp_names)]
    script = ['#!bin/bash'] + cmds
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        logger.info(f'Writing commands to tmp script...')
    with open(f'{tmp_dir}/{date}_bpnet_grid_search_commands_tmp.sh', "w") as outfile:
        outfile.write("\n".join(script))

    #Run script
    if not manually_run:
        logger.info(f'Beginning hyperparameter grid search with {len(cmds)} trials...')
        for cmd in cmds:
            subprocess.run(cmd, shell = True)
    logger.info("Finished!")

    return None

bpnet_train_as_grid_search(dataspec = options.dataspec,
                          config = options.config,
                          output_directory = options.output_directory,
                          memfrac_gpu = options.memfrac_gpu,
                          filters = options.filters,
                          conv_kernel_size = options.conv_kernel_size,
                          tconv_kernel_size = options.tconv_kernel_size,
                          n_dil_layers = options.n_dil_layers,
                          loss_weight = options.loss_weight,
                          learning_rate = options.learning_rate,
                          seq_width = options.seq_width,
                          manually_run=options.manually_run,
                          tmp_dir = options.tmp_dir,
                          num_workers = options.num_workers)
