
"""
Melanie Weilert
Jul 2022
Purpose: After training, collect fold-interation test chromosome data for performance plotting.
"""

# Setup
import warnings;warnings.filterwarnings("ignore")
from tensorflow.python.util import deprecation; deprecation._PRINT_DEPRECATION_WARNINGS = False

#Modules
import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import plotnine
from glob import glob
from plotnine import *
from itertools import product, compress
from pybedtools import BedTool
from keras import backend as K

from datetime import datetime

from bpnet.utils import read_json, create_tf_session
from bpnet.dataspecs import DataSpec
from bpnet.datasets import StrandedProfile
from bpnet.extractors import StrandedBigWigExtractor
from bpnet.BPNet import BPNetSeqModel
from bpnet.metrics import eval_profile

#Setup
os.chdir('/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/analysis/')

#Define datasets
combined_path_dict = {
    'Zld': {'positive': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_zld_mbt_nexus_combined_positive.bw',
            'negative': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_zld_mbt_nexus_combined_negative.bw'},
    'Dl': {'positive': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_dl_mbt_nexus_combined_positive.bw',
           'negative': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_dl_mbt_nexus_combined_negative.bw'},
    'Twi': {'positive': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_twi_mbt_nexus_combined_positive.bw',
            'negative': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_twi_mbt_nexus_combined_negative.bw'},
    'Cad': {'positive': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_cad_mbt_nexus_combined_positive.bw',
            'negative': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_cad_mbt_nexus_combined_negative.bw'},
    'Bcd': {'positive': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_bcd_mbt_nexus_combined_positive.bw',
            'negative': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_bcd_mbt_nexus_combined_negative.bw'},
    'GAF': {'positive': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_gaf_mbt_nexus_combined_positive.bw',
            'negative': f'/l/Zeitlinger/ZeitlingerLab/Manuscripts/Zelda_and_Nucleosomes/Analysis/data/bw/dm6/nexus/combined/orer_gaf_mbt_nexus_combined_negative.bw'}
}

print([[os.path.exists(v2) for k2,v2 in v1.items()] for k1,v1 in combined_path_dict.items()])

rep_path_dict = {
    'Zld': {
        '1': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_zld_mbt_nexus_1_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_zld_mbt_nexus_1_negative.bw'
        },
        '2': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_zld_mbt_nexus_3_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_zld_mbt_nexus_3_negative.bw'
        }
    },
    'Dl': {
        '1': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_dl_mbt_nexus_1_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_dl_mbt_nexus_1_negative.bw'
        },
        '2': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_dl_mbt_nexus_3_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_dl_mbt_nexus_3_negative.bw'
        }
    },
    'Twi': {
        '1': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_twi_mbt_nexus_1_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_twi_mbt_nexus_1_negative.bw'
        },
        '2': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_twi_mbt_nexus_3_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_twi_mbt_nexus_3_negative.bw'
        }
    },
    'Cad': {
        '1': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_cad_mbt_nexus_2_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_cad_mbt_nexus_2_negative.bw'
        },
        '2': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_cad_mbt_nexus_3_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_cad_mbt_nexus_3_negative.bw'
        }
    },
    'Bcd': {
        '1': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_bcd_mbt_nexus_2_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_bcd_mbt_nexus_2_negative.bw'
        },
        '2': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_bcd_mbt_nexus_3_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_bcd_mbt_nexus_3_negative.bw'
        }
    },
    'GAF': {
        '1': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_gaf_mbt_nexus_1_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_gaf_mbt_nexus_1_negative.bw'
        },
        '2': {
            'positive': f'../data/bw/dm6/nexus/individual/orer_gaf_mbt_nexus_4_positive.bw',
            'negative': f'../data/bw/dm6/nexus/individual/orer_gaf_mbt_nexus_4_negative.bw'
        }
    }
}


print([[[os.path.exists(v3) for k3,v3 in v2.items()] for k2,v2 in v1.items()] for k1,v1 in rep_path_dict.items()])

#Collect data
data_dict = {}
for fold in range(1,4):
    name = f'fold{fold}'
    data_dict[name]={}

    #Load the model
    K.clear_session()
    model = BPNetSeqModel.from_mdir(f'bpnet/models/optimized_model/fold{fold}')

    #Get data loader specs for validation set
    dataspec = DataSpec.load(f'bpnet/dataspec/dataspec.yml')
    tasks = list(dataspec.task_specs.keys())

    #Get the configurations generated by the model
    config = read_json(f'bpnet/models/optimized_model/fold{fold}/config.gin.json')
    excl_chromosomes = config['bpnet_data.exclude_chr']
    valid_chromosomes = config['bpnet_data.valid_chr']
    test_chromosomes = config['bpnet_data.test_chr']
    print(f'Fold{fold} exluded chromosomes:', excl_chromosomes)
    print(f'Fold{fold} validation chromosomes:', valid_chromosomes)
    print(f'Fold{fold} test chromosomes:', test_chromosomes)

    #Specify the custom object for loading stranded profiles from the test chromosome.
    dl = StrandedProfile(dataspec,
                         incl_chromosomes = test_chromosomes,
                         excl_chromosomes = excl_chromosomes,
                         peak_width=config['bpnet_data.peak_width'],
                         seq_width=config['bpnet_data.seq_width'],
                         shuffle=False)
    data = dl.load_all(num_workers=8)

    #Get regions for test chromosomes
    regions_df = pd.DataFrame([list(data['metadata']['range']['chr']),
                                  list(data['metadata']['range']['start']),
                                  list(data['metadata']['range']['end']),
                                  list(data['metadata']['range']['strand'])]).transpose()
    regions_bed = BedTool.from_dataframe(regions_df)

    #Get obs  values from the test chromosom
    y_obs = data['targets']

    #Get predicted values from the test chromosome
    y_pred_seqmodel = model.seqmodel.predict(data['inputs']['seq']) #Provides counts and logits
    y_pred = {task: y_pred_seqmodel[f'{task}/profile'] * np.exp(y_pred_seqmodel[f'{task}/counts'][:, np.newaxis])
                    for task in model.tasks} #Provides counts*profile

    # Extract replicate values for each region for comparison
    rep1_obs = {t: np.stack([StrandedBigWigExtractor(bigwig_file = rep_path_dict[t]['1']['positive']).extract(regions_bed),
                             StrandedBigWigExtractor(bigwig_file = rep_path_dict[t]['1']['negative']).extract(regions_bed)], axis = 2)
                for t in tqdm(tasks)}
    rep2_obs = {t: np.stack([StrandedBigWigExtractor(bigwig_file = rep_path_dict[t]['2']['positive']).extract(regions_bed),
                             StrandedBigWigExtractor(bigwig_file = rep_path_dict[t]['2']['negative']).extract(regions_bed)], axis = 2)
                for t in tqdm(tasks)}

    #Generate metapeak from profiles
    y_mp_2d = {t: np.mean(y_obs[f'{t}/profile'], axis = 0) for t in tasks}
    y_mp = {t: np.stack([y_mp_2d[t] for i in range(y_obs[f'{t}/profile'].shape[0])], axis = 0) for t in tasks}

    #Take the absolute value of the replicates introduced because ChIP-nexus bigwigs have the strand information encoded as axes.
    rep1_obs = {k: np.abs(v) for k,v in rep1_obs.items()}
    rep2_obs = {k: np.abs(v) for k,v in rep2_obs.items()}

    data_dict[name]['y_obs'] = y_obs
    data_dict[name]['y_pred'] = y_pred
    data_dict[name]['rep1_obs'] = rep1_obs
    data_dict[name]['rep2_obs'] = rep2_obs
    data_dict[name]['y_mp'] = y_mp
    data_dict[name]['test_chrom'] = test_chromosomes
    data_dict[name]['valid_chrom'] = valid_chromosomes
    data_dict[name]['exclude_chrom'] = excl_chromosomes

#Write data as .pkl file
import pickle
with open('bpnet/pkl/obs_vs_pred_metrics.pkl', 'wb') as handle:
    pickle.dump(data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
