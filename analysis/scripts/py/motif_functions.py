"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions to modifying/generating motifs not covered by `basepairmodels.CLI`
"""

import os
import json
import pandas as pd
import numpy as np

def resize_coordinates(coords_df, width = 1, fix = 'start', chrom_column = 'chrom',
                       start_column = 'start', end_column = 'end', strand_column = 'strand'):
    """
    Purpose: Given a pd.df of coordinates (chrom, start, end, strand), defined by the inputs return resized coordinates.
        + Follow same convention in R of GenomicRanges::resize
        + The resizing will occur in a 5' -> 3' direction relative to the strand of each region.
    Inputs:
         + coords_df: pd.df of coordinates (chrom, start, end, strand)
         + width: how wide should the coordinates be? These are 0-based coordinates.
         + fix:
             +if 'start', extend fragments downstream while anchoring the START coordinate (relative to the orientation of the coordinate.)
             +if 'center', extend fragments upstream and downstream (relative to the center of the coordinate.)
             +if 'end', extend fragments upstream while anchoring the END coordinate (relative to the orientation of the coordinate.)
        + [parameter]_column: the column name of the pd.df that represents the column
    Output: pd.df of coordinates that are resized
    """
    import warnings
    warnings.filterwarnings("ignore")

    coords_df = coords_df.reset_index(drop = True)

    #Define metadata
    metadata_cols = coords_df.columns
    filter_cols = [chrom_column, start_column, end_column, strand_column]
    add_cols = [t for t in metadata_cols if t not in filter_cols]

    #Define new object
    new_starts = []
    new_ends = []

    #Anchor the 5' -> 3' oriented START coordinate and resize
    if fix=='start':
        for i,row in coords_df.iterrows():
            if row[strand_column] == '-':
                new_start = row[end_column] - width
                new_end = row[end_column]
            else:
                new_start = row[start_column]
                new_end = row[start_column] + width
            new_starts.append(new_start)
            new_ends.append(new_end)

    #Anchor the 5' -> 3' oriented END coordinate and resize
    elif fix=='end':
        for i,row in coords_df.iterrows():
            if row[strand_column] == '-':
                new_start = row[start_column]
                new_end = row[start_column] + width
            else:
                new_start = row[end_column] - width
                new_end = row[end_column]
            new_starts.append(new_start)
            new_ends.append(new_end)

    #Anchor the 5' -> 3' oriented CENTER coordinate and resize
    elif fix=='center':

        #Decide on upstream and downstream width allocations if odd width:
        if width % 2 == 0:
            width_up = width // 2
            width_down = width // 2
        else:
            width_up = width // 2
            width_down = (width // 2) + 1

        #Find centers, then assign starts and ends
        new_centers = []
        for i,row in coords_df.iterrows():
            if row[strand_column] == '-':
                center = row.end - ((row.end - row.start) // 2)
                new_start = center - width_down
                new_end = center + width_up
            else:
                center = row.start + ((row.end - row.start) // 2)
                new_start  = center - width_up
                new_end  = center + width_down
            new_starts.append(new_start)
            new_ends.append(new_end)
    else:
        print('fix parameter is not correct. Please refer to documentation.')


    new_coords_df = pd.DataFrame([list(coords_df[chrom_column]), new_starts, new_ends, list(coords_df[strand_column])]).transpose()
    new_coords_df.columns = filter_cols
    new_coords_df = pd.concat([new_coords_df, coords_df[add_cols]], axis = 1)
    return(new_coords_df)




def import_modisco_seqlets_txt(seqlets_txt_path):
    import pandas as pd
    from Bio import SeqIO
    """
    Purpose: Based on `basepairmodels` modisco seqlet mapping import the seqlets.txt and write to pd.df
    Output: pd.df of coordinates
    """
    seqs = []
    names = []

    #Parse through .fastq-like format
    with open(seqlets_txt_path) as seqlets:
        for s in SeqIO.parse(seqlets, "fasta"):
            seqs.append(str(s.seq))
            names.append(str(s.id))

    #Write to pd.df and add metadata columns
    seqlets_df = pd.DataFrame([names, seqs]).transpose()
    seqlets_df.columns = ['id', 'seq']
    seqlets_df['region_id'] = [int(s.split(':')[0].replace('example', '')) for s in list(seqlets_df['id'])]
    seqlets_df['modisco_motif_start'] = [int(s.split(':')[1].split('-')[0]) for s in list(seqlets_df['id'])]
    seqlets_df['modisco_motif_end'] = [int(s.split(':')[1].split('-')[1]) for s in list(seqlets_df['id'])]

    return(seqlets_df)

def extract_seqs_from_df(coords_df, fasta_path, chrom_column = 'chrom', start_column = 'start', end_column = 'end'):
    """
    Purpose: Given a pd.df of coordinates (chrom, start, end), defined by the inputs, return genomic sequences.
    Output: list of sequences
    """
    import pysam
    fasta_ref = pysam.FastaFile(fasta_path)
    seqs = []
    for i,row in coords_df.iterrows():
        seq = fasta_ref.fetch(row[chrom_column], row[start_column], row[end_column]).upper()
        seqs.append(seq)
    return(seqs)

def assign_strand_to_coordinates(coords_df, fasta_path, chrom_column = 'chrom', start_column = 'start', end_column = 'end'):
    """
    Purpose: Assign strandiness from a seq of coordinates by checking with the genomic sequences
    Input: coords_df with (chrom, start, end) columns
    Output: list of strandiness of each coordinate
    """
    from Bio.Seq import Seq

    #Extract sequences from genome
    motif_seqs = extract_seqs_from_df(coords_df = coords_df, fasta_path = fasta_path, chrom_column = chrom_column,
                                      start_column = start_column, end_column = end_column)

    #Check sequences for match and strandiness
    match = []
    mismatch_idx = []
    extracted_mismatch_seq = []
    basepairmodel_mismatch_seq = []
    for i in range(coords_df.shape[0]):

        #Define sequences
        coord_seq = Seq(coords_df.loc[i, 'seq']).upper()
        coord_seq_rev = coord_seq.reverse_complement().upper()
        extract_seq = motif_seqs[i]
        if coord_seq==extract_seq:
            match_status = '+'
        elif coord_seq_rev== extract_seq:
            match_status = '-'
        else:
            match_status = '*'
            # mismatch_idx.append(i)
            # extracted_mismatch_seq.append(extract_seq)
            # basepairmodel_mismatch_seq.append(coord_seq)
            print(f'idx {i}: extracted {extract_seq} != basepairmodel {coord_seq}')
        match.append(match_status)

    #assert '*' not in set(match), 'One of the sequences does not match. Check modisco<->.bed relationship.'

    return(match)


def map_modisco_coordinates_to_genome(seqlets_df, input_peak_path, fasta_path, input_seqlen, output_seqlen,
                                      assign_strand = True, existing_strand_col = 'strand',
                                      region_id_col = 'region_id', window_motif_start_col = 'modisco_motif_start', window_motif_end_col = 'modisco_motif_end',
                                      peak_col_names = ['chrom', 'st', 'end', 'name', 'score', 'strand', 'signalValue', 'p', 'q', 'summit', 'start']):
    """
    Purpose: Given a set of coords_df from seqlets.txt, map the genomic coordinates
    Inputs:
        + seqlets_df: pd.df of motif coordinates mapped across a WINDOW (requires chrom, [[window_motif_start_col]], [[window_motif_end_col]], [region_id_col])
        + input_peak_path: path to .narrowpeak (9 columns) or other coordinate set defining the window coordinates that the motifs were derived from
        + fasta_path: path to .fasta file matching genome of motifs
        + input_seqlen: input sequence length of model
        + output_seqlen: output sequence length of motif
        + assign_strand: If True (default), assign strandiness based on `seq` column in seqlets_df object.
        + existing_strand_col: If assign_strand==False, then specify which column contains strand information of motifs
        + region_id_col: Name of column that `seqlets_df` and `input_peak_path.index` share to link genomic coordinates to an index.
        + window_motif_start_col: Name of column that `seqlets_df` contain to designate window-based start
        + window_motif_end_col: Name of column that `seqlets_df` contain to designate window-based end
        + peak_col_names: Columns of the `input_peak_path` .narrowPeak. Defaults to the output 10-column (NOT NARROWPEAK) of shap_scores.peaks_valid_scores.bed output.
    Output: pd.df of motif coordinates from TF-MoDISco

    Note: input_peak_path was designed to be the coordinate output of 'shap_scores' with the extra column, defined by default peak_col_names
    """
    import sys
    sys.path.insert(0, f'/n/projects/mw2098/shared_code/basepairmodels/functions')
    from motifs import import_modisco_seqlets_txt, assign_strand_to_coordinates

    #Import peaks
    peaks_df = pd.read_csv(input_peak_path, sep='\t', header=None, names=peak_col_names)
    assert list(set(list(peaks_df.end - peaks_df.start)))==[input_seqlen], 'Input peak coordinates do not match `input_seqlen`'
    peaks_df[region_id_col] = peaks_df.index
    peaks_df['genomic_window_start'] = peaks_df.start + (input_seqlen - output_seqlen)//2
    peaks_df['genomic_window_end'] = peaks_df.end - (input_seqlen - output_seqlen)//2

    #Format seqlets such that they don't have redundant coordinate columns with peaks_df
    seqlets_df = seqlets_df.rename(columns={window_motif_start_col: 'modisco_window_start', window_motif_end_col: 'modisco_window_end'})

    #Merge modisco coordinates with peak coordinates
    coords_df = seqlets_df.merge(peaks_df, on=region_id_col, how='left').reset_index(drop = True)

    #Resolve window-based and genomic-based coordinates
    coords_df['genomic_motif_start'] = coords_df['genomic_window_start'] + coords_df['modisco_window_start']
    coords_df['genomic_motif_end'] =   coords_df['genomic_window_start'] + coords_df['modisco_window_end']

    #Confirm strandiness
    if assign_strand:
        coords_df['genomic_motif_strand'] = assign_strand_to_coordinates(coords_df = coords_df, fasta_path = fasta_path,
                                                           start_column = 'genomic_motif_start', end_column = 'genomic_motif_end')
    else:
        coords_df['genomic_motif_strand'] = coords_df[existing_strand_col]

    #Extract meaningful columns
    motifs_df = coords_df[['chrom','genomic_motif_start', 'genomic_motif_end', 'genomic_motif_strand', region_id_col, 'genomic_window_start', 'genomic_window_end', 'modisco_window_start', 'modisco_window_end']]
    return(motifs_df)


"""
Palindromic motifs will overlap on the + and - strand, but not necessarily be the exact same coordinates.
Motifs with the same
Keep the motif with a higher contribution score.
Given:
    dfi: list of motifs with CWM-scanned convention from `bpnet cwm-scan` [requires `row_idx`, `pattern_start`, `pattern_end` and `pattern_name`]
    pattern: string of pattern to filter for under `pattern_name` column of dfi
    overlap_threshold: %age of overlap that returns an unacceptable motif.
Returns: list of `row_idx` values to remove that are redundant.
"""

def mark_overlapping_motifs(dfi, pattern, motif_group_column,
                            overlap_threshold = .7,
                            motif_index_column = 'index',
                            motif_start_column = 'start', motif_end_column = 'end',
                            motif_name_column = 'name', motif_length_column = 'width',
                            motif_signal_column = 'agg_sim'):
    import itertools
    remove_row_idx = np.array([])

    def getOverlap(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    i = dfi[dfi[motif_name_column]==pattern]


    #for each region
    for idx in i[motif_group_column].unique():
        i_across_r = i[i[motif_group_column]==idx]

        #Define windows to check for overlaps
        window_interv_dict = {row[motif_index_column]: (row[motif_start_column],row[motif_end_column])
                              for idx, row in i_across_r.iterrows()}

        #Define all combinations of windows to test
        overlap_combs = list(itertools.combinations(window_interv_dict.keys(),2))

        #Test each combination, mark the one with a lower contribution score for removal
        for comb in overlap_combs:
            ov_count = getOverlap(window_interv_dict[comb[0]], window_interv_dict[comb[1]])
            if ov_count/i[motif_length_column].iloc[0]: #If there is sufficient, overlap, select comb to remove.
                if i[motif_signal_column][i[motif_index_column]==comb[0]].iloc[0] > i[motif_signal_column][i[motif_index_column]==comb[1]].iloc[0]:
                    remove_row_idx = np.append(remove_row_idx, comb[1])
                else:
                    remove_row_idx = np.append(remove_row_idx, int(comb[0]))

    return remove_row_idx


def write_fasta(seqs, fasta_file, wrap=80):
    """Write sequences to a fasta file.

    Parameters
    ----------
    seqs : list of sequences
        Sequences indexed by sequence id.
    fasta_file : str
        Path to write the sequences to.
    wrap: int
        Number of AA/NT before the line is wrapped.
    """
    with open(fasta_file, 'w') as f:
        for gid, gseq in enumerate(seqs):
            f.write('>{}\n'.format(gid))
            for i in range(0, len(gseq), wrap):
                f.write('{}\n'.format(gseq[i:i + wrap]))
