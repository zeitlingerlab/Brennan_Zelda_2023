#Melanie Weilert
#August 2020
#Purpose: Functions pertaining motifs surrounding BPNet.

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(plyranges)
library(testit)
library(readr)

########################################################################
# Input:
#   + `dfi` can either be a GRanges or a data.frame
#     + `dfi` is the format of the output of `bpnet cwm-scan`. It requires a `pattern_name` column for filtering
#   + `pattern_name` is the name of the motif to filter for
# Output: GRangesObject with same columns as `dfi.df` with the palindromic motifs filtered.
########################################################################

remove_palindromic_motifs_from_bpnet_instances<-function(dfi, pattern_to_filter, 
                                                         chrom_name = 'example_chrom', 
                                                         start_name = 'pattern_start_abs', 
                                                         end_name = 'pattern_end_abs',
                                                         value_name = 'contrib_weighted',
                                                         motif_name = 'pattern_name',
                                                         region_name = 'example_idx',
                                                         starts.in.df.are.0based = F){

  if(class(dfi)[1]=="GRanges"){
    dfi.gr<-dfi
  } else{
    #Convert to GRanges
    dfi.gr<-makeGRangesFromDataFrame(dfi, keep.extra.columns = T, starts.in.df.are.0based = starts.in.df.are.0based,
                                     seqnames.field = chrom_name, start.field = start_name, end.field = end_name)
  }

  #Subset by motifs
  dfi_filt.gr<-dfi.gr %>% dplyr::filter(!!rlang::sym(motif_name)==pattern_to_filter)

  #Create islands of overlapping motifs
  dfi_islands.gr<-GenomicRanges::reduce(dfi_filt.gr, ignore.strand = T)
  dfi_islands.gr$island_idx<-1:length(dfi_islands.gr)

  #Remap island ids to the motifs
  dfi_vs_islands.ov<-findOverlaps(dfi_filt.gr, dfi_islands.gr, ignore.strand = T)
  testit::assert("Motifs overlap redundantly...", (dfi_vs_islands.ov@from %>% unique() %>% length) == (dfi_vs_islands.ov@from %>% length))
  dfi_filt.gr$island_idx<-dfi_islands.gr$island_idx[dfi_vs_islands.ov@to]

  #For each island id, select motif with the best importance
  dfi_unique.gr <- dfi_filt.gr %>% as.data.frame %>%
    dplyr::group_by(island_idx, !!rlang::sym(region_name)) %>% #Group by both island_id and example_idx to keep the same motifs mapped across different task-annotated windows
    dplyr::slice_max(order_by = !!rlang::sym(value_name), n = 1) %>% as.data.frame() %>%
    makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = F) %>%
    plyranges::select(-island_idx)

  message(length(dfi_filt.gr)-length(dfi_unique.gr), " palindromic motifs removed of ", length(dfi_filt.gr))
  return(dfi_unique.gr)
}



########################################################################
# Input:
#   + `dfi` can either be a GRanges or a data.frame
#     + `dfi` is the format of the output of `bpnet cwm-scan`. It requires a `pattern_name` column for filtering
# Output: GRangesObject with same columns as `dfi.df` with the palindromic motifs filtered.
########################################################################

remove_redundant_motifs_from_bpnet_instances<-function(dfi,
                                                       chrom_name = 'example_chrom', 
                                                       start_name = 'pattern_start_abs', 
                                                       end_name = 'pattern_end_abs',
                                                       value_name = 'contrib_weighted',
                                                       region_name = 'example_idx',
                                                       starts.in.df.are.0based = F){
  if(class(dfi)[1]=="GRanges"){
    dfi.gr<-dfi
  } else{
    #Convert to GRanges
    dfi.gr<-makeGRangesFromDataFrame(dfi, keep.extra.columns = T, starts.in.df.are.0based = starts.in.df.are.0based,
                                     seqnames.field = chrom_name, start.field = start_name, end.field = end_name)
  }
  
  #Subset by motifs
  dfi_filt.gr<-dfi.gr
  
  #Create islands of overlapping motifs
  dfi_islands.gr<-GenomicRanges::reduce(dfi_filt.gr, ignore.strand = T)
  dfi_islands.gr$island_idx<-1:length(dfi_islands.gr)
  
  #Remap island ids to the motifs
  dfi_vs_islands.ov<-findOverlaps(dfi_filt.gr, dfi_islands.gr, ignore.strand = T)
  testit::assert("Motifs overlap redundantly...", (dfi_vs_islands.ov@from %>% unique() %>% length) == (dfi_vs_islands.ov@from %>% length))
  dfi_filt.gr$island_idx<-dfi_islands.gr$island_idx[dfi_vs_islands.ov@to]
  
  #For each island id, select motif with the best importance
  dfi_unique.gr <- dfi_filt.gr %>% as.data.frame %>%
    dplyr::group_by(island_idx, !!rlang::sym(region_name)) %>% #Group by both island_id and example_idx to keep the same motifs mapped across different task-annotated windows
    dplyr::slice_max(order_by = !!rlang::sym(value_name), n = 1) %>% as.data.frame() %>%
    makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = F) %>%
    plyranges::select(-island_idx)
  
  message(length(dfi_filt.gr)-length(dfi_unique.gr), " redundant motifs removed of ", length(dfi_filt.gr))
  return(dfi_unique.gr)
}

########################################################################
# Purpose: Filter out redundant example_idx windows.
# When BPNet was run, it was trained on peaks for each of the factors in the model. However, some regions of the genome have peaks for multiple TFs, meaning that some `example_idx` windows cover the same regions. This results in redundantly mapped motifs. Here, we will filter out those redundant by passing two criteria in order of priority:
#
# For a set of overlapping windows:
#
# 1. Keep the window with the highest number of motifs mapped (not just patterns of interest)
# 2. If (1) is the same, keep the window with the most central set of motifs in the window.
#
# Input: `instances.df` is the format of the output of `bpnet cwm-scan`. It requires a `pattern_name` column for filtering.
# Output: a vector of `example_idx`s that are unique and prioritized.

# Input:
#   + `dfi.df` is the format of the output of `bpnet cwm-scan`. It requires a `pattern_name` column for filtering
#   + `pattern_name` is the name of the motif to filter for
# Output: Vector of `example_idx` that should be kept.
########################################################################

mark_nonredundant_bpnet_windows<-function(instances, window_length = 1000){

  if(class(instances)=="GRanges"){
    instances.gr<-instances
  } else{
    instances.gr<-makeGRangesFromDataFrame(instances, keep.extra.columns = T, starts.in.df.are.0based = T,
                                           seqnames.field = 'example_chrom', start.field = 'pattern_start_abs', end.field = 'pattern_end_abs')
  }
  #Create GRanges of windows from information
  windows.gr<-GRanges(seqnames(instances.gr), IRanges(start = instances.gr$example_start, end = instances.gr$example_end), strand = "*")

  #Create islands of overlapping windows
  islands.gr<-GenomicRanges::reduce(windows.gr, ignore.strand = T)
  islands.gr$island_idx<-1:length(islands.gr)

  #Remap island ids to the motifs
  instances_vs_islands.ov<-findOverlaps(instances.gr, islands.gr, ignore.strand = T)
  testit::assert("Windows overlap rerundantly...", instances_vs_islands.ov@from %>% unique() %>% length == instances_vs_islands.ov@from %>% length)
  instances.gr$island_idx<-islands.gr$island_idx[instances_vs_islands.ov@to]

  #Get example indexes with the most frequent and most central motifs
  example_idx_filt.df <- instances.gr %>% as.data.frame %>%
    dplyr::group_by(island_idx, example_idx) %>%
    dplyr::summarize(motif_count = n(),   #Get the central score and the motif count
                     central_score = sum(abs(pattern_start-(window_length/2))) + sum(abs(pattern_end-(window_length/2)))) %>%
    dplyr::group_by(island_idx) %>%
    slice_max(order_by = motif_count, n = 1) %>%  #get top motif counts
    slice_max(order_by = min(central_score)) #get top by most centered motifs

  testit::assert("There is not a 1:1 ratio of island:example after filtering",
         example_idx_filt.df$example_idx %>% unique %>% length == example_idx_filt.df$island_idx %>% unique %>% length)

  instances_filt.gr<-instances.gr[instances.gr$example_idx %in% example_idx_filt.df$example_idx]

  return(instances_filt.gr$example_idx %>% unique)
}




########################################################################
# Purpose: Resolve similar motif overlaps in an inter-pattern fashion

# Many of the mapped motifs contain similarities such that they will map over the same sets of coordinate regions. 
# To resolve these, we need to approach each comparison with a hierarchical, individualized regimen. 
# Below, the steps will be listed for each group.
# 
# Here, we define a function for differentiating a set of motifs in a collected GRanges by a column. 
# If in that column: x > y, then we select x as the preferential motif. 
# We can use contribution scores, pattern lengths, etc. as long as the preferential feature is LARGER, then this comparison can be made.
#
# Input: 
# + motif_x: string of `pattern_name` column to subset by for first motif
# + motif_y: string of `pattern_name` column to subset by for second motif
# + dfi.gr: Granges object that contains both motifs and is derived from `bpnet cwm-scan`
# + diff_column: column name to differentiate by. If x>y in this column, then x will be kept and y will be removed.
# + id_column: column name that designated motif unique ids
# + ignore.strand = T
# Output: dfi_filt.gr: a filtered GRanges object with redundant motifs removed

########################################################################

differentiate_overlapping_motifs<-function(motif_x, motif_y, dfi.gr, diff_column, id_column, ignore.strand = T){
  #Subset motifs
  x.gr<-dfi.gr %>% plyranges::filter(pattern_name == motif_x)
  y.gr<-dfi.gr %>% plyranges::filter(pattern_name == motif_y)
  
  #Find overlaps between the motifs
  ov.df<-findOverlaps(x.gr, y.gr, ignore.strand = ignore.strand) %>% as.data.frame()
  
  #Determine which possess a match score that is better
  ov.df$diff_column_x<-x.gr@elementMetadata[[diff_column]][ov.df$queryHits]
  ov.df$diff_column_y<-y.gr@elementMetadata[[diff_column]][ov.df$subjectHits]
  ov.df$is_x_a_match<-ov.df$diff_column_x>ov.df$diff_column_y
  
  #Mark motif IDs that need to be removed
  ids_to_remove<-c(x.gr@elementMetadata[[id_column]][ov.df$queryHits[!ov.df$is_x_a_match]],
                   y.gr@elementMetadata[[id_column]][ov.df$subjectHits[ov.df$is_x_a_match]])
  
  dfi_filt.gr<-dfi.gr[!(dfi.gr@elementMetadata[[id_column]] %in% ids_to_remove)]
  message(length(ids_to_remove), ' regions filtered. ', length(dfi_filt.gr), ' motifs returned.')
  return(dfi_filt.gr)
}


#Kudos to Yue for making this on 02/01/2021
#Given a string sequence 's', change it into a 1-hot encoded matrix
one_hot_encode_DNA <- function(s, n = c("A", "C", "G", "T")) {
    # Construct matrix
    s <- toupper(s)
    seq_len <- nchar(s)
    seq_split <- unlist(strsplit(x = s, split = ""))
    seq_mat <- matrix(data = rep(0, seq_len * length(n)), nrow = 4)
    rownames(seq_mat) <- n
    colnames(seq_mat) <- seq_split
    # Encode
    for (i in n) {
        seq_mat[rownames(seq_mat) == i, colnames(seq_mat) == i] <- 1
    }
    return(seq_mat)
}
