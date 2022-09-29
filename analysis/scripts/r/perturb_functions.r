#Melanie Weilert
#July 2020
#Purpose: Functions pertaining to perturbation collections.

########################################################################
#Function defining any unique motif pair set.

# Inputs:
#   name_x: string that matches a motif `pattern_name` in the `dfi` given. Will be treated as motif_x in a motif pair.
#   name_y: string that matches a motif `pattern_name` in the `dfi` given. Will be treated as motif_y in a motif pair.
#   dfi: data frame of motifs to find pairs across. Should be in the format of `bpnet cwm-scan` command output.
#   remove_overlapping: If T, remove all motif pairs where the motifs would overlap
#   cores: Number of cores to run when removing overlapping motifs
# Outputs:
#   pairs.df: data.frame object that contains motif information of each motif pair that shares an `example_idx` window.

########################################################################

find_motif_pairs<-function(name_x, name_y, dfi, 
                           merge_by = "example_idx",
                           motif_column = 'pattern_name',
                           motif_unique_column = 'pattern_name_unique',
                           motif_start_column = 'pattern_start',
                           motif_end_column = 'pattern_end',
                           remove_overlapping = T, cores = 8){
  
  #subset by motifs
  dfi_x<-dfi %>% dplyr::filter(!!rlang::sym(motif_column)==name_x)
  dfi_y<-dfi %>% dplyr::filter(!!rlang::sym(motif_column)==name_y)
  
  #join unique motif pairs IF they share the same `example_idx`
  motif_unique_column_x<-paste0(motif_unique_column, '.x')
  motif_unique_column_y<-paste0(motif_unique_column, '.y')
  motif_column_x<-paste0(motif_column, '.x')
  motif_column_y<-paste0(motif_column, '.y') 
  motif_start_column_x<-paste0(motif_start_column, '.x')
  motif_start_column_y<-paste0(motif_start_column, '.y') 
  motif_end_column_x<-paste0(motif_end_column, '.x')
  motif_end_column_y<-paste0(motif_end_column, '.y') 
  
  pairs_xy.df<-dplyr::inner_join(x = dfi_x, y = dfi_y, by = merge_by)%>%
    dplyr::filter( !!rlang::sym(motif_unique_column_x) != !!rlang::sym(motif_unique_column_y))
  
  if(nrow(pairs_xy.df)==0){
    return(NULL)
  }else{
    #add a motif pair identifier
    pairs_xy.df$pair_name<-paste0(pairs_xy.df[[motif_column_x]], "_", pairs_xy.df[[motif_column_y]])
    
    #Pairwise comparisons require x->y and y->x, account for this for motifs that are different
    if(name_x==name_y){pairs.df<-pairs_xy.df} else{
      pairs_yx.df<-dplyr::inner_join(x = dfi_y, y = dfi_x, by = merge_by)%>%
        dplyr::filter( !!rlang::sym(motif_unique_column_x) != !!rlang::sym(motif_unique_column_y))
      pairs_yx.df$pair_name<-paste0(pairs_yx.df[[motif_column_x]], "_", pairs_yx.df[[motif_column_y]])
      pairs.df<-rbind(pairs_xy.df, pairs_yx.df)
    }
    
    if(remove_overlapping){
      #Check if motif pairs are overlapping.
      is_overlapping<-mclapply(1:nrow(pairs.df), function(x){
        any(pairs.df[[motif_start_column_x]][x]:pairs.df[[motif_end_column_x]][x] %in% pairs.df[[motif_start_column_y]][x]:pairs.df[[motif_end_column_y]][x])
      }, mc.cores = cores) %>% unlist
      message(sum(is_overlapping), " overlapping motifs removed. ", sum(!is_overlapping), " motifs remain.")
      pairs.df<-pairs.df[!is_overlapping,]
    }
    
    pairs.df$pair_id<-1:nrow(pairs.df)
    return(pairs.df)
  }
  
  
}

########################################################################
#Extracting perturbations from generate_perturbations***.py and saving to data.frame

# For any unique pair set, conduct this code below to extract pairwise perturbations. In the event of a motif A and motif B measured by TF A, we will have the following metrics:
#
# 1. WT (dA_ref): height/sum of TF A [at motif A] when motif A and B are present
# 2. dA_b: height/sum of TF A [at motif A] when motif A is present.
# 3. dA_a: height/sum of TF A [at motif A] when motif B is present
# 4. dA_ab: height/sum of TF A [at motif A] when no motif is present.
#
# [at motif A] is only pertinent in the case that the measurements are found across the motif-window.

# Inputs:
#   pairs.df: data.frame from the `find_motif_pairs` function above.
#   motif_specs.list: list of format `motif_specs.list=list(motif: list(color = ?, tf = ?))`. This will map which motif/task association is meaningful.
#   perturb_prefix: filepath prefix showing where the chromosome-separated perturbation.tsv.gz files are
#     + note that if you import the perturbations pre-emptively and this input is a list, then it will automatically use the list (separated by chromosomes)
#   cores: Number of cores to run when extracting perturbations across motif pairs
# Outputs:
#   summary.df: data.frame object that contains pairwise perturbations across each motif pair. The indexes of the motif_pair should match the summary.df.
########################################################################

#The function below seeks to extract the perturbation valued needed to fulfill the whole-window measurements.
#The perturb_prefix is based off of chromosome-separated perturbations that were derived using the
#code from `/n/projects/mw2098/shared_code/bpnet/bpnet_generate_whole_window_perturbations.py`.

extract_whole_window_paired_perturbs<-function(pairs.df, motif_specs.list, perturb_prefix, 
                                               chrom_column = 'example_chrom.x', 
                                               region_index_column = 'example_idx',
                                               motifA_column = 'pattern_name.x',
                                               motifA_unique_column = 'pattern_name_unique.x', 
                                               motifA_index_column = 'motif_id.x',
                                               motifA_center_column = 'pattern_center.x',
                                               motifB_column = 'pattern_name.y',
                                               motifB_unique_column = 'pattern_name_unique.y', 
                                               motifB_index_column = 'motif_id.y',
                                               motifB_center_column = 'pattern_center.y',
                                               cores = 8){

  #For each chromosome, find motif pairs and subset
  summary.df<-lapply(unique(pairs.df[[chrom_column]]), function(chr){
    message("Processing perturbations from: ", chr)
    #Subset inputs by chromosome
    if(type(perturb_prefix)=='list'){
      perturbs_across_chr.df<-perturb_prefix[[chr]]
    } else if(type(perturb_prefix)=='character'){
      perturbs_across_chr.df<-read_tsv(paste0(perturb_prefix, chr, ".tsv.gz")) %>% as.data.frame
    }
    pairs_across_chr.df<-pairs.df %>% dplyr::filter(!!rlang::sym(chrom_column) == chr)
    
    summary_across_chr.df<-mclapply(pairs_across_chr.df$pair_id, function(pid){
      pair.df<-pairs_across_chr.df %>% dplyr::filter(pair_id==pid)
      
      #Define different parameters for subsetting the perturbs
      umotifA<-pair.df[[motifA_unique_column]]
      motifA<-pair.df[[motifA_column]]
      TFA<-motif_specs.list[[motifA]]$tf
      TFA_cols<-paste0(TFA, c("/pc", "/pred_sum"))
      
      umotifB<-pair.df[[motifB_unique_column]]
      motifB<-pair.df[[motifB_column]]

      #Height of TFA at motifA when both motifs are present
      WT.df<-perturbs_across_chr.df %>% dplyr::filter(!!rlang::sym(region_index_column) == pair.df[[region_index_column]], mut=="Reference") %>%
        dplyr::select(TFA_cols)
      colnames(WT.df)<-gsub(TFA, "WT", colnames(WT.df))
      
      #Height of TFA at motifA when only motif A is present
      dA_b.df <- perturbs_across_chr.df %>% dplyr::filter(!!rlang::sym(region_index_column) == pair.df[[region_index_column]], mut==umotifB) %>%
        dplyr::select(TFA_cols)
      colnames(dA_b.df)<-gsub(TFA, "dA_b", colnames(dA_b.df))
      
      #Height of TFA at motifA when only motif B is present
      dA_a.df <- perturbs_across_chr.df %>% dplyr::filter(!!rlang::sym(region_index_column) == pair.df[[region_index_column]], mut==umotifA) %>%
        dplyr::select(TFA_cols)
      colnames(dA_a.df)<-gsub(TFA, "dA_a", colnames(dA_a.df))
      
      #Height of TFA at motifA when neither motif is present
      dA_ab.df <- perturbs_across_chr.df %>% dplyr::filter(!!rlang::sym(region_index_column) == pair.df[[region_index_column]], 
                                                           (mut==paste0(umotifA, "_", umotifB) | mut==paste0(umotifB, "_", umotifA))) %>%
        dplyr::select(TFA_cols)
      colnames(dA_ab.df)<-gsub(TFA, "dA_ab", colnames(dA_ab.df))
      
      df<-data.frame(pair_id=pid, motifA = motifA, umotifA = umotifA, motifB = motifB, umotifB = umotifB)
      df<-cbind(df, WT.df, dA_b.df, dA_a.df, dA_ab.df)
      
      #Record additional helpful information
      df$row_idx_A<-pair.df[[motifA_index_column]]
      df$row_idx_B<-pair.df[[motifB_index_column]]
      df$strand_A<-pair.df$strand.x
      df$strand_B<-pair.df$strand.y
      df$distance_b_minus_a<-pair.df[[motifB_center_column]] - pair.df[[motifA_center_column]]

      return(df)
    }, mc.cores = cores) %>% rbindlist
  }) %>% rbindlist

  return(summary.df)
}


#The function below seeks to extract the paired perturbation valued needed to fulfill the motif-window measurements.
#The perturb_prefix is based off of chromosome-separated perturbations that were derived using the
#code from `/n/projects/mw2098/shared_code/bpnet/bpnet_generate_motif_window_perturbations.py`.

extract_motif_window_paired_perturbs<-function(pairs.df, motif_specs.list, perturb_prefix, 
                                               chrom_column = 'example_chrom.x', 
                                               region_index_column = 'example_idx',
                                               motifA_column = 'pattern_name.x',
                                               motifA_unique_column = 'pattern_name_unique.x', 
                                               motifA_index_column = 'motif_id.x',
                                               motifA_center_column = 'pattern_center.x',
                                               motifB_column = 'pattern_name.y',
                                               motifB_unique_column = 'pattern_name_unique.y', 
                                               motifB_index_column = 'motif_id.y',
                                               motifB_center_column = 'pattern_center.y',
                                               cores = 8){

  #For each chromosome, find motif pairs and subset
  summary.df<-lapply(unique(pairs.df[[chrom_column]]), function(chr){
    message("Processing perturbations from: ", chr)
    #Subset inputs by chromosome
    if(type(perturb_prefix)=='list'){
      perturbs_across_chr.df<-perturb_prefix[[chr]]
    } else if(type(perturb_prefix)=='character'){
      perturbs_across_chr.df<-read_tsv(paste0(perturb_prefix, chr, ".tsv.gz")) %>% as.data.frame
    }
    pairs_across_chr.df<-pairs.df %>% dplyr::filter(!!rlang::sym(chrom_column) == chr)

    summary_across_chr.df<-mclapply(pairs_across_chr.df$pair_id, function(pid){
      pair.df<-pairs_across_chr.df %>% dplyr::filter(pair_id==pid)

      #Define different parameters for subsetting the perturbs
      umotifA<-pair.df[[motifA_unique_column]]
      motifA<-pair.df[[motifA_column]]
      TFA<-motif_specs.list[[motifA]]$tf
      TFA_cols<-paste0(TFA, c("/pc", "/pred_sum", "/pred_max"))

      umotifB<-pair.df[[motifB_unique_column]]
      motifB<-pair.df[[motifB_column]]

      #Height of TFA at motifA when both motifs are present
      WT.df<-perturbs_across_chr.df %>% dplyr::filter(!!rlang::sym(region_index_column) == pair.df[[region_index_column]], 
                                                      motif==umotifA, mut=="Reference") %>%
        dplyr::select(TFA_cols)
      colnames(WT.df)<-gsub(TFA, "WT", colnames(WT.df))

      #Height of TFA at motifA when only motif A is present
      dA_b.df <- perturbs_across_chr.df %>% dplyr::filter(!!rlang::sym(region_index_column) == pair.df[[region_index_column]], 
                                                          motif==umotifA, mut==umotifB) %>%
        dplyr::select(TFA_cols)
      colnames(dA_b.df)<-gsub(TFA, "dA_b", colnames(dA_b.df))

      #Height of TFA at motifA when only motif B is present
      dA_a.df <- perturbs_across_chr.df %>% dplyr::filter(!!rlang::sym(region_index_column) == pair.df[[region_index_column]], 
                                                          motif==umotifA, mut==umotifA) %>%
        dplyr::select(TFA_cols)
      colnames(dA_a.df)<-gsub(TFA, "dA_a", colnames(dA_a.df))

      #Height of TFA at motifA when neither motif is present
      dA_ab.df <- perturbs_across_chr.df %>% dplyr::filter(!!rlang::sym(region_index_column) == pair.df[[region_index_column]], 
                                                           motif==umotifA,
                                                           (mut==paste0(umotifA, "_", umotifB) | mut==paste0(umotifB, "_", umotifA))) %>%
        dplyr::select(TFA_cols)
      colnames(dA_ab.df)<-gsub(TFA, "dA_ab", colnames(dA_ab.df))

      df<-data.frame(pair_id=pid, motifA = motifA, umotifA = umotifA, motifB = motifB, umotifB = umotifB)
      df<-cbind(df, WT.df, dA_b.df, dA_a.df, dA_ab.df)

      #Record additional helpful information
      df$row_idx_A<-pair.df[[motifA_index_column]]
      df$row_idx_B<-pair.df[[motifB_index_column]]
      df$strand_A<-pair.df$strand.x
      df$strand_B<-pair.df$strand.y
      df$distance_b_minus_a<-pair.df[[motifB_center_column]] - pair.df[[motifA_center_column]]

      return(df)
    }, mc.cores = cores) %>% rbindlist
  }) %>% rbindlist

  return(summary.df)
}

#The function below seeks to extract the single perturbation valued needed to fulfill the motif-window measurements.
#The perturb_prefix is based off of chromosome-separated perturbations that were derived using the
#code from `/n/projects/mw2098/shared_code/bpnet/bpnet_generate_motif_window_perturbations.py`.

extract_motif_window_single_perturbs<-function(instances.df, motif_specs.list, perturb_prefix, cores = 8){
  
  message('Not updated for basepairmodel motifs.df column formatting! TODO!')
  
  #For each chromosome, find motif pairs and subset
  summary.df<-lapply(unique(instances.df$example_chrom), function(chr){
    message("Processing perturbations from: ", chr)
    #Subset inputs by chromosome
    perturbs_across_chr.df<-read_tsv(paste0(perturb_prefix, chr, ".tsv.gz")) %>% as.data.frame
    inst_across_chr.df<-instances.df %>% dplyr::filter(example_chrom == chr)

    summary_across_chr.df<-mclapply(inst_across_chr.df$row_idx, function(ridx){
      inst.df<-inst_across_chr.df %>% dplyr::filter(row_idx==ridx)

      #Define different parameters for subsetting the perturbs
      umotifA<-inst.df$pattern_name_unique
      motifA<-inst.df$pattern_name
      TFA<-motif_specs.list[[motifA]]$tf
      TFA_cols<-paste0(TFA, c("/pc", "/pred_sum", "/pred_max"))

      #Height of TFA at motifA when both motifs are present
      WT.df<-perturbs_across_chr.df %>% dplyr::filter(example_idx == inst.df$example_idx, motif==umotifA, mut=="Reference") %>%
        dplyr::select(TFA_cols)
      colnames(WT.df)<-gsub(TFA, "WT", colnames(WT.df))

      #Height of TFA at motifA when only motif B is present
      dA_a.df <- perturbs_across_chr.df %>% dplyr::filter(example_idx == inst.df$example_idx, motif==umotifA, mut==umotifA) %>%
        dplyr::select(TFA_cols)
      colnames(dA_a.df)<-gsub(TFA, "dA_a", colnames(dA_a.df))

      df<-data.frame(row_idx=ridx, motifA = motifA, umotifA = umotifA)
      df<-cbind(df, WT.df, dA_a.df)

      #Record additional helpful information
      df$strand_A<-inst.df$strand
      return(df)
    }, mc.cores = cores) %>% rbindlist
  }) %>% rbindlist

  return(summary.df)
}
