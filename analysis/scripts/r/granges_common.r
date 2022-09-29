library(GenomicRanges)
library(rtracklayer)
library(dplyr)

#20180525_modification: See commented out code with Views(). Upgrading to R version 5.0 will cause RangesList to become deprecated and GRangesList is not compatible with this function. In order to circumvent, a newer updated version of splitting a GRanges into GRangesList was added. 
#20180625_modification: added check_chromosome_boundary function in order to determine which peaks interfere with the chromosome boundaries of the assemblies you want to obtain a metapeak from. This will help you with exo_metapeak_matrix function errors.

apply_seqlengths <- function(gr, genome=Mmusculus) {
  seqlengths(gr) <- seqlengths(genome)[seqlevels(gr)]
  gr
}

check_coverage_argument <- function(cvg, regions=NULL) {
  if(class(cvg) == "character") {
    if(is.null(regions)) {
      cvg <- rtracklayer::import(cvg, as="RleList")
    } else {
      stopifnot(file.exists(cvg))
      cvg <- rtracklayer::import(cvg, which=regions, as="RleList")
    }
  } 
  cvg
}

regionApply <- function(regions, cvg, func, ...) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(unlist(
           viewApply(
             #Views(cvg, as(regions, "GRangesList")),
             Views(cvg, split(regions, seqnames(regions))),
             function(x) { func(as.numeric(x)) },
             simplify=FALSE
           ), use.names=FALSE), use.names=FALSE)
  ans
}

regionApply_list <- function(regions, cvg, func) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewApply(
             #Views(cvg, as(regions, "GRangesList")),
             Views(cvg, split(regions, seqnames(regions))),
             function(x) { func(as.numeric(x)) },
             simplify=FALSE
           ), use.names=FALSE)
  ans
}

regionSums <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewSums(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMeans <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMeans(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionWhichMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewWhichMaxs(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionWhichMins <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
    viewWhichMins(
      #Views(cvg, as(regions, "GRangesList"))
      Views(cvg, split(regions, seqnames(regions)))
    ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMaxs <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMaxs(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

regionMins <- function(regions, cvg) {
  cvg <- check_coverage_argument(cvg, regions)
  seqlevels(regions) <- names(cvg)
  oo <- order(regions)
  regions <- regions[oo]
  ans <- unlist(
           viewMins(
             #Views(cvg, as(regions, "GRangesList"))
             Views(cvg, split(regions, seqnames(regions)))
           ), use.names=FALSE)
  ans[oo] <- ans  # restore original order
  ans
}

total_signal <- function(cov) {
  if(class(cov) == "character") {
    cache.file <- paste0(cov, ".ts.rds")
    if(file.exists(cache.file)) {
      return(readRDS(cache.file))
    } else {
      cov <- check_coverage_argument(cov)
      ts <- sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x)))))
      saveRDS(ts, file=cache.file)
      return(ts)
    }
  } else {
    return(sum(as.numeric(sapply(cov, function(x) sum(as.numeric(x))))))
  }
}

nexus_regionSums <- function(gr, sample){
  signal <- regionSums(gr, sample$pos) + abs(regionSums(gr, sample$neg))
  signal
}

#FXN: Check Chromosome Boundaries

#Inputs: (1) gr=regions of interest for exo_metapeak_matrix, (2)) resize_boundary=the extensions created by the "downstream" input command and "upstream" input command on exo_metapeak_matrix.

#Outputs: Returns an vector of indices that oversteps the designated chromosome boundaries and are therefore incompatible with matrix computation. Use indices WRT "gr" input.

#Example: Input: Granges called "example.ranges" of 100 regions from the dm6 genome that we want to plot +/- 500bp from the center.
#         Use Case: library(BSgenome.Dmelanogaster.UCSC.dm6)  <-call library of genome you want to check across
#                   bad_spots<-check_chromosome_boundaries(gr=example.ranges, genome=BSgenome.Dmelanogaster.UCSC.dm6, resize_boundary=1001) <-apply the function to find if you have any bad regions
#                   if(length(bad_spots>0){example.ranges<-example.ranges[-bad_spots]} <- IF the bad regions exist, remove them from your GRanges. ****Make sure to always wrap this in an IF statement. Calling a variable of length=0 as a negative will erase the entire contents of the variable.

check_chromosome_boundaries<-function(gr, resize_boundary, genome=NULL){
  
  if(!is.null(genome)){
    seqlevels(gr)<-seqlevels(genome)
    seqinfo(gr)<-seqinfo(genome)
  }
  
  key<-data.frame(seqlength_filler=seqlengths(gr), seqnames=names(seqlengths(gr)))
  gr_merged<-merge(as.data.frame(gr), key, by.x="seqnames", by.y="seqnames") %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  gr_merged$end_distance<-gr_merged$seqlength_filler-(end(gr_merged)+resize_boundary) #distance from end
  gr_merged$start_distance<-start(gr_merged)-resize_boundary #distance from start
  
  flagged_indexes<-which(gr_merged$end_distance<=0 | gr_merged$start_distance<=0)
  if(length(flagged_indexes)>0){
    gr_merged<-gr_merged[-flagged_indexes]
  }
  message("Returning GRanges with metadata. We removed ", length(flagged_indexes), " regions.")
  return(gr_merged)
}

#FXN: Augment GRanges

#Purpose: Intended to provide augmentation for BPNet training of a GRanges. Augmentation is when you offset and distort your data in order for the network to get a more robust sense of its scope. In this case, translation across a genome is the only augmentation this supports. 

#Inputs: 
#(1) gr=regions of interest, 
#(2) augmentation_repeat_limits = the range of repeated regions for variable augmentation. If you desire only 1 region ever, please select a scalar input of 1.
#(3) augmentation_window_limits = the range of offset centers that you translate the windows across.
#(4) desired_window_width = how wide you would like the GRanges extended to post-augmentation for BPNet training.
#(5) output_bed_file = if not NA, then writes GRanges to an output BED file at the specified path.
#(6) seed = if not NA, then selects seed to control randomization for reproducibility

#Outputs: Returns a GRanges from the function with the desired augmentation applied. If output_bed_file is not NA, also write this GRanges to specified BED path.

augment_granges<-function(gr, augmentation_repeat_limits=1:5, augmentation_window_limits=-100:100, desired_window_width=1000, output_bed_file=NA, seed=NA){
  #set seed if desired
  if(!is.na(seed)){ 
    set.seed(seed)
  }
  
  #Determine augmentation repeats
  augmentation_repeats<-sample(augmentation_repeat_limits, length(gr), replace=T)
  gr_idx<-rep(1:length(gr), augmentation_repeats)
  gr_exp<-gr[gr_idx] %>% resize(width=1, fix="start")
  gr_aug<-gr_exp %>% IRanges::shift(shift=sample(augmentation_window_limits, length(gr_exp), replace=T))
  
  #Resize to desired window width
  gr_final<-gr_aug %>% resize(width=desired_window_width, fix="center")
  
  #write to BED if desired
  if(!is.na(output_bed_file)){ 
    export(gr_final, output_bed_file)
    return(gr_final)
  }
  else{return(gr_final)}
}



