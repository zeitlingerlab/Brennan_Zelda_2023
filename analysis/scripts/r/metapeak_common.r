library(GenomicRanges)
library(reshape2)
library(lattice)
library(magrittr)
library(IRanges)
library(data.table)

#20180525_modification: See commented out code with Views(). Upgrading to R version 5.0 will cause RangesList to become deprecated and GRangesList is not compatible with this function. In order to circumvent, a newer updated version of splitting a GRanges into GRangesList was added.

read_matrix <- function(gr, cov, reverse_reads=FALSE, df=FALSE, nu=25) {
  if(class(cov) == "character") {
    #cov <- import.bw(cov, which=gr, as="RleList") #deprecated command, rtracklayer will recognize values
    cov <- rtracklayer::import.bw(cov, which=gr, as="RleList")
  }
  transform_function <- if(reverse_reads) { rev } else { identity }
  gr <- granges(gr) # remove metadata to speed up RleList data retrieval in Views()
  o <- order(gr)
  gr <- gr[o]

  if(df==TRUE){
	reads.list <- regionApply_list(gr,cov, as.numeric )
  	reads.list <-lapply(reads.list, function(reads) { approx(reads, n=nu)$y })
	reads.m <- matrix(as.numeric(unlist(reads.list, use.name=FALSE)), nrow=length(reads.list), byrow=T)
	if(reverse_reads == T)reads.m <- reads.m[, ncol(reads.m):1]
  }else{
    reads.list <- regionApply(gr, cov, as.numeric)
	  #rl <- as(gr, "GRangesList")
	  rl <- split(gr, seqnames(gr))
	  rl <- rl[which(lapply(rl, function(x)length(x))!=0)] #eliminate chromosomes that are not used
	  view <- RleViewsList(rleList=cov[names(rl)], rangesList=rl)
	  reads.list <- viewApply(view, function(x) { transform_function(as.numeric(x)) })
	  reads.m <- matrix(unlist(sapply(reads.list, as.numeric)), nrow=length(gr), byrow=TRUE)
  }
  reads.m[o, ] <- reads.m
  reads.m
}

standard_metapeak_matrix <- function(regions.gr, sample.cov, upstream=100, downstream=100, keep_region_coordinates = FALSE, diff_length = FALSE, approx_nu=25) {
  if(diff_length == FALSE){
    testit::assert("Coordinates in GRanges are variable length!", length(unique(width(regions.gr)))==1)
    if(keep_region_coordinates){
      reads <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1])
    } else {
      regions.gr <- resize(regions.gr, width=downstream)
      regions.gr <- resize(regions.gr, width=upstream + width(regions.gr), fix="end")

      reads <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1])
    }
  }else{
    reads <- matrix(nrow=length(regions.gr), ncol=approx_nu)
  }
  i_p <- which(as.logical(strand(regions.gr) == "+" | strand(regions.gr) == "*"))
  i_n <- which(as.logical(strand(regions.gr) == "-"))

  message("There are ", length(i_p), " positive granges and ", length(i_n), " negative granges")

  # if(class(sample.cov) == "character") {
  #   sample.cov <- import(sample.cov, which=regions.gr)
  # }
  if(diff_length == FALSE){
	  if(length(i_p) > 0) reads[i_p, ] <- read_matrix(regions.gr[i_p], sample.cov)
		if(length(i_n) > 0) reads[i_n, ] <- read_matrix(regions.gr[i_n], sample.cov, reverse_reads=TRUE)
	}else{
		if(length(i_p) > 0)reads[i_p, ] <- read_matrix(regions.gr[i_p], sample.cov, df=diff_length, nu=approx_nu)
		if(length(i_n) > 0)reads[i_n, ]<- read_matrix(regions.gr[i_n], sample.cov, reverse_reads=TRUE, df=diff_length,  nu=an)
	}
  reads
}

exo_metapeak_matrix <- function(regions.gr, sample, upstream=100, downstream=100, keep_region_coordinates = FALSE) {
  if(keep_region_coordinates){
    reads <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1])
  } else {
    regions.gr <- GenomicRanges::resize(regions.gr, width=downstream)
    regions.gr <- GenomicRanges::trim(resize(regions.gr, width=upstream + width(regions.gr), fix="end")) #resize regions upstream on top of the previous GRanges file
    regions.gr <- regions.gr[width(regions.gr) == upstream+downstream] #regions must be equal to upstream+downstream
  }

  i_p <- which(as.vector(strand(regions.gr) == "+") | as.vector(strand(regions.gr) == "*")) #create i_p for where regions are +
  i_n <- which(as.vector(strand(regions.gr) == "-")) #create i_n for where regions are -

  message("There are ", length(i_p), " positive granges and ", length(i_n), " negative granges")

  reads.p <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1]) #make matrix size (num of regions)x(width of regions, def:200)
  reads.n <- reads.p #same empty matrix for neg as pos strand

	  if(length(i_p) > 0) {
	    reads.p[i_p, ] <- read_matrix(regions.gr[i_p], sample$pos, df=FALSE, nu=NULL)
	    reads.n[i_p, ] <- abs(read_matrix(regions.gr[i_p], sample$neg, df=FALSE, nu=NULL))
	  } #writes in all the "positive" GRanges as they are listed in that order

	  if(length(i_n) > 0) {
	    reads.p[i_n, ] <- abs(read_matrix(regions.gr[i_n], sample$neg, reverse_reads=TRUE, df=FALSE, nu=NULL))
	    reads.n[i_n, ] <- read_matrix(regions.gr[i_n], sample$pos, reverse_reads=TRUE, df=FALSE, nu=NULL)
	  } #accounts for orientation and writes in the "negative" GRanges in their spots, but reverses them

  list(pos=reads.p, neg=reads.n)
}

standard_metapeak <- function(gr, sample, upstream=100, downstream=100, sample_name=NA, smooth=NA, different_length=FALSE, approx_n=25) {
  message("standard metapeak: ", sample_name)
  if(different_length==FALSE){
	  reads <- standard_metapeak_matrix(gr, sample, upstream, downstream, diff_length=FALSE, approx_nu=NULL)
	  reads.df <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
	                        reads=colMeans(reads),
	                        sample_name=sample_name)
	  if(!is.na(smooth)) reads.df$reads <- as.numeric(runmean(Rle(reads.df$reads), k=smooth, endrule="constant"))

	}else{
	 reads <- standard_metapeak_matrix(gr, sample, dl= T, an=approx_n)
  	 reads.df <- data.frame(tss_distance=1:ncol(reads),
                        reads=colMeans(reads),
                        sample_name=sample_name, stringsAsFactors = FALSE)
	}

  reads.df
}

exo_metapeak <- function(gr, sample, upstream=100, downstream=100, sample_name=NA, smooth=NA) {
  message("exo metapeak: ", sample_name)
  reads.list <- exo_metapeak_matrix(gr, sample, upstream, downstream) #get list of reads from function exo_metapeak_matrix

  reads.p <- reads.list$pos
  reads.n <- reads.list$neg

  df.p <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                     reads=colMeans(reads.p),
                     strand="+")

  df.n <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                     reads=colMeans(reads.n) *(-1),
                     strand="-")

  if(!is.na(smooth)) {
    df.n$reads <- as.numeric(runmean(Rle(df.n$reads), k=smooth, endrule="constant"))
    df.p$reads <- as.numeric(runmean(Rle(df.p$reads), k=smooth, endrule="constant"))
  }

  reads.df <- rbind(df.p, df.n)
  reads.df$sample_name <- sample_name
  reads.df$sample <- paste(reads.df$sample_name, reads.df$strand)
  reads.df
}

standard_metapeak_across_large_region_set<-function(gr, sample, upstream = 100, downstream = 100, cores = 8, chunk_size = 50000){
  chunks<-seq(1, length(gr), chunk_size)
  chunks.df<-mclapply(chunks, function(x){
    upper_limit<-min(x+(chunk_size-1), length(gr)); idx<-x:upper_limit
    df<-standard_metapeak(gr = gr[idx], sample = sample, upstream = upstream, downstream = downstream)
    n_chunk<-length(idx)
    df$reads_sum<-df$reads * n_chunk #get sum of reads across regions rather than mean
    return(df)
  }, mc.cores = cores) %>%
    rbindlist() %>%
    dplyr::group_by(tss_distance) %>%
    #once all regions are found, find mean
    dplyr::summarize(reads = sum(reads_sum, na.rm = T)/length(gr))
  colnames(chunks.df)<-c("position", "reads")
  return(chunks.df)
}

#The functions below rely on the load_bigwig function that deals with multiple samples. If you would like to use this strategy, refer to sample_common.R code.

get_exo_metapeak <- function(gr, sample, upstream=100, downstream=101, smooth=NA, sample_format = "merged", sample_name = NA){
  if(sample_format == "merged"){
    if(is.na(sample_name)){
      sample_name <-sample
      }
    sample_path = load_bigwig(sample)
    metapeak <- exo_metapeak(gr, sample_path,upstream=upstream, downstream=downstream, sample_name=sample_name, smooth=smooth)
  }
  if(sample_format == "separate"){
    if(is.na(sample_name)){
      sample_name <-sample
    }
    sample_path1 <- load_bigwig(sample, sample_format = "separate")[[1]]
    sample_path2 <- load_bigwig(sample, sample_format = "separate")[[2]]

    gr.ex <- resize(gr, upstream, "end") %>% resize(., downstream + upstream, "start")

    cov1 <- list(pos=rtracklayer::import.bw(sample_path1$pos, which=gr.ex, as="RleList"), neg = rtracklayer::import.bw(sample_path1$neg, which=gr.ex, as="RleList"))
    cov2 <- list(pos=rtracklayer::import.bw(sample_path2$pos, which=gr.ex, as="RleList"), neg = rtracklayer::import.bw(sample_path2$neg, which=gr.ex, as="RleList"))

    cov <- list(pos = cov1$pos + cov2$pos, neg = cov1$neg + cov2$neg)
    metapeak <- exo_metapeak(gr, cov, upstream=upstream, downstream=downstream, sample_name =sample_name, smooth=smooth)
  }
  if(sample_format == "data"){
    metapeak <- exo_metapeak(gr, sample, upstream=upstream, downstream=downstream, sample_name =sample_name, smooth=smooth)
  }

  metapeak
}

get_standard_metapeak <- function(gr, sample, upstream=100, downstream=101, smooth=NA){
  sample_path = rtracklayer::import(sample) #returns problems with getSeqlevels (try just immediately going to standard_metapeak)
  metapeak <- standard_metapeak(gr, sample_path, upstream=upstream, downstream=downstream, sample_name=sample, smooth=smooth)
  metapeak
}

get_exo_matrix <- function(gr, sample, upstream=100, downstream=101, sample_format = "merged"){
  if(sample_format == "merged"){
    sample_path = load_bigwig(sample)
    exo_matrix <- exo_metapeak_matrix(gr, sample_path,upstream=upstream, downstream=downstream)
  }
  if(sample_format == "separate"){
    sample_path1 <- load_bigwig(sample, sample_format = "separate")[[1]]
    sample_path2 <- load_bigwig(sample, sample_format = "separate")[[2]]

    gr.ex <- resize(gr, upstream, "end") %>% resize(., downstream + upstream, "start")

    cov1 <- list(pos=rtracklayer::import.bw(sample_path1$pos, which=gr.ex, as="RleList"), neg = rtracklayer::import.bw(sample_path1$neg, which=gr.ex, as="RleList"))
    cov2 <- list(pos=rtracklayer::import.bw(sample_path2$pos, which=gr.ex, as="RleList"), neg = rtracklayer::import.bw(sample_path2$neg, which=gr.ex, as="RleList"))

    cov <- list(pos = cov1$pos + cov2$pos, neg = cov1$neg + cov2$neg)

    exo_matrix <-exo_metapeak_matrix(gr, cov,upstream=upstream, downstream=downstream)
  }
  if(sample_format == "data"){
    exo_matrix <-exo_metapeak_matrix(gr, sample,upstream=upstream, downstream=downstream)
  }
  exo_matrix
}
