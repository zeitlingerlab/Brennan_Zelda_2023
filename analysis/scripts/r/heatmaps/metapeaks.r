library(GenomicRanges)
library(reshape2)
library(lattice)

read_matrix <- function(gr, cov, reverse_reads=FALSE) {
  transform_function <- if(reverse_reads) { rev } else { identity }
  o <- order(gr)
  gr <- gr[o]
  rl <- as(gr, "RangesList")
  view <- RleViewsList(rleList=cov[names(rl)], rangesList=rl)
  reads.list <- viewApply(view, function(x) { transform_function(as.numeric(x)) })
  reads.m <- matrix(unlist(sapply(reads.list, as.numeric)), nrow=length(gr), byrow=TRUE)
  reads.m[o, ] <- reads.m
  reads.m
}

standard_metapeak_matrix <- function(regions.gr, sample.cov, upstream=100, downstream=100) {
  regions.gr <- resize(regions.gr, width=downstream)
  regions.gr <- resize(regions.gr, width=upstream + width(regions.gr), fix="end")
  
  failed_resize <- which(width(regions.gr) != upstream + downstream)
  if(length(failed_resize) > 0) {
    stop("The following regions.gr elements cannot be resized due to chromosome boundaries: ", paste0(failed_resize, collapse=", "))
  }
  
  regions.p <- regions.gr[strand(regions.gr) == "+" | strand(regions.gr) == "*"]
  regions.n <- regions.gr[strand(regions.gr) == "-"]
  
  reads <- NULL
  
  if(class(sample.cov) == "character") {
    sample.cov <- import.bw(sample.cov, which=regions.gr, as="RleList")
  }
  
  if(length(regions.p) > 0) reads <- rbind(reads, read_matrix(regions.p, sample.cov))
  if(length(regions.n) > 0) reads <- rbind(reads, read_matrix(regions.n, sample.cov, reverse_reads=TRUE))

  stopifnot(nrow(reads) == length(regions.gr))
  reads
}

standard_metapeak <- function(gr, sample.cov, upstream=100, downstream=100, sample_name=NA, smooth=NA) {
  message("standard metapeak: ", sample_name)
  
  reads <- standard_metapeak_matrix(gr, sample.cov, upstream, downstream)
  
  reads.df <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                        reads=colMeans(reads), 
                        sample_name=sample_name)
  if(!is.na(smooth)) reads.df$reads <- as.numeric(runmean(Rle(reads.df$reads), k=smooth, endrule="constant"))
  reads.df  
}

enrichment_metapeak <- function(gr, sample.cov, bg.cov, upstream=100, downstream=100, sample_name=NA, smooth=NA) {
  message("enrichment metapeak: ", sample_name)
  
  reads.ip <- standard_metapeak_matrix(gr, sample.cov, upstream, downstream)
  reads.bg <- standard_metapeak_matrix(gr, bg.cov, upstream, downstream)
  
  ts.ip <- total_signal(sample.cov)
  ts.bg <- total_signal(bg.cov)
  
  reads.df <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                        reads=colMeans(reads.ip) / ts.ip, 
                        background=colMeans(reads.bg) / ts.bg,
                        sample_name=sample_name)

  reads.df$enrichment <- reads.df$reads / reads.df$background

  if(!is.na(smooth)) {
    reads.df$reads      <- as.numeric(runmean(Rle(reads.df$reads), k=smooth, endrule="constant"))
    reads.df$background <- as.numeric(runmean(Rle(reads.df$background), k=smooth, endrule="constant"))
    reads.df$enrichment <- reads.df$reads / reads.df$background
  }
  
  reads.df  
}

base_frequencies <- function(gr, upstream, downstream, genome=Dmelanogaster) {

  gr <- resize(gr, width=downstream)
  gr <- resize(gr, width=upstream + width(gr), fix="end")

  m <- consensusMatrix(getSeq(genome, gr))[1:4, ]
  m <- m / colSums(m)
  df.bases <- as.data.frame(t(m))
  df.bases$tss_distance <- (-1 * upstream):(downstream - 1)
  df.bases <- melt(df.bases, id.var="tss_distance")
  names(df.bases)[2:3] <- c("base", "frequency")
  df.bases
}


