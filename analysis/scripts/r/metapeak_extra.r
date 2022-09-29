library(GenomicRanges)
library(reshape2)
library(lattice)
library(magrittr)
library(IRanges)

#20180525_modification: See commented out code with Views(). Upgrading to R version 5.0 will cause RangesList to become deprecated and GRangesList is not compatible with this function. In order to circumvent, a newer updated version of splitting a GRanges into GRangesList was added.

read_matrix <- function(gr, cov, reverse_reads=FALSE, df=F, nu=25) {
  if(class(cov) == "character") {
      cov <- import.bw(cov, which=gr, as="RleList")
    }
  transform_function <- if(reverse_reads) { rev } else { identity }
  o <- order(gr)
  gr <- gr[o]

  if(df==T){
	reads.list <- regionApply_list(gr,cov, as.numeric )
  	reads.list <-lapply(reads.list, function(reads) { approx(reads, n=nu)$y })
	reads.m <- matrix(as.numeric(unlist(reads.list, use.name=F)), nrow=length(reads.list), byrow=T)
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

standard_metapeak_matrix <- function(regions.gr, sample.cov, upstream=100, downstream=100, diff_length = F, approx_nu=25) {
  if(diff_length == F){
	    regions.gr <- resize(regions.gr, width=downstream)
	    regions.gr <- resize(regions.gr, width=upstream + width(regions.gr), fix="end")

	    reads <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1])
	}else{
		reads <- matrix(nrow=length(regions.gr), ncol=approx_nu)
	}
  i_p <- which(strand(regions.gr) == "+" | strand(regions.gr) == "*")
  i_n <- which(strand(regions.gr) == "-")

  message("There are ", length(i_p), " positive granges and ", length(i_n), " negative granges")

  # if(class(sample.cov) == "character") {
  #   sample.cov <- import(sample.cov, which=regions.gr)
  # }
  if(diff_length == F){
	  	if(length(i_p) > 0) reads[i_p, ] <- read_matrix(regions.gr[i_p], sample.cov)
		if(length(i_n) > 0) reads[i_n, ] <- read_matrix(regions.gr[i_n], sample.cov, reverse_reads=TRUE)
	}else{
		if(length(i_p) > 0)reads[i_p, ] <- read_matrix(regions.gr[i_p], sample.cov, df=diff_length, nu=approx_nu)
		if(length(i_n) > 0)reads[i_n, ]<- read_matrix(regions.gr[i_n], sample.cov, reverse_reads=TRUE, df=diff_length,  nu=an)
	}
  reads
}

exo_metapeak_matrix <- function(regions.gr, sample, upstream=100, downstream=100) {
  regions.gr <- resize(regions.gr, width=downstream) #resize regions downstream to discrete value
  regions.gr <- trim(resize(regions.gr, width=upstream + width(regions.gr), fix="end")) #resize regions upstream on top of the previous GRanges file
  regions.gr <- regions.gr[width(regions.gr) == upstream+downstream] #regions must be equal to upstream+downstream

  i_p <- which(strand(regions.gr) == "+" | strand(regions.gr) == "*") #create i_p for where regions are +
  i_n <- which(strand(regions.gr) == "-") #create i_n for where regions are -

  message("There are ", length(i_p), " positive granges and ", length(i_n), " negative granges")

  reads.p <- matrix(nrow=length(regions.gr), ncol=width(regions.gr)[1]) #make matrix size (num of regions)x(width of regions, def:200)
  reads.n <- reads.p #same empty matrix for neg as pos strand

	  if(length(i_p) > 0) {
	    reads.p[i_p, ] <- read_matrix(regions.gr[i_p], sample$pos, df=F, nu=NULL)
	    reads.n[i_p, ] <- abs(read_matrix(regions.gr[i_p], sample$neg, df=F, nu=NULL))
	  } #writes in all the "positive" GRanges as they are listed in that order

	  if(length(i_n) > 0) {
	    reads.p[i_n, ] <- abs(read_matrix(regions.gr[i_n], sample$neg, reverse_reads=TRUE, df=F, nu=NULL))
	    reads.n[i_n, ] <- read_matrix(regions.gr[i_n], sample$pos, reverse_reads=TRUE, df=F, nu=NULL)
	  } #accounts for orientation and writes in the "negative" GRanges in their spots, but reverses them



  list(pos=reads.p, neg=reads.n)
}

standard_metapeak <- function(gr, sample, upstream=100, downstream=100, sample_name=NA, smooth=NA, different_length=F,approx_n=25) {
  message("standard metapeak: ", sample_name)
  if(different_length==F){
    reads <- standard_metapeak_matrix(gr, sample, upstream, downstream, diff_length=F, approx_nu=NULL)
    reads.df <- data.frame(tss_distance=(-1 * upstream):(downstream - 1),
                           reads=colMeans(reads),
                           sample_name=sample_name)
    if(!is.na(smooth)) reads.df$reads <- as.numeric(runmean(Rle(reads.df$reads), k=smooth, endrule="constant"))
    
  }else{
    reads <- standard_metapeak_matrix(gr, sample, dl= T, an=approx_n)
    reads.df <- data.frame(tss_distance=1:ncol(reads),
                           reads=colMeans(reads),
                           sample_name=sample_name)
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
    reads.df
}


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

    cov1 <- list(pos=import.bw(sample_path1$pos, which=gr.ex, as="RleList"), neg = import.bw(sample_path1$neg, which=gr.ex, as="RleList"))
    cov2 <- list(pos=import.bw(sample_path2$pos, which=gr.ex, as="RleList"), neg = import.bw(sample_path2$neg, which=gr.ex, as="RleList"))

    cov <- list(pos = cov1$pos + cov2$pos, neg = cov1$neg + cov2$neg)
    metapeak <- exo_metapeak(gr, cov, upstream=upstream, downstream=downstream, sample_name =sample_name, smooth=smooth)
  }
  if(sample_format == "data"){
    metapeak <- exo_metapeak(gr, sample, upstream=upstream, downstream=downstream, sample_name =sample_name, smooth=smooth)
  }

  metapeak
}

get_standard_metapeak <- function(gr, sample, upstream=100, downstream=101, smooth=NA){
  sample_path = import(sample)
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

    cov1 <- list(pos=import.bw(sample_path1$pos, which=gr.ex, as="RleList"), neg = import.bw(sample_path1$neg, which=gr.ex, as="RleList"))
    cov2 <- list(pos=import.bw(sample_path2$pos, which=gr.ex, as="RleList"), neg = import.bw(sample_path2$neg, which=gr.ex, as="RleList"))

    cov <- list(pos = cov1$pos + cov2$pos, neg = cov1$neg + cov2$neg)

    exo_matrix <-exo_metapeak_matrix(gr, cov,upstream=upstream, downstream=downstream)
  }
  if(sample_format == "data"){
    exo_matrix <-exo_metapeak_matrix(gr, sample,upstream=upstream, downstream=downstream)
  }
  exo_matrix
}




########################################################################################################################
#multiple_grange_metapeak_function
########################################################################################################################
get_multi_gr_standard_metapeak <- function(gr.list = gr.list, bw = bw, upstream = 500, downstream = 500, smooth = NA){
    mclapply(names(gr.list), function(x){
    gr <- gr.list[[x]]
    standard_metapeak(gr,bw, upstream = upstream, downstream = downstream, smooth = smooth, sample_name = x)
    },mc.cores = 10) %>% dplyr::bind_rows()
}


get_multi_gr_exo_metapeak <- function(gr.list = gr.list, bw = bw, upstream = 500, downstream = 500, smooth = NA){
    mclapply(names(gr.list), function(x){
    gr <- gr.list[[x]]
    exo_metapeak(gr,bw, upstream = upstream, downstream = downstream, smooth = smooth, sample_name = x)
    },mc.cores = 10) %>% dplyr::bind_rows()
}


########################################################################################################################
#multiple_bigwig_metapeak_function
########################################################################################################################
get_multi_bw_standard_metapeak <- function(gr = gr, bw.list = bw.list, upstream = 500, downstream = 500,smooth = NA){
    mclapply(names(bw.list), function(x){
    bw <- bw.list[[x]]
    standard_metapeak(gr,bw, upstream = upstream, downstream = downstream, smooth = smooth, sample_name = x)
    },mc.cores = 10) %>% dplyr::bind_rows()
}


get_multi_bw_exo_metapeak <- function(gr = gr, bw.list = bw.list, upstream = 500, downstream = 500,smooth = NA){
    mclapply(names(bw.list), function(x){
    bw <- bw.list[[x]]
    exo_metapeak(gr,bw, upstream = upstream, downstream = downstream, smooth = smooth, sample_name = x)
    },mc.cores = 10) %>% dplyr::bind_rows()
}

########################################################################################################################
#multiple_grange_and_bigwig_metapeak_function
########################################################################################################################
get_multi_gr_bw_standard_metapeak <- function(gr.list = gr.list, bw.list = bw.list, upstream = 500, downstream = 500, smooth = NA){
  mclapply(names(gr.list), function(x){
    gr <- gr.list[[x]]
    get_multi_bw_standard_metapeak(gr,bw.list, upstream = upstream, downstream = downstream, smooth = smooth) %>% mutate(grange_name = x) # The ample_name will be the bw name
  },mc.cores = 10) %>% dplyr::bind_rows() %>% rename(bw_name = sample_name)
}


get_multi_gr_bw_exo_metapeak <- function(gr.list = gr.list, bw = bw.list, upstream = 500, downstream = 500, smooth = NA){
  mclapply(names(gr.list), function(x){
    gr <- gr.list[[x]]
    get_multi_bw_exo_metapeak(gr,bw.list, upstream = upstream, downstream = downstream, smooth = smooth) %>%mutate(grange_name = x)# The ample_name will be the bw name
  },mc.cores = 10) %>% bind_rows() %>% rename(bw_name = sample_name)
}



########################################################################################################################
#metapeak_plotting_function
########################################################################################################################
plot_standard_metapeak <- function(metapeak, name = name, xlab = xlab, ylab = ylab, ncol = NA){
    x <- ggplot(metapeak, aes(x=tss_distance, y=reads)) +
         geom_area(position="identity", alpha = 0.2) +
         geom_smooth(method = "loess", span = 0.3, color = "dodgerblue4")+
         ggtitle(name) +
            facet_wrap(~sample_name, scales = "fixed", ncol = ncol) +
            theme(legend.position = "none") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black")) +
            xlab(xlab) +
            ylab(ylab)
    x
        }

plot_exo_metapeak<- function(metapeak, pos.col, neg.col, name = name, xlab = xlab, ylab = ylab, ncol = NA){
    x <- ggplot(metapeak, aes(x=tss_distance, y=reads, fill=strand)) +
          geom_area(position="identity") +
          scale_fill_manual(values=c(pos.col, neg.col)) +
          ggtitle(name) +
            facet_wrap(~sample_name, scales = "fixed", ncol = ncol) +
            theme(legend.position = "none") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black")) +
            xlab(xlab) + ylab(ylab)

            x
                }
