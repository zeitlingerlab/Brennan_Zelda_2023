suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-f", "--file"),
              type="character",
              help="Path to granges file"),
  make_option(c("-n", "--name"),
              type="character",
              help="name of the output data"),
  make_option(c("-m", "--normalization"),
              default = "T",
              type="character",
              help="parameter to indicate if normalize to the total reads, default is true"),
  make_option(c("-t", "--type"),
              default = "atacseq",
              type="character",
              help="experiment type, atacseq or mnase"))

opt <- parse_args(OptionParser(option_list=option_list))

total_reads_normalization <- function(gr){
  total_reads <- length(gr)
  gr.cov <- coverage(gr) / total_reads * 1000000
  gr.cov
}

if(opt$type == "atacseq"){
  message("sample being processed is atacseq")
  message("reading the granges file")
  gr <- readRDS(opt$file)
  if(opt$normalization == "T"){
    message("generating the normalized coverage files")
    gr.cov <- total_reads_normalization(gr)
  }else{
    message("generating the non-normalized coverage files")
    gr.cov <- coverage(gr)
  }
  
  message("exporting atac bigwig files")
  export(gr.cov, paste0(opt$name, ".bw"))
  
}

if(opt$type == "mnase"){
  message("sample being processed is mnase")
  message("reading the granges file")
  gr <- get(load(opt$file))
  if(opt$normalization == "T"){
    message("generating the normalized coverage files")
    gr.cov <- total_reads_normalization(gr)
  }else{
    message("generating the non-normalized coverage files")
    gr.cov <- coverage(gr)
  }
  
  message("exporting mnase bigwig files")
  export(gr.cov, paste0(opt$name, ".bw"))
}
