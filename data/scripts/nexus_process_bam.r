suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

option_list <- list(
  make_option(c("-f", "--file"),
              type="character",
              help="Path of BAM file to process"),
  make_option(c("-n", "--name"),
              type="character",
              help="Name for resulting R object"),
  make_option(c("-p", "--paired"),
              action="store_true",
              default=FALSE,
              help="Paired-end mode"),
  make_option(c("-u", "--unique"),
 		type="logical",
		default=TRUE,
		help="remove duplicates based on random barcode"),
  make_option(c("-o", "--output_filtered"),
              type="logical",
              default=F,
              help="logical indicating whether to output a _filtered.bam in addition to unfiltered"))

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(rtracklayer))

id <- paste0("[", opt$name, "] ")

stopifnot(file.exists(opt$file))

if(opt$paired) {
  message(id, "Reading paired-end BAM: ", opt$file)
  bam <- readGAlignmentPairs(opt$file, param=ScanBamParam(what="qname"))
  bam.gr <- granges(bam)
  bam.gr$barcode <- gsub("^(.*)_pair_.*$", "\\1", mcols(first(bam))$qname)
} else {
  message(id, "Reading BAM: ", opt$file)
  bam.gr <- granges(readGAlignments(opt$file, use.names=TRUE))
  bam.gr$barcode <- names(bam.gr)
  names(bam.gr) <- NULL
}

message(id, pn(length(bam.gr)), " fragments")

if(opt$unique){
	message(id, "Removing barcode duplicates...")
	bam_uniq.gr <- GenomicRanges::split(bam.gr, bam.gr$barcode) %>%
	               unique %>%
	               unlist(use.names=FALSE)
}else{
	bam_uniq.gr <- bam.gr
}


message(id, "Saving ", pn(length(bam_uniq.gr)), " fragments...")
saveRDS(bam_uniq.gr, file=paste(opt$name, ".granges.rds", sep=""))

if(opt$output_filtered==T){
  export(bam_uniq.gr, paste(opt$name, ".bam",sep=""), format = "bam")
}
