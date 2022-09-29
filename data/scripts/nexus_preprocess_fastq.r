suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-f", "--file"),
              type="character",
              help="Path of FASTQ file to process"),
  make_option(c("-t", "--trim"),
              type="integer",
              default=0,
              help="Pre-trim all reads to this length before processing"),
  make_option(c("-k", "--keep"),
              type="integer",
              default=18,
              help="Minimum number of bases required after barcode to keep read"),
  make_option(c("-b", "--barcodes"),
              type="character",
              default="CTGA",
              help="Barcode sequences (comma-separated) that follow random barcode"),
  make_option(c("-r", "--randombarcode"),
             type="integer",
             default=5,
             help="Number of bases at the start of each read used for random barcode"),
  make_option(c("-c", "--chunksize"),
              type="integer",
              default=1000,
              help="Number of reads to process at once (in thousands)"),
  make_option(c("-o", "--output"),
              type="character",
              help="Output FASTQ file (gzip compressed)"),
  make_option(c("-p", "--processors"),
              type="character",
              default=2,
              help="Number of simultaneous processing cores to utilize"))


opt <- parse_args(OptionParser(option_list=option_list))

if(file.exists(opt$output)) stop("Output file ", opt$output, " already exists.")

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(BiocParallel))

pn <- function(value) {
  prettyNum(value, big.mark=",")
}

process_chunk <- function(fq_chunk, opt) {

  if(opt$trim > 0) { #Trim data if an option
    fq_chunk <- narrow(fq_chunk, start=1, end=pmin(width(fq_chunk), opt$trim))
  }

  output_file <- opt$output

  barcodes <- strsplit(opt$barcodes, split=",")[[1]] #split barcode strings
  if(length(unique(str_length(barcodes))) != 1) {
    stop("Fixed barcodes must all be the same length: ", barcodes)
  }
  barcode_start <- opt$randombarcode + 1 #define start
  barcode_end   <- barcode_start + str_length(barcodes[1]) - 1 #define end

  fq.lockfile <- paste0(output_file, ".lock")


  barcode_reads <- narrow(sread(fq_chunk), start=barcode_start, end=barcode_end) #define regions of barcode_reads

  bc_matches <- as.list(rep(NA, times=length(barcodes))) #allocate
  names(bc_matches) <- barcodes #name

  for(barcode in barcodes) { #calculate barcode matches for each barcode component and store in bc_matches
    matches <- elementNROWS(vmatchPattern(barcode, barcode_reads, fixed=FALSE)) == 1
    n_count <- elementNROWS(vmatchPattern("N", barcode_reads, fixed=TRUE))
    bc_matches[[barcode]] <- which(matches == TRUE & n_count <= 1) #where a match occurs <= 1 N's are present
  }

  # don't allow a read to match multiple barcodes
  previous_matches <- c()
  for(i in seq_along(bc_matches)) {
    bc_matches[[i]] <- bc_matches[[i]][!bc_matches[[i]] %in% previous_matches]
    previous_matches <- c(previous_matches, bc_matches[[i]])
  }
  total_matches <- sum(elementNROWS(bc_matches))

  message("[", opt$file, "] ", pn(length(fq_chunk)), " reads with ", pn(total_matches), " barcode matches")

  if(total_matches > 0) {
    for(barcode in barcodes) {
      fq.matched <- fq_chunk[bc_matches[[barcode]]]
      if(length(fq.matched) == 0) next

      # Reject reads that are too short
      fq.matched <- fq.matched[width(sread(fq.matched)) >= barcode_end + opt$keep]

      # Keep random barcode
      random_bc  <- substr(sread(fq.matched), 1, opt$randombarcode)

      fq.matched <- narrow(fq.matched, start=barcode_end + 1, width=width(fq.matched) - (barcode_end + 1))

      fq.new <- ShortReadQ(sread   = sread(fq.matched),
                           quality = quality(fq.matched),
                           id      = BStringSet(paste0(random_bc, "_", barcode)))

      lock_status <- system(paste("lockfile", "-1", fq.lockfile), intern=FALSE)
      if(lock_status > 0) stop("lockfile command failed.")
      output_mode <- ifelse(file.exists(output_file), "a", "w")
      writeFastq(fq.new, file=output_file, compress=TRUE, mode=output_mode)
      file.remove(fq.lockfile)
    }
  }
  TRUE
}

yieldHelper <- function() {
  fq <- yield(fqstream, withIds=FALSE)
  if(length(fq) > 0) {
    fq
  } else {
    NULL
  }
}

# Restrict FastqStreamer threading
nothing <- .Call(ShortRead:::.set_omp_threads, 1L)

fqstream <- FastqStreamer(opt$file, n=opt$chunksize * 1000)
bpparam <- MulticoreParam(workers=opt$processors, stop.on.error=TRUE, log=FALSE)

results <- bpiterate(ITER=yieldHelper, FUN=process_chunk, opt=opt, BPPARAM=bpparam)
close(fqstream)
