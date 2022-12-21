#!/bin/env Rscript

# Call from command line

library(geneBodyCoverage)
library(getopt)




# Read command line arguments ----
if(! interactive()){
  spec <- matrix(c(
    'bam_path',  'b', 1, "character", "Path to BAM file",
    'bed_path',  'r', 1, "character", "Path to reference BED12 file",
    'endedness', 'e', 1, "integer",   "From which end to start counting (3 or 5)",
    'delimiter', 'd', 1, "character", "Delimiter of the bed file (e.g. space or tab)",
    'chrom_list','c', 2, "character", "Name of the list of chromosomes, or blank (elegans, hg38 accepted)",
    'out_path',  'o', 1, "character", "Path to save output (must not exists, if not -f)",
    'nproc',     'p', 2, "integer",   "Number of threads",
    'future_mem','m', 2, "integer",   "Max memory for future globals",
    'force',     'f', 0, "logical",   "Set to overwrite existing output",
    'min_length','l', 2, "integer",   "Minimum transcript length, any transcript shorter will be ignored (default: 100 bp)",
    'help',      'h', 0, "logical",   "Print this help"
  ), byrow=TRUE, ncol=5)

  opt <- getopt(spec)
} else{
  # Options for experimenting
  opt <- list(
    bam_path = "tests/testthat/fixtures/simple.bam",
    bed_path = "tests/testthat/fixtures/c_elegans.PRJNA13758.WS281_filtered.bed",
    endedness = 3L,
    delimiter = "\t",
    chrom_list = "elegans",
    out_path = tempfile("out_gbc.tsv")
  )
}


if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}


cat("Starting: ",date(),"\n\n")


# check options ----

if(file.exists(opt$out_path) && is.null(opt$force)){
  stop("Output file already exists. Use --force to overwrite.")
}

if(is.null(opt$nproc)) opt$nproc <- 1
if(!is.numeric(opt$nproc)) stop("Error: nproc must be numeric, got ", opt$nproc)
if(opt$nproc > 16) stop("nproc too high")

if(! opt$endedness %in% c(3L, 5L)){
  stop("Endedness must be 3 or 5, not ", opt$endedness)
}

if(nchar(opt$delimiter) != 1L){
  stop("Invalid delimiter: ", opt$delimiter," has ",nchar(opt$delimiter)," characters.")
}

if(! is.null(opt$chrom_list)){
  if(opt$chrom_list == "elegans"){
    opt$chrom_list <- as.roman(1:5) |>
      as.character() |>
      c("X","MtDNA")
  } else if(opt$chrom_list == "hg38"){
    opt$chrom_list <- c("hg38_5","hg38_13","hg38_18")
  } else{
    stop("Unrecognized chromosome list: ", opt$chrom_list)
  }
}

# allow 1.5 GB = 1610612736
if(is.null(opt$future_mem)){
  options(future.globals.maxSize = 1610612736)
} else{
  options(future.globals.maxSize = opt$future_mem)
}

if(is.null(opt$min_length)){
  opt$min_length <- 100
}





# transcript coordinates bin breaks (we will count reads in each bin, starting from 3' end)
# Note for Parse: the longest transcript in chromosomes 5,13,18 appears to be 21kb, and there is a 5' bias.

bins <- c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,
          1700,1800,1900,2000,2100,2200,2300,2400,2500,2700,2900,3100,3400,3700,
          4000,4500,5000,5500,6000,6500,7000,8000,9000,10000,Inf)



cat("Using ",length(bins), " bins:\n  ", bins, "\n\n")




# Process reference ----

cat("Read bed12 file with tx and exon coordinates from file: ", opt$bed_path, "\n")
gene_struct_raw <- readr::read_delim(opt$bed_path, delim = opt$delimiter,
                          col_names = c("chr","transcript_start","transcript_end","transcript_name","score",
                                        "strand","thick_start","thick_end","rgb",
                                        "nb_exons","exons_lengths","exons_starts"),
                          col_types = readr::cols(
                            chr =              readr::col_factor(opt$chrom_list), #1
                            transcript_start = readr::col_integer(),              #2 on the + strand (even if tx on -)
                            transcript_end =   readr::col_integer(),              #3
                            transcript_name =  readr::col_character(),            #4
                            score =            readr::col_factor("0"),            #5 ignored
                            strand =           readr::col_factor(c("+","-")),     #6
                            thick_start =      readr::col_integer(),              #7 ignored
                            thick_end =        readr::col_integer(),              #8 ignored
                            rgb =              readr::col_factor("0"),            #9 ignored
                            nb_exons =         readr::col_integer(),              #10
                            exons_lengths =    readr::col_character(),            #11
                            exons_starts =     readr::col_character()             #12 on the + strand (even if tx on -)
                          ))
readr::stop_for_problems(gene_struct_raw)

gene_struct <- gene_struct_raw |>
  dplyr::mutate(exons_lengths = strsplit(exons_lengths, ",") |>
           purrr::map(as.integer),
         exons_starts = strsplit(exons_starts, ",") |>
           purrr::map(as.integer))


# Read alignments
cat("\nRead alignments from file: ", opt$bam_path, "\n\n")
bam_content <- GenomicAlignments::readGAlignments(opt$bam_path)







# Count reads in bins ----






if(opt$nproc == 1 || ! requireNamespace("furrr", quietly = TRUE)){

  if(opt$nproc > 1){
    cat("Warning: Package {furrr} not available, can not use multiple processors!\n")
  }

  cat("\nFor each transcript, getting bins coordinates (using a single thread)\n")

  tx_struct <- gene_struct |>
    dplyr::filter(transcript_end -  transcript_start > opt$min_length) |>
    dplyr::mutate(bins_gr = purrr::pmap(list(chr,
                                             strand,
                                             transcript_start,
                                             transcript_end,
                                             exons_lengths,
                                             exons_starts),
                                        bins_to_granges,
                                        bins,
                                        opt))

  cat("Counting alignments (using a single thread)\n")

  bin_quantifs <- tx_struct |>
    dplyr::select(bins_gr,
                  chr,
                  strand) |>
    purrr::pmap(bins_coverage,
                bam_content)

} else{

  future::plan(multicore, workers = opt$nproc)


  cat("\nFor each transcript, getting bins coordinates (using ",opt$nproc," threads)\n")
  tx_struct <- gene_struct |>
    dplyr::filter(transcript_end -  transcript_start > opt$min_length) |>
    dplyr::mutate(bins_gr = furrr::future_pmap(list(bins_in_genomic_coords,
                                             chr,
                                             strand,
                                             transcript_start,
                                             transcript_end,
                                             exons_lengths,
                                             exons_starts),
                                        bins_to_granges,
                                        bins,
                                        opt))

  cat("Count alignments with ",opt$nproc," threads...\n")
  bin_quantifs <- tx_struct |>
    dplyr::select(bins_gr,
                  chr,
                  strand) |>
    furrr::future_pmap(bins_coverage,
                       bam_content)
}





# Save result ----
cat("Done.\n")
tx_quantif <- tx_struct |>
  dplyr::mutate(quantification = bin_quantifs)

cat("Saving results at ", opt$out_path, "\n")
tx_quantif |>
  dplyr::select(transcript_name,
                chr, transcript_start, transcript_end, strand,
                quantification) |>
  saveRDS(file = opt$out_path)

cat("  saved.\n")

q(status=0)

