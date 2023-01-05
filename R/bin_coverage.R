


countUniqueOverlappers <- function(x, y){
  GenomicAlignments::findOverlaps(x, y) |>
    S4Vectors::to() |>
    unique() |>
    length()
}



get_bins_coverage <- function(bins_gr, chr, strand,bam_content){


  # Count reads overlapping each bin

  # Note benchmarking:
  # without filtering, takes 21s (alloc 18.6GB)
  # filter on chromosome and strand, 2s (alloc 1.5 GB)
  # filter on chrom, strand, and start/end, 4.7s (alloc 700MB)
  bam_subset <- bam_content[GenomeInfoDb::seqnames(bam_content) == as.character(chr) &
                              BiocGenerics::strand(bam_content) == as.character(strand)]

  sapply(bins_gr,
         \(.x) countUniqueOverlappers(.x,
                                      bam_subset))
}



#' Bins coverage
#'
#' @param bins_gr Bins positions in genomic coordinates as a GRanges
#' @param chr,strand Genomic position of the transcript (to filter the BAM)
#' @param bam_content Content of the BAM file containing reads
#'
#' @return Count of reads in bins
#' @export
bins_coverage <- function(bins_gr, chr, strand, bam_content){
  tryCatch(get_bins_coverage(bins_gr, chr, strand, bam_content),
           error = \(.x) NA)
}



#' Breakpoints coverage
#'
#' @param breaks_gr Breakpoint positions as a GRanges
#' @param reads Content of the BAM file containing reads
#'
#' @return Count of reads on breakpoints
#' @export
breaks_coverage <- function(breaks_gr, reads){
  GenomicRanges::countOverlaps(breaks_gr, reads)
}

