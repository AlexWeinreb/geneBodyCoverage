last_element <- function(x) x[length(x)]


coords_to_unspliced <- function(coords_in_spliced, strand, exons_lengths, exons_starts){

  if((strand == "+" && opt$endedness == 3L) |
     (strand == "-" && opt$endedness == 5L)){
    # count from right to left (where BED is left to right)
    tx_end <- last_element(exons_lengths) + last_element(exons_starts)
    exons_starts <- tx_end - rev(exons_lengths + exons_starts)
    exons_lengths <- rev(exons_lengths)
  }

  cumlength <- cumsum(exons_lengths)
  which_exon <- map_int(coords_in_spliced,
                        ~ min(which(.x <= cumlength)))

  exons_starts[which_exon] -
    c(0,cumlength)[which_exon] +
    coords_in_spliced
}


coords_to_genomic <- function(.tx_st, .tx_end, .bins, .strand){
  `if`((.strand == "+" && opt$endedness == 5L) |
         (.strand == "-" && opt$endedness == 3L),
       .tx_st + .bins + 1,
       .tx_end - .bins)
}


coords_as_granges <- function(tx_chr, tx_bins, tx_strand,
                              tx_start, tx_ex_starts, tx_ex_len){

  if((tx_strand == "+" && opt$endedness == 5L) |
     (tx_strand == "-" && opt$endedness == 3L)){
    # left to right
    tx_gr <- GRanges(seqnames = tx_chr,
                     ranges = IRanges(start = tx_bins[-length(tx_bins)] + 1,
                                      end = tx_bins[-1]),
                     strand = tx_strand)
  } else{
    tx_gr <- GRanges(seqnames = tx_chr,
                     ranges = IRanges(start = tx_bins[-1] + 1,
                                      end = tx_bins[-length(tx_bins)]),
                     strand = tx_strand)
  }


  # make holes for introns
  exons_gr <- GRanges(seqnames = tx_chr,
                      ranges = IRanges(start = tx_start +
                                         tx_ex_starts + 1,
                                       width = tx_ex_len),
                      strand = tx_strand)

  tx_gr |>
    as("GRangesList") |>
    sapply(\(.x) intersect(.x, exons_gr))
}







