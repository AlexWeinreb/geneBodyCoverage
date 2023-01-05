#' Select last element of vector
#'
#' @param x Vector
#'
#' @return The last element of `x`
last_element <- function(x){
  x[length(x)]
}


#' Convert coordinates from spliced to unspliced
#'
#' @param coords_in_spliced Bins positions in spliced coordinates (where introns count 0 bp)
#' @param tx_strand Genomic strand of this transcript
#' @param exons_lengths,exons_starts Exons lengths and starts
#' @param opt List containing endedness (do we count from 3' or from 5')
#'
#' @return Bins positions in unspliced coordinates (i.e. the transcript starts at 0, but we add space for the introns)
#'
coords_to_unspliced <- function(coords_in_spliced, tx_strand, exons_lengths,
                                exons_starts, opt){

  if((tx_strand == "+" && opt$endedness == 3L) |
     (tx_strand == "-" && opt$endedness == 5L)){
    # count from right to left (where BED is left to right)
    tx_end <- last_element(exons_lengths) + last_element(exons_starts)
    exons_starts <- tx_end - rev(exons_lengths + exons_starts)
    exons_lengths <- rev(exons_lengths)
  }

  cumlength <- cumsum(exons_lengths)
  which_exon <- purrr::map_int(coords_in_spliced,
                        ~ min(which(.x <= cumlength)))

  exons_starts[which_exon] -
    c(0,cumlength)[which_exon] +
    coords_in_spliced
}


#' Convert unspliced to genomic coordinates
#'
#' @param tx_st Genomic position of transcript start
#' @param tx_end Genomic position of transcript end
#' @param bins Bins positions in unspliced coordinates
#' @param tx_strand Genomic strand of transcript
#' @param opt List containing endedness (do we count from 3' or from 5')
#'
#' @return Bins positions in genomic coordinates
coords_to_genomic <- function(tx_st, tx_end, bins, tx_strand, opt){
  `if`((tx_strand == "+" && opt$endedness == 5L) |
         (tx_strand == "-" && opt$endedness == 3L),
       tx_st + bins + 1,
       tx_end - bins)
}


#' Coordinates to GRanges
#'
#' @param tx_chr,tx_strand,tx_start Transcript genomic position
#' @param bins Bins positions in genomic coordinates
#' @param tx_ex_starts,tx_ex_len Starts and lengths of the exons
#' @param opt List containing endedness (do we count from 3' or from 5')
#'
#' @return A list of GRranges objects, each element of the list is one bin, There might be several Ranges within a single element if that bins spans several exons (making holes for introns)
#'
#' @importFrom methods as
coords_as_bin_granges <- function(tx_chr, bins, tx_strand,
                              tx_start, tx_ex_starts, tx_ex_len,
                              opt){

  if((tx_strand == "+" && opt$endedness == 5L) |
     (tx_strand == "-" && opt$endedness == 3L)){
    # left to right
    tx_gr <- GenomicRanges::GRanges(seqnames = tx_chr,
                     ranges = IRanges::IRanges(start = bins[-length(bins)] + 1,
                                      end = bins[-1]),
                     strand = tx_strand)
  } else{
    tx_gr <- GenomicRanges::GRanges(seqnames = tx_chr,
                     ranges = IRanges::IRanges(start = bins[-1] + 1,
                                      end = bins[-length(bins)]),
                     strand = tx_strand)
  }


  # make holes for introns
  exons_gr <- GenomicRanges::GRanges(seqnames = tx_chr,
                      ranges = IRanges::IRanges(start = tx_start +
                                         tx_ex_starts + 1,
                                       width = tx_ex_len),
                      strand = tx_strand)

  tx_gr |>
    as("GRangesList") |>
    sapply(\(.x) GenomicRanges::intersect(.x, exons_gr))
}



#' Convert bins to a list of GRanges objects in genomic coordinates
#'
#' As input, we take the breaks (the bins limits), as output, we get GRanges that cover the bins (excluding exons).
#' Note this is different from `breaks_to_granges()` which returns a small GRange around each breakpoint.
#'
#' @param tx_chr,tx_strand,tx_start,tx_end,spliced_tx_width Genomic position of the transcript
#' @param exons_lengths,exons_starts Structure of the transcript
#' @param bins Breaks for the bins
#' @param opt List containing endedness (do we count from 3' or from 5')
#'
#' @return A list of GRranges objects, each element of the list is one bin, There might be several Ranges within a single element if that bins spans several exons (making holes for introns)
#' @export
#'
#' @examples
#' bins_to_granges(tx_chr = "I", tx_strand = "+",
#'   tx_start = 100, tx_end = 200, spliced_tx_width = 60,
#'   exons_lengths = c(10,50),exons_starts = c(0,150),
#'   bins = c(0,5,15,30),
#'   opt = list(endedness = 5L))
bins_to_granges <- function(tx_chr, tx_strand, tx_start, tx_end, spliced_tx_width,
                            exons_lengths, exons_starts,
                            bins,
                            opt){

  bins_in_spliced <- bins[bins <= spliced_tx_width]

  bins_in_unspliced <- coords_to_unspliced(bins_in_spliced, tx_strand, exons_lengths,
                                  exons_starts, opt)

  bins_in_genomic <- coords_to_genomic(tx_start, tx_end, bins_in_unspliced, tx_strand, opt)

  bins_as_granges <- coords_as_bin_granges(tx_chr, bins_in_genomic, tx_strand,
                                       tx_start, exons_starts, exons_lengths,
                                       opt)

  bins_as_granges
}





#' Convert bins to a list of GRanges objects in genomic coordinates
#'
#' As input, we take the breakpoints, as output, we get GRanges that cover these points (1 bp-wide each).
#'  Note this is different from `bins_to_granges()` which returns the bins between the breakpoints.
#'
#' @param tx_chr,tx_strand,tx_start,tx_end,spliced_tx_width Genomic position of the transcript
#' @param exons_lengths,exons_starts Structure of the transcript
#' @param breakpoints Coordinates of the breakpoints (in the spliced transcript)
#' @param opt List containing endedness (do we count from 3' or from 5')
#'
#' @return A list of GRranges objects, each element of the list is one breakpoint.
#' @export
#'
#' @examples
#' breaks_to_genomic(tx_chr = "I", tx_strand = "+",
#'   tx_start = 100, tx_end = 200, spliced_tx_width = 60,
#'   exons_lengths = c(10,50), exons_starts = c(0,150),
#'   breakpoints = c(0,5,15,30),
#'   opt = list(endedness = 5L))
breaks_to_genomic <- function(tx_chr, tx_strand, tx_start, tx_end, spliced_tx_width,
                            exons_lengths, exons_starts,
                            breakpoints,
                            opt){

  breaks_in_spliced <- breakpoints[breakpoints <= spliced_tx_width]

  breaks_in_unspliced <- coords_to_unspliced(breaks_in_spliced, tx_strand, exons_lengths,
                                           exons_starts, opt)

  breaks_in_genomic <- coords_to_genomic(tx_start, tx_end, breaks_in_unspliced, tx_strand, opt)

  breaks_in_genomic
}




