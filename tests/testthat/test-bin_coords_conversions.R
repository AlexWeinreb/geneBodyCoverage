


# devtools::load_all()
# library(testthat)

bins <- c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,
          1700,1800,1900,2000,2100,2200,2300,2400,2500,2700,2900,3100,3400,3700,
          4000,4500,5000,5500,6000,6500,7000,8000,9000,10000,Inf)

chroms <- as.roman(1:5) |>
  as.character() |>
  c("X","MtDNA")

gene_struct <- readr::read_delim(test_path("fixtures","c_elegans.PRJNA13758.WS281_filtered.bed"),
                          delim = "\t",
                          col_names = c("chr","transcript_start","transcript_end","transcript_name","score",
                                        "strand","thick_start","thick_end","rgb",
                                        "nb_exons","exons_lengths","exons_starts"),
                          col_types = readr::cols(
                            chr = readr::col_factor(chroms),      #1
                            transcript_start = readr::col_integer(),        #2 on the + strand (even if tx on -)
                            transcript_end = readr::col_integer(),          #3
                            transcript_name = readr::col_character(),       #4
                            score = readr::col_factor("0"),                 #5 ignore
                            strand = readr::col_factor(c("+","-")),         #6
                            thick_start = readr::col_integer(),             #7 ignore
                            thick_end = readr::col_integer(),               #8 ignore
                            rgb = readr::col_factor("0"),                   #9 ignore
                            nb_exons = readr::col_integer(),                #10
                            exons_lengths = readr::col_character(),         #11
                            exons_starts = readr::col_character()           #12 on the + strand (even if tx on -)
                          )) |>
  dplyr::mutate(exons_lengths = strsplit(exons_lengths, ",") |>
                  purrr::map(as.integer),
                exons_starts = strsplit(exons_starts, ",") |>
                  purrr::map(as.integer),
                spliced_tx_width = purrr::map_int(exons_lengths, sum))





# gene-by-gene rather than function by function

# Run all the transcript bins conversions
prepare_transcript <- function(tx_name, endedness){
  opt <- list(endedness = endedness)
  xx <- gene_struct |>
    dplyr::filter(transcript_name == tx_name) |>
    as.list()


  bins_in_spliced_coords <- purrr::map(xx$spliced_tx_width,
                                       ~bins[bins <= .x])[[1]]

  bins_in_unspliced_coords <- coords_to_unspliced(bins_in_spliced_coords,
                                                  xx$strand,
                                                  xx$exons_lengths[[1]],
                                                  xx$exons_starts[[1]],
                                                  opt)

  bins_in_genomic_coords <- coords_to_genomic(xx$transcript_start,
                                              xx$transcript_end,
                                              bins_in_unspliced_coords,
                                              xx$strand,
                                              opt)

  bins_gr <- coords_as_bin_granges(xx$chr,
                               bins_in_genomic_coords,
                               xx$strand,
                               xx$transcript_start,
                               xx$exons_starts[[1]],
                               xx$exons_lengths[[1]],
                               opt)
  breaks_gr <- coords_as_min_granges(xx$chr,
                                     xx$strand,
                                     bins_in_genomic_coords)

  list(bins_in_spliced_coords = bins_in_spliced_coords,
       bins_in_unspliced_coords = bins_in_unspliced_coords,
       bins_in_genomic_coords = bins_in_genomic_coords,
       bins_gr = bins_gr,
       breaks_gr = breaks_gr)
}






# Test genes ----

#~ MTCE.3.1, 5' ----
local({

  res <- prepare_transcript(tx_name = "MTCE.3.1",endedness =  5L)


  expect_identical(res$bins_in_spliced_coords,
                   c(0,100,200,300,400))
  expect_identical(res$bins_in_unspliced_coords,
                   c(0,100,200,300,400))
  expect_identical(res$bins_in_genomic_coords,
                   c(113,213,313,413,513))


  my_test_bin_gr <- purrr::map2(c(114,214,314,414),
                     c(213,313,413,513),
                     ~ GenomicRanges::GRanges(seqnames = factor("MtDNA", chroms),
                              ranges = IRanges::IRanges(.x, .y),
                              strand = "+"))

  expect_identical(res$bins_gr,
                   my_test_bin_gr)

  gr_full <- bins_to_granges(tx_chr = factor("MtDNA", levels = chroms),
                             tx_strand = factor("+", levels=c("+","-")),
                             tx_start = 112, tx_end = 549, spliced_tx_width = 437,
                             exons_lengths = c(437), exons_starts = c(0),
                             bins = bins, opt = list(endedness = 5L))

  expect_identical(gr_full, my_test_bin_gr)


  my_test_min_gr <- GenomicRanges::GRanges(seqnames = factor("MtDNA", chroms),
                                           ranges = IRanges::IRanges(start = c(113,213,313,413,513),
                                                                     width = 1L),
                                           strand = "+")
  expect_identical(res$breaks_gr, my_test_min_gr)

  full_min_gr <- breaks_to_granges(tx_chr = factor("MtDNA", levels = chroms),
                    tx_strand = factor("+", levels=c("+","-")),
                    tx_start = 112, tx_end = 549, spliced_tx_width = 437,
                    exons_lengths = c(437), exons_starts = c(0),
                    breakpoints = bins, opt = list(endedness = 5L))

  expect_identical(full_min_gr, my_test_min_gr)
})








#~ MTCE.3.1, 3' ----

local({
  res <- prepare_transcript(tx_name = "MTCE.3.1",endedness =  3L)

  expect_identical(res$bins_in_spliced_coords,
                   c(0,100,200,300,400))
  expect_identical(res$bins_in_unspliced_coords,
                   c(0,100,200,300,400))
  expect_identical(res$bins_in_genomic_coords,
                   c(549,449,349,249,149))


  my_test_bin_gr <- purrr::map2(c(450,350,250,150),
                     c(549,449,349,249),
                     ~GenomicRanges::GRanges(seqnames = factor("MtDNA", levels = chroms),
                              ranges = IRanges::IRanges(.x, .y),
                              strand = "+"))



  expect_identical(res$bins_gr,
                   my_test_bin_gr)

  gr_full <- bins_to_granges(tx_chr = factor("MtDNA", levels = chroms),
                             tx_strand = factor("+", levels=c("+","-")),
                             tx_start = 112, tx_end = 549,
                             exons_lengths = c(437), exons_starts = c(0), spliced_tx_width = 437,
                             bins = bins, opt = list(endedness = 3L))

  expect_identical(gr_full, my_test_bin_gr)

  my_test_min_gr <- GenomicRanges::GRanges(seqnames = factor("MtDNA", chroms),
                                           ranges = IRanges::IRanges(start = c(549,449,349,249,149),
                                                                     width = 1L),
                                           strand = "+")
  expect_identical(res$breaks_gr, my_test_min_gr)

  full_min_gr <- breaks_to_granges(tx_chr = factor("MtDNA", levels = chroms),
                                   tx_strand = factor("+", levels=c("+","-")),
                                   tx_start = 112, tx_end = 549,
                                   exons_lengths = c(437), exons_starts = c(0), spliced_tx_width = 437,
                                   breakpoints = bins, opt = list(endedness = 3L))

  expect_identical(full_min_gr, my_test_min_gr)
})






#~ B0348.5b.1, 5' ----
local({
  res <- prepare_transcript(tx_name = "B0348.5b.1",endedness =  5L)




  expect_identical(res$bins_in_spliced_coords,
                   c(0,100,200,300,400,500,600,700,800,900,1000))
  expect_identical(res$bins_in_unspliced_coords,
                   c(0,100,200,300,400,557,657,757,857,957,1057))
  expect_identical(res$bins_in_genomic_coords,
                   c(5536,5636,5736,5836,5936,6093,6193,6293,6393,6493,6593))


  my_test_bin_gr <- purrr::map2(c(5537,5637,5737,5837,5937,6094,6194,6294,6394,6494),
                     c(5636,5736,5836,5936,6093,6193,6293,6393,6493,6593),
                     ~GenomicRanges::GRanges(seqnames = factor("V", levels = chroms),
                              ranges = IRanges::IRanges(.x, .y),
                              strand = "+"))

  my_test_bin_gr[[5]] <- GenomicRanges::GRanges(seqnames = factor("V", levels = chroms),
          ranges = IRanges::IRanges(c(5937,6024),c(5966,6093)),
          strand = "+")


  expect_identical(res$bins_gr,
                   my_test_bin_gr)


  gr_full <- bins_to_granges(tx_chr = factor("V", levels = chroms),
                             tx_strand = factor("+", levels=c("+","-")),
                             tx_start = 5535, tx_end = 6634, spliced_tx_width = 1042,
                             exons_lengths = c(431, 611), exons_starts = c(0, 488),
                             bins = bins, opt = list(endedness = 5L))

  expect_identical(gr_full, my_test_bin_gr)

  my_test_min_gr <- GenomicRanges::GRanges(seqnames = factor("V", chroms),
                                           ranges = IRanges::IRanges(start = c(5536,5636,5736,5836,5936,6093,6193,6293,6393,6493,6593),
                                                                     width = 1L),
                                           strand = "+")
  expect_identical(res$breaks_gr, my_test_min_gr)

  full_min_gr <- breaks_to_granges(tx_chr = factor("V", levels = chroms),
                                   tx_strand = factor("+", levels=c("+","-")),
                                   tx_start = 5535, tx_end = 6634, spliced_tx_width = 1042,
                                   exons_lengths = c(431, 611), exons_starts = c(0, 488),
                                   breakpoints = bins, opt = list(endedness = 5L))

  expect_identical(full_min_gr, my_test_min_gr)
})








#~ B0348.5b.1, 3' ----

local({
  res <- prepare_transcript(tx_name = "B0348.5b.1",endedness =  3L)


  expect_identical(res$bins_in_spliced_coords,
                   c(0,100,200,300,400,500,600,700,800,900,1000))
  expect_identical(res$bins_in_unspliced_coords,
                   c(0,100,200,300,400,500,600,757,857,957,1057))
  expect_identical(res$bins_in_genomic_coords,
                   c(6634,6534,6434,6334,6234,6134,6034,5877,5777,5677,5577))


  my_test_bin_gr <- purrr::map2(c(6535,6435,6335,6235,6135,6035,5878,5778,5678,5578),
                     c(6634,6534,6434,6334,6234,6134,6034,5877,5777,5677),
                     ~GenomicRanges::GRanges(seqnames = factor("V", levels = chroms),
                              ranges = IRanges::IRanges(.x, .y),
                              strand = "+"))

  my_test_bin_gr[[7]] <- GenomicRanges::GRanges(seqnames = factor("V", levels = chroms),
                             ranges = IRanges::IRanges(c(5878,6024),c(5966,6034)),
                             strand = "+")


  expect_identical(res$bins_gr,
                   my_test_bin_gr)

  gr_full <- bins_to_granges(tx_chr = factor("V", levels = chroms),
                             tx_strand = factor("+", levels=c("+","-")),
                             tx_start = 5535, tx_end = 6634, spliced_tx_width = 1042,
                             exons_lengths = c(431, 611), exons_starts = c(0, 488),
                             bins = bins, opt = list(endedness = 3L))

  expect_identical(gr_full, my_test_bin_gr)

  my_test_min_gr <- GenomicRanges::GRanges(seqnames = factor("V", chroms),
                                           ranges = IRanges::IRanges(start = c(6634,6534,6434,6334,6234,6134,6034,5877,5777,5677,5577),
                                                                     width = 1L),
                                           strand = "+")
  expect_identical(res$breaks_gr, my_test_min_gr)

  full_min_gr <- breaks_to_granges(tx_chr = factor("V", levels = chroms),
                                   tx_strand = factor("+", levels=c("+","-")),
                                   tx_start = 5535, tx_end = 6634, spliced_tx_width = 1042,
                                   exons_lengths = c(431, 611), exons_starts = c(0, 488),
                                   breakpoints = bins, opt = list(endedness = 3L))

  expect_identical(full_min_gr, my_test_min_gr)
})







#~ Y74C9A.6, 5' ----
local({
  res <- prepare_transcript(tx_name = "Y74C9A.6",endedness =  5L)




  expect_identical(res$bins_in_spliced_coords,
                   c(0,100))
  expect_identical(res$bins_in_unspliced_coords,
                   c(0,100))
  expect_identical(res$bins_in_genomic_coords,
                   c(3909,3809))


  my_test_bin_gr <- purrr::map2(c(3810),
                     c(3909),
                     ~GenomicRanges::GRanges(seqnames = factor("I", levels = chroms),
                              ranges = IRanges::IRanges(.x, .y),
                              strand = "-"))


  expect_identical(res$bins_gr,
                   my_test_bin_gr)

  gr_full <- bins_to_granges(tx_chr = factor("I", levels = chroms),
                             tx_strand = factor("-", levels=c("+","-")),
                             tx_start = 3746, tx_end = 3909, spliced_tx_width = 163,
                             exons_lengths = c(163), exons_starts = c(0),
                             bins = bins, opt = list(endedness = 5L))

  expect_identical(gr_full, my_test_bin_gr)

  my_test_min_gr <- GenomicRanges::GRanges(seqnames = factor("I", chroms),
                                           ranges = IRanges::IRanges(start = c(3909,3809),
                                                                     width = 1L),
                                           strand = "-")
  expect_identical(res$breaks_gr, my_test_min_gr)

  full_min_gr <- breaks_to_granges(tx_chr = factor("I", levels = chroms),
                                   tx_strand = factor("-", levels=c("+","-")),
                                   tx_start = 3746, tx_end = 3909, spliced_tx_width = 163,
                                   exons_lengths = c(163), exons_starts = c(0),
                                   breakpoints = bins, opt = list(endedness = 5L))

  expect_identical(full_min_gr, my_test_min_gr)
})








#~ Y74C9A.6, 3' ----

local({
  res <- prepare_transcript(tx_name = "Y74C9A.6",endedness =  3L)


  expect_identical(res$bins_in_spliced_coords,
                   c(0,100))
  expect_identical(res$bins_in_unspliced_coords,
                   c(0,100))
  expect_identical(res$bins_in_genomic_coords,
                   c(3747,3847))


  my_test_bin_gr <- purrr::map2(c(3748),
                     c(3847),
                     ~GenomicRanges::GRanges(seqnames = factor("I", levels = chroms),
                              ranges = IRanges::IRanges(.x, .y),
                              strand = "-"))


  expect_identical(res$bins_gr,
                   my_test_bin_gr)

  gr_full <- bins_to_granges(tx_chr = factor("I", levels = chroms),
                             tx_strand = factor("-", levels=c("+","-")),
                             tx_start = 3746, tx_end = 3909, spliced_tx_width = 163,
                             exons_lengths = c(163), exons_starts = c(0),
                             bins = bins, opt = list(endedness = 3L))

  expect_identical(gr_full, my_test_bin_gr)

  my_test_min_gr <- GenomicRanges::GRanges(seqnames = factor("I", chroms),
                                           ranges = IRanges::IRanges(start = c(3747,3847),
                                                                     width = 1L),
                                           strand = "-")
  expect_identical(res$breaks_gr, my_test_min_gr)

  full_min_gr <- breaks_to_granges(tx_chr = factor("I", levels = chroms),
                                   tx_strand = factor("-", levels=c("+","-")),
                                   tx_start = 3746, tx_end = 3909, spliced_tx_width = 163,
                                   exons_lengths = c(163), exons_starts = c(0),
                                   breakpoints = bins, opt = list(endedness = 3L))

  expect_identical(full_min_gr, my_test_min_gr)
})








#~ F23F1.10.1, 5' ----
local({

  res <- prepare_transcript(tx_name = "F23F1.10.1",endedness =  5L)




  expect_identical(res$bins_in_spliced_coords,
                   c(0,100,200,300))
  expect_identical(res$bins_in_unspliced_coords,
                   c(0,100,361,461))
  expect_identical(res$bins_in_genomic_coords,
                   c(41982,41882,41621,41521))


  my_test_bin_gr <- purrr::map2(c(41883,41622,41522),
                     c(41982,41882,41621),
                     ~GenomicRanges::GRanges(seqnames = factor("II", levels = chroms),
                              ranges = IRanges::IRanges(.x, .y),
                              strand = "-"))


  my_test_bin_gr[[2]] <- GenomicRanges::GRanges(seqnames = factor("II", levels = chroms),
                             ranges = IRanges::IRanges(c(41622,41844),c(41682,41882)),
                             strand = "-")

  expect_identical(res$bins_gr,
                   my_test_bin_gr)

  gr_full <- bins_to_granges(tx_chr = factor("II", levels = chroms),
                             tx_strand = factor("-", levels=c("+","-")),
                             tx_start = 41470, tx_end = 41982, spliced_tx_width = 351,
                             exons_lengths = c(212, 139), exons_starts = c(0, 373),
                             bins = bins, opt = list(endedness = 5L))

  expect_identical(gr_full, my_test_bin_gr)

  my_test_min_gr <- GenomicRanges::GRanges(seqnames = factor("II", chroms),
                                           ranges = IRanges::IRanges(start = c(41982,41882,41621,41521),
                                                                     width = 1L),
                                           strand = "-")
  expect_identical(res$breaks_gr, my_test_min_gr)

  full_min_gr <- breaks_to_granges(tx_chr = factor("II", levels = chroms),
                                   tx_strand = factor("-", levels=c("+","-")),
                                   tx_start = 41470, tx_end = 41982, spliced_tx_width = 351,
                                   exons_lengths = c(212, 139), exons_starts = c(0, 373),
                                   breakpoints = bins, opt = list(endedness = 5L))

  expect_identical(full_min_gr, my_test_min_gr)
})








#~ F23F1.10.1, 3' ----

local({
  res <- prepare_transcript(tx_name = "F23F1.10.1",endedness =  3L)


  expect_identical(res$bins_in_spliced_coords,
                   c(0,100,200,300))
  expect_identical(res$bins_in_unspliced_coords,
                   c(0,100,200,461))
  expect_identical(res$bins_in_genomic_coords,
                   c(41471,41571,41671,41932))


  my_test_bin_gr <- purrr::map2(c(41472,41572,41672),
                     c(41571,41671,41932),
                     ~GenomicRanges::GRanges(seqnames = factor("II", levels = chroms),
                              ranges = IRanges::IRanges(.x, .y),
                              strand = "-"))


  my_test_bin_gr[[3]] <- GenomicRanges::GRanges(seqnames = factor("II", levels = chroms),
                             ranges = IRanges::IRanges(c(41672,41844),c(41682,41932)),
                             strand = "-")


  expect_identical(res$bins_gr,
                   my_test_bin_gr)


  gr_full <- bins_to_granges(tx_chr = factor("II", levels = chroms),
                             tx_strand = factor("-", levels=c("+","-")),
                             tx_start = 41470, tx_end = 41982, spliced_tx_width = 351,
                             exons_lengths = c(212, 139), exons_starts = c(0, 373),
                             bins = bins, opt = list(endedness = 3L))

  expect_identical(gr_full, my_test_bin_gr)

  my_test_min_gr <- GenomicRanges::GRanges(seqnames = factor("II", chroms),
                                           ranges = IRanges::IRanges(start = c(41471,41571,41671,41932),
                                                                     width = 1L),
                                           strand = "-")
  expect_identical(res$breaks_gr, my_test_min_gr)

  full_min_gr <- breaks_to_granges(tx_chr = factor("II", levels = chroms),
                                   tx_strand = factor("-", levels=c("+","-")),
                                   tx_start = 41470, tx_end = 41982, spliced_tx_width = 351,
                                   exons_lengths = c(212, 139), exons_starts = c(0, 373),
                                   breakpoints = bins, opt = list(endedness = 3L))

  expect_identical(full_min_gr, my_test_min_gr)
})





#~ Y66H1A.8a.1: should be too short because spliced transcript under 100 bp ----
local({

  gene_struct |>
    dplyr::filter(transcript_name == "Y66H1A.8a.1") |>
    dplyr::pull(spliced_tx_width) |>
    expect_identical(60L)


})


