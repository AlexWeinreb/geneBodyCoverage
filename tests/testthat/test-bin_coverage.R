
# Initializations ----

bam_content <- GenomicAlignments::readGAlignments(file = test_path("fixtures", "simple.bam"),
                                                  index = test_path("fixtures", "simple.bam"))

## Prepare bins

bins <- c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,
          1700,1800,1900,2000,2100,2200,2300,2400,2500,2700,2900,3100,3400,3700,
          4000,4500,5000,5500,6000,6500,7000,8000,9000,10000,Inf)

gene_struct <- readr::read_delim(test_path("fixtures","c_elegans.PRJNA13758.WS281_filtered.bed"),
                                 delim = "\t",
                                 col_names = c("chr","transcript_start","transcript_end","transcript_name","score",
                                               "strand","thick_start","thick_end","rgb",
                                               "nb_exons","exons_lengths","exons_starts"),
                                 col_types = readr::cols(
                                   chr = readr::col_factor(as.roman(1:5) |>
                                                             as.character() |>
                                                             c("X","MtDNA")),      #1
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
                  purrr::map(as.integer))

prepare_transcript <- function(tx_name, endedness){
  opt <- list(endedness = endedness)
  xx <- gene_struct |>
    dplyr::filter(transcript_name == tx_name) |>
    dplyr::mutate(tx_width = purrr::map_int(exons_lengths, sum)) |>
    as.list()

  bins_in_spliced_coords <- purrr::map(xx$tx_width,
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

  bins_gr <- coords_as_granges(xx$chr,
                               bins_in_genomic_coords,
                               xx$strand,
                               xx$transcript_start,
                               xx$exons_starts[[1]],
                               xx$exons_lengths[[1]],
                               opt)
  bins_gr
}

# Tests ----


#~ MTCE.3.1 ----

bins_coverage(prepare_transcript(tx_name = "MTCE.3.1",endedness =  5L),
              "MtDNA", "+", bam_content) |>
  expect_identical(c(2L,1L,1L,0L))

bins_coverage(prepare_transcript(tx_name = "MTCE.3.1",endedness =  3L),
              "MtDNA", "+", bam_content) |>
  expect_identical(c(0L,1L,1L,2L))




#~ B0348.5b.1 ----

bins_coverage(prepare_transcript(tx_name = "B0348.5b.1",endedness =  5L),
              "V", "+",
              bam_content) |>
  expect_identical(c(0L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L))

bins_coverage(prepare_transcript(tx_name = "B0348.5b.1",endedness =  3L),
              "V", "+",
              bam_content) |>
  expect_identical(c(0L,0L,0L,0L,0L,1L,1L,1L,2L,1L))



#~ Y74C9A.6 ----

bins_coverage(prepare_transcript(tx_name = "Y74C9A.6",endedness =  5L),
              "I", "-",
              bam_content) |>
  expect_identical(c(2L))

bins_coverage(prepare_transcript(tx_name = "Y74C9A.6",endedness =  3L),
              "I", "-",
              bam_content) |>
  expect_identical(c(3L))



#~ F23F1.10.1 ----

bins_coverage(prepare_transcript(tx_name = "F23F1.10.1",endedness =  5L),
              "II", "-",
              bam_content) |>
  expect_identical(c(1L, 2L, 0L))

bins_coverage(prepare_transcript(tx_name = "F23F1.10.1",endedness =  3L),
              "II", "-",
              bam_content) |>
  expect_identical(c(0L, 2L, 1L))



# Different bins ----

bins <- c(0, 10, 20,60, 150, 300, 3000)

bins_coverage(prepare_transcript(tx_name = "F23F1.10.1",endedness =  5L),
              "II", "-",
              bam_content) |>
  expect_identical(c(0L, 0L, 0L, 1L, 2L))

bins_coverage(prepare_transcript(tx_name = "F23F1.10.1",endedness =  3L),
              "II", "-",
              bam_content) |>
  expect_identical(c(0L, 0L, 0L, 0L, 2L))



