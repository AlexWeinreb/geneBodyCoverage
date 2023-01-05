
# library(testthat)


test_that("Breaks from 5' gives expected result",
          {
            tmp <- file.path(tempdir(), "out.tsv")
            on.exit(file.remove(tmp))

            # Breaks, 5' ----
            cmd <- paste0("Rscript ", shQuote(system.file("genebody_coverage.R", package = "geneBodyCoverage"))," ",
                          "--bam_path=fixtures/simple.bam ",
                          "--bed_path=fixtures/c_elegans.PRJNA13758.WS281_filtered.bed ",
                          "--endedness=5 ",
                          "--chrom_list=elegans ",
                          "--out_path=",tmp," ",
                          "--nproc=1 ")
            cmd
            out <- system(cmd)

            system("Rscript C:\\Users\\Alexis Weinreb\\Projects\\other\\geneBodyCoverage\\inst\\genebody_coverage.R --bam_path=fixtures/simple.bam --bed_path=fixtures/c_elegans.PRJNA13758.WS281_filtered.bed --endedness=5 --chrom_list=elegans --out_path=C:\\Users\\ALEXIS~1\\AppData\\Local\\Temp\\RtmpcVyqI8/out.tsv --nproc=1 ")
            system("Rscript C:\\Users\\Alexis\ Weinreb\\Projects\\other\\geneBodyCoverage\\inst\\genebody_coverage.R --bam_path=fixtures/simple.bam --bed_path=fixtures/c_elegans.PRJNA13758.WS281_filtered.bed --endedness=5 --chrom_list=elegans --out_path=C:\\Users\\ALEXIS~1\\AppData\\Local\\Temp\\RtmpcVyqI8/out.tsv --nproc=1 ")
            expect_identical(out, 0L)

            expect_identical(readr::read_tsv(tmp, skip = 1L, show_col_types = FALSE),
                             readr::read_tsv("fixtures/results_script/breaks_5prime.tsv", skip = 1L, show_col_types = FALSE))

          })


test_that("Breaks from 3' gives expected result",
          {
            tmp <- file.path(tempdir(), "out.tsv")
            on.exit(file.remove(tmp))

            # Breaks, 3' ----
            cmd <- paste0("Rscript ", shQuote(system.file("genebody_coverage.R", package = "geneBodyCoverage"))," ",
                          "--bam_path=fixtures/simple.bam ",
                          "--bed_path=fixtures/c_elegans.PRJNA13758.WS281_filtered.bed ",
                          "--endedness=3 ",
                          "--chrom_list=elegans ",
                          "--out_path=",tmp," ",
                          "--nproc=1 ")
            cmd
            out <- system(cmd)

            expect_identical(out, 0L)

            expect_identical(readr::read_tsv(tmp, skip = 1L, show_col_types = FALSE),
                             readr::read_tsv("fixtures/results_script/breaks_3prime.tsv", skip = 1L, show_col_types = FALSE))

          })


test_that("Bins from 5' gives expected result",
          {
            tmp <- file.path(tempdir(), "out.tsv")
            on.exit(file.remove(tmp))

            # Bins, 5' ----
            cmd <- paste0("Rscript ", shQuote(system.file("genebody_coverage.R", package = "geneBodyCoverage"))," ",
                          "--bam_path=fixtures/simple.bam ",
                          "--bed_path=fixtures/c_elegans.PRJNA13758.WS281_filtered.bed ",
                          "--endedness=5 ",
                          "--within_bin ",
                          "--chrom_list=elegans ",
                          "--out_path=",tmp," ",
                          "--nproc=1 ")
            cmd
            out <- system(cmd)

            expect_identical(out, 0L)

            expect_identical(readr::read_tsv(tmp, skip = 1L, show_col_types = FALSE),
                             readr::read_tsv("fixtures/results_script/bins_5prime.tsv", skip = 1L, show_col_types = FALSE))

          })






test_that("Bins from 3' gives expected result",
          {
            tmp <- file.path(tempdir(), "out.tsv")
            on.exit(file.remove(tmp))

            # Bins, 3' ----
            cmd <- paste0("Rscript ", shQuote(system.file("genebody_coverage.R", package = "geneBodyCoverage"))," ",
                          "--bam_path=fixtures/simple.bam ",
                          "--bed_path=fixtures/c_elegans.PRJNA13758.WS281_filtered.bed ",
                          "--endedness=3 ",
                          "--within_bin ",
                          "--chrom_list=elegans ",
                          "--out_path=",tmp," ",
                          "--nproc=1 ")
            cmd
            out <- system(cmd)

            expect_identical(out, 0L)

            expect_identical(readr::read_tsv(tmp, skip = 1L, show_col_types = FALSE),
                             readr::read_tsv("fixtures/results_script/bins_3prime.tsv", skip = 1L, show_col_types = FALSE))

          })




