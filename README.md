
# geneBodyCoverage


Compute transcript coverage from a BAM file and transcripts annotations.

## Installation

Install the development version of geneBodyCoverage with:

``` r
remotes::install_github("AlexWeinreb/geneBodyCoverage")
```

## Inputs

The annotation is in the same format as RSeQC, i.e. a BED file generated from a genePred file. This can be obtained from a GTF using the [kentUtils toolset](https://hgdownload.soe.ucsc.edu/admin/exe/):

```bash
gtfToGenePred /path/canonical_geneset.gtf /path/canonical_geneset.genePred

genePredToBed /path/canonical_geneset.genePred /path/canonical_geneset.genePred.bed
```

The BAM file containing RNA-Seq reads should be sorted and indexed.


## Bins Mode and Breaks Mode

There are two ways to count the reads from a set of breakpoints: first, we might take each of these breakpoints and count the number of reads that align on top of it (Breaks Mode). Alternatively, we can take the breakpoints as limits of bins, and count the number of reads that fall anywhere within the bin (Bins Mode). RSeQC seems to use the first approach (i.e. Breaks Mode).

Here is an illustrative fictional gene, with 2 exons, on the + strand, covered by 3 reads:

```{text}
coordinates        0    5         10   15
exon structure     ##########------########
breakpoints        |    |          |    |  
reads                 >>>>
                        >>>>>------>>>>>>>
                           >>------>
```

The table of number of reads per breakpoint is:

  breakpoint   number of reads
------------ -----------------
         0                   0
         5                   2
        10                   2
        15                   1

Whereas the number of reads falling within the bins defined by the breakpoints are:

  bin start     bin end  number of reads
----------- ----------- ----------------
        1           5                  2
        6          10                  3
       11          15                  1



Using the package within R, we can use two sets of functions to quantify what we wish:

To count reads on top of the breakpoints, first, the set of breakpoints needs to be transformed into genomic coordinates with `breaks_to_genomic()`. Second, we can count the number of reads on each of these breakpoints with `breaks_coverage()`.

To count reads within bins, first, the set of breakpoints is transformed into bins (as GRanges objects) with the function `bins_to_granges()`. Second, we can count the number of reads that fall within each bin using `bins_coverage()`.






## Wrapping script

For use on a cluster or as part as non-R pipelines, a script is provided within the package, with the name `genebody_coverage.R`. It can be called directly from the command line, for example with:

```bash
Rscript.exe genebody_coverage.R \
  --bam_path=../tests/testthat/fixtures/simple.bam \
  --bed_path=../tests/testthat/fixtures/c_elegans.PRJNA13758.WS281_filtered.bed \
  --endedness=5 \
  --chrom_list=elegans \
  --out_path=c:/Users/Alexis\ Weinreb/Projects/tests/tmp/gbc_out/bins_output.rds
```

Use switch `--help` to get more details on the available options:

```bash
$ Rscript.exe genebody_coverage.R --help
Usage: inst/genebody_coverage.R [-[-bam_path|b] <character>] [-[-bed_path|r] <character>] [-[-endedness|e] <integer>] [-[-within_bin|w]] [-[-delimiter|d] [<character>]] [-[-chrom_list|c] [<character>]] [-[-out
_path|o] <character>] [-[-nproc|p] [<integer>]] [-[-future_mem|m] [<integer>]] [-[-force|f]] [-[-min_length|l] [<integer>]] [-[-help|h]]
    -b|--bam_path      Path to BAM file
    -r|--bed_path      Path to reference BED12 file
    -e|--endedness     From which end to start counting (3 or 5)
    -w|--within_bin    Count reads within bins (default: count reads on the breaks)
    -d|--delimiter     Delimiter of the bed file (e.g. space or tab, default: tab)
    -c|--chrom_list    Name of the list of chromosomes, or blank (elegans, hg38 accepted)
    -o|--out_path      Path to save output (must not exists, if not -f)
    -p|--nproc         Number of threads
    -m|--future_mem    Max memory for future globals
    -f|--force         Set to overwrite existing output
    -l|--min_length    Minimum transcript length, any transcript shorter will be ignored (default: 100 bp)
    -h|--help          Print this help
```

Note in particular that if `--within_bin` is set, the quantification will be in Bins Mode. If this switch is not set, the quantification will be in Breaks Mode.

This script is installed in the root directory of the package, you can use `path.package("geneBodyCoverage")` to find the path. In the source package, you can find it in the `inst/` directory. To use it in your workflows, you may want to copy (or symlink) it in a convenient location, e.g. `~/bin`.



## Note on performance

In Bins Mode, the script runs reasonably fast for a few million reads, but may struggle when the size of the genome and number of reads becomes large. It may be useful to subset a few chromosomes.

One way to subset the bed file to keep only chromosomes 5, 13 and 18:
```bash
grep -e "^5 " -e "^13 " -e "18 " Homo_sapiens.gp.bed > Homo_sapiens_subset.gp.bed
```

And using [samtools view](http://www.htslib.org/doc/samtools-view.html) to subset the bam file, indicating the chromosomes to keep as regions:

```bash
samtools view \
  -bo subset_chrom5_13_18.bam \
  full.bam \
  5 13 18
```


## See also

The function [genebody_coverage.py](https://rseqc.sourceforge.net/#genebody-coverage-py) from the RSeQC package, which computes coverage on percentiles of the transcripts length. [Written in Python](https://github.com/MonashBioinformaticsPlatform/RSeQC/blob/master/rseqc/modules/geneBody_coverage.py).
