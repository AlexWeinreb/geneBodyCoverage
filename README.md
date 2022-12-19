
# geneBodyCoverage


Compute transcript coverage from a BAM file and transcripts annotations.

## Installation

Writing in progress. Probably doesn't install yet.

You will be able to install the development version of geneBodyCoverage with:

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

## Note on performance

The script runs reasonably fast for a few million reads, but may struggle when the size of the genome becomes large. It may be useful to subset a few chromosomes.

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
