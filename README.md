
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rps10package-rename and add info on package

<!-- badges: start -->
<!-- badges: end -->

Package readme in progress

## Installation

You can install the development version of rps10package from
[GitHub](https://github.com/) with:

*TODO-still working on documentation*

``` r
#temporary installation method
devtools::load_all("~/rps10package")
#> ℹ Loading rps10package
#> Loading required package: dada2
#> 
#> Loading required package: Rcpp
#> 
#> Loading required package: devtools
#> 
#> Loading required package: usethis
#> 
#> Loading required package: ShortRead
#> 
#> Loading required package: BiocGenerics
#> 
#> 
#> Attaching package: 'BiocGenerics'
#> 
#> 
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> 
#> 
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
#>     table, tapply, union, unique, unsplit, which.max, which.min
#> 
#> 
#> Loading required package: BiocParallel
#> 
#> Loading required package: Biostrings
#> 
#> Loading required package: S4Vectors
#> 
#> Loading required package: stats4
#> 
#> 
#> Attaching package: 'S4Vectors'
#> 
#> 
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> 
#> 
#> Loading required package: IRanges
#> 
#> Loading required package: XVector
#> 
#> Loading required package: GenomeInfoDb
#> 
#> 
#> Attaching package: 'Biostrings'
#> 
#> 
#> The following object is masked from 'package:base':
#> 
#>     strsplit
#> 
#> 
#> Loading required package: Rsamtools
#> 
#> Loading required package: GenomicRanges
#> 
#> Loading required package: GenomicAlignments
#> 
#> Loading required package: SummarizedExperiment
#> 
#> Loading required package: MatrixGenerics
#> 
#> Loading required package: matrixStats
#> 
#> 
#> Attaching package: 'MatrixGenerics'
#> 
#> 
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> 
#> 
#> Loading required package: Biobase
#> 
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> 
#> 
#> Attaching package: 'Biobase'
#> 
#> 
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> 
#> 
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> 
#> 
#> Loading required package: dplyr
#> 
#> 
#> Attaching package: 'dplyr'
#> 
#> 
#> The following object is masked from 'package:ShortRead':
#> 
#>     id
#> 
#> 
#> The following objects are masked from 'package:GenomicAlignments':
#> 
#>     first, last
#> 
#> 
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> 
#> 
#> The following object is masked from 'package:matrixStats':
#> 
#>     count
#> 
#> 
#> The following objects are masked from 'package:GenomicRanges':
#> 
#>     intersect, setdiff, union
#> 
#> 
#> The following objects are masked from 'package:Biostrings':
#> 
#>     collapse, intersect, setdiff, setequal, union
#> 
#> 
#> The following object is masked from 'package:GenomeInfoDb':
#> 
#>     intersect
#> 
#> 
#> The following object is masked from 'package:XVector':
#> 
#>     slice
#> 
#> 
#> The following objects are masked from 'package:IRanges':
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> 
#> 
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> 
#> 
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, union
#> 
#> 
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> 
#> 
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> 
#> 
#> Loading required package: purrr
#> 
#> 
#> Attaching package: 'purrr'
#> 
#> 
#> The following object is masked from 'package:ShortRead':
#> 
#>     compose
#> 
#> 
#> The following object is masked from 'package:GenomicRanges':
#> 
#>     reduce
#> 
#> 
#> The following object is masked from 'package:XVector':
#> 
#>     compact
#> 
#> 
#> The following object is masked from 'package:IRanges':
#> 
#>     reduce
#> 
#> 
#> Loading required package: furrr
#> 
#> Loading required package: future
#> 
#> Loading required package: tidyr
#> 
#> 
#> Attaching package: 'tidyr'
#> 
#> 
#> The following object is masked from 'package:S4Vectors':
#> 
#>     expand
#> 
#> 
#> Loading required package: stringr
#> 
#> Loading required package: metacoder
#> 
#> This is metacoder version 0.3.6 (stable)
#> 
#> 
#> Attaching package: 'metacoder'
#> 
#> 
#> The following object is masked from 'package:ShortRead':
#> 
#>     reverse
#> 
#> 
#> The following objects are masked from 'package:Biostrings':
#> 
#>     complement, reverse
#> 
#> 
#> The following object is masked from 'package:XVector':
#> 
#>     reverse
#> 
#> 
#> The following object is masked from 'package:IRanges':
#> 
#>     reverse
#> 
#> 
#> Loading required package: phyloseq
#> 
#> 
#> Attaching package: 'phyloseq'
#> 
#> 
#> The following object is masked from 'package:metacoder':
#> 
#>     filter_taxa
#> 
#> 
#> The following object is masked from 'package:SummarizedExperiment':
#> 
#>     distance
#> 
#> 
#> The following object is masked from 'package:Biobase':
#> 
#>     sampleNames
#> 
#> 
#> The following object is masked from 'package:GenomicRanges':
#> 
#>     distance
#> 
#> 
#> The following object is masked from 'package:IRanges':
#> 
#>     distance
#> 
#> 
#> Loading required package: readr
#> 
#> Loading required package: ggplot2
#> 
#> 
#> Attaching package: 'ggplot2'
#> 
#> 
#> The following object is masked from 'package:metacoder':
#> 
#>     map_data
#> 
#> 
#> Loading required package: gridExtra
#> 
#> 
#> Attaching package: 'gridExtra'
#> 
#> 
#> The following object is masked from 'package:dplyr':
#> 
#>     combine
#> 
#> 
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> 
#> 
#> The following object is masked from 'package:BiocGenerics':
#> 
#>     combine
#> 
#> 
#> Loading required package: sessioninfo
document()
#> ℹ Updating rps10package documentation
#> ℹ Loading rps10package

# install.packages("devtools")
# devtools::install_github("grunwaldlab/rps10_metabarcoding_tool")
```

## Specify directory paths

Before you start, put your read fastq reads into a directory where you
would like to do your analysis. Please note: **xx** files and
directories will be created.

Please ensure that you have enough storage to save intermediate files
and reformatted directories on your computer to proceed.

If you are using the Silva of Unite databases, an ordinary personal
computer may not have enough memory.  
***TODO-by default, intermediate files will be saved to a temporary
directory. Provide instruction to user on how to override this.***

Prepare and input metadata and primer files according to the templates
provided. Using the provided templates, paste sample names in the first
column. In the second column, user will specify the barcode they have
used. Currently, user can specify ‘its’ or ‘rps10’. These primer names
are case sensitive.

The metadata file should look like this:

***TODO add screenshot***

You will also need to make a sheet containing info on the primers used
and forward and reverse primer sequences.

***TODO add screenshot***

If barcodes are pooled, user will include each sample name twice, and
under primer_name, either its or rps10 barcode.

***TODO add screenshot of pooled barcode scenario***

Include any other metadata columns you’d like to use during later steps
of the analysis. Formatting of names does not matter, but ensure that
there are no blank left cells.

``` r
#Change as needed based on where test dataset is located
directory_path<-"~/rps10package/raw_data/rps10_ITS" ##choose a directory for all downstream steps
primer_path <-file.path(directory_path, "primer_info.csv") ##modify .csv name or keep this name
metadata_path <-file.path(directory_path,"metadata.csv") ##modify .csv name or keep this name. The sample_name in the metadata sheet needs to match the first part (before first underscore), of the zipped raw FASTQ files
#Specify the location of where cutadapt is located
cutadapt_path<-"/opt/homebrew/bin/cutadapt"
```

## Reorganize data tables for downstream steps and obtain primer counts across all sample reads

The sample names, primer sequences and other metadata are reorganized in
preparation for running Cutadapt to remove primers.

``` r
data_tables <-
  prepare_reads(
    directory_path,
    primer_path,
    metadata_path,
    maxN = 0,
    multithread = TRUE
  )
#> Rows: 2 Columns: 3
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (3): primer_name, forward, reverse
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 4 Columns: 4
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (4): sample_name, primer_name, well, organism
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 16 Columns: 7
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (3): primer_name, orientation, sequence
#> dbl (4): S1_R1, S1_R2, S2_R1, S2_R2
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> # A tibble: 16 × 7
#>    primer_name orientation sequence               S1_R1 S1_R2 S2_R1 S2_R2
#>    <chr>       <chr>       <chr>                  <dbl> <dbl> <dbl> <dbl>
#>  1 rps10       forward     GTTGGTTAGAGYARAAGACT   35634     0 34601     1
#>  2 its         forward     CTTGGTCATTTAGAGGAAGTAA 48291     1 33485     0
#>  3 rps10       f_compt     CAACCAATCTCRTYTTCTGA       0     0     0     0
#>  4 its         f_compt     GAACCAGTAAATCTCCTTCATT     0     0     0     0
#>  5 rps10       f_rev       TCAGAARAYGAGATTGGTTG       0     0     0     0
#>  6 its         f_rev       AATGAAGGAGATTTACTGGTTC     0     0     0     0
#>  7 rps10       f_rc        AGTCTTYTRCTCTAACCAAC       0     8     0     9
#>  8 its         f_rc        TTACTTCCTCTAAATGACCAAG     0  2290     0  5316
#>  9 rps10       reverse     ATRYYTAGAAAGAYTYGAACT      0 35033     1 34233
#> 10 its         reverse     GCTGCGTTCTTCATCGATGC       1 47471     0 33214
#> 11 rps10       r_compt     TAYRRATCTTTCTRARCTTGA      0     0     0     0
#> 12 its         r_compt     CGACGCAAGAAGTAGCTACG       0     0     0     0
#> 13 rps10       r_rev       TCAAGYTYAGAAAGATYYRTA      0     0     0     0
#> 14 its         r_rev       CGTAGCTACTTCTTGCGTCG       0     0     0     0
#> 15 rps10       r_rc        AGTTCRARTCTTTCTARRYAT      4     0     4     0
#> 16 its         r_rc        GCATCGATGAAGAACGCAGC    3225     0  6732     0
```

<img src="man/figures/README-create data tables-1.png" width="100%" />
