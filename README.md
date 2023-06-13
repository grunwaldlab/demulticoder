
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rps10package

<!-- badges: start -->
<!-- badges: end -->

The goal of rps10package is to …

## Installation

You can install the development version of rps10package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("grunwaldlab/rps10_metabarcoding_tool")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(rps10package)
#> Loading required package: dada2
#> Loading required package: Rcpp
#> Loading required package: devtools
#> Loading required package: usethis
#> Loading required package: ShortRead
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
#>     table, tapply, union, unique, unsplit, which.max, which.min
#> Loading required package: BiocParallel
#> Loading required package: Biostrings
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: XVector
#> Loading required package: GenomeInfoDb
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
#> Loading required package: Rsamtools
#> Loading required package: GenomicRanges
#> Loading required package: GenomicAlignments
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
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
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:ShortRead':
#> 
#>     id
#> The following objects are masked from 'package:GenomicAlignments':
#> 
#>     first, last
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> The following object is masked from 'package:matrixStats':
#> 
#>     count
#> The following objects are masked from 'package:GenomicRanges':
#> 
#>     intersect, setdiff, union
#> The following objects are masked from 'package:Biostrings':
#> 
#>     collapse, intersect, setdiff, setequal, union
#> The following object is masked from 'package:GenomeInfoDb':
#> 
#>     intersect
#> The following object is masked from 'package:XVector':
#> 
#>     slice
#> The following objects are masked from 'package:IRanges':
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, union
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: purrr
#> 
#> Attaching package: 'purrr'
#> The following object is masked from 'package:ShortRead':
#> 
#>     compose
#> The following object is masked from 'package:GenomicRanges':
#> 
#>     reduce
#> The following object is masked from 'package:XVector':
#> 
#>     compact
#> The following object is masked from 'package:IRanges':
#> 
#>     reduce
#> Loading required package: furrr
#> Loading required package: future
#> Loading required package: tidyr
#> 
#> Attaching package: 'tidyr'
#> The following object is masked from 'package:S4Vectors':
#> 
#>     expand
#> Loading required package: stringr
#> Loading required package: metacoder
#> This is metacoder version 0.3.6 (stable)
#> 
#> Attaching package: 'metacoder'
#> The following object is masked from 'package:ShortRead':
#> 
#>     reverse
#> The following objects are masked from 'package:Biostrings':
#> 
#>     complement, reverse
#> The following object is masked from 'package:XVector':
#> 
#>     reverse
#> The following object is masked from 'package:IRanges':
#> 
#>     reverse
#> Loading required package: phyloseq
#> 
#> Attaching package: 'phyloseq'
#> The following object is masked from 'package:metacoder':
#> 
#>     filter_taxa
#> The following object is masked from 'package:SummarizedExperiment':
#> 
#>     distance
#> The following object is masked from 'package:Biobase':
#> 
#>     sampleNames
#> The following object is masked from 'package:GenomicRanges':
#> 
#>     distance
#> The following object is masked from 'package:IRanges':
#> 
#>     distance
#> Loading required package: readr
#> Loading required package: ggplot2
#> 
#> Attaching package: 'ggplot2'
#> The following object is masked from 'package:metacoder':
#> 
#>     map_data
#> Loading required package: gridExtra
#> 
#> Attaching package: 'gridExtra'
#> The following object is masked from 'package:dplyr':
#> 
#>     combine
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> The following object is masked from 'package:BiocGenerics':
#> 
#>     combine
#> Loading required package: sessioninfo
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
