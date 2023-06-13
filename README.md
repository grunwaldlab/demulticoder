
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
document()
#> Warning: [AssignTaxonomy.R:7] @param requires name and description
#> Warning: [AssignTaxonomy.R:8] @param requires name and description
#> Warning: [AssignTaxonomy.R:196] @param requires name and description
#> Warning: [AssignTaxonomy.R:197] @param requires name and description
#> Warning: [AssignTaxonomy.R:208] @param requires name and description
#> Warning: [AssignTaxonomy.R:209] @param requires name and description
#> Warning: [AssignTaxonomy.R:211] @param requires name and description
#> Warning: [AssignTaxonomy.R:229] @param requires name and description
#> Warning: [AssignTaxonomy.R:230] @param requires name and description
#> Warning: [AssignTaxonomy.R:231] @param requires name and description
#> Warning: [AssignTaxonomy.R:250] @param requires name and description
#> Warning: [AssignTaxonomy.R:275] @param requires name and description
#> Warning: [AssignTaxonomy.R:276] @param requires name and description
#> Warning: [AssignTaxonomy.R:285] @param requires name and description
#> Warning: [AssignTaxonomy.R:286] @param requires name and description
#> Warning: [AssignTaxonomy.R:296] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:102] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:103] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:122] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:123] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:199] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:257] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:306] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:307] @return requires a value
#> Warning: [MakeASVAbundMatrix.R:366] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:377] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:378] @param requires name and description
#> Warning: [MakeASVAbundMatrix.R:398] @param requires name and description

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
