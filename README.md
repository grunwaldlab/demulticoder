
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rps10package-working on name

## Description

<!-- badges: start -->
<!-- badges: end -->

Package readme in progress

## Installation

You can install the development version of rps10package from
[GitHub](https://github.com/) with:

*TODO-still working on documentation*

``` r
#temporary installation method
#ADJUST PATHS FOR MOMENT
devtools::load_all("/Users/masudermann/rps10package")
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

## Documentation

### Specify directory paths

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
directory_path<-"/Users/masudermann/rps10package/raw_data/rps10_ITS" ##choose a directory for all downstream steps
primer_path <-file.path(directory_path, "primer_info.csv") ##modify .csv name or keep this name
metadata_path <-file.path(directory_path,"metadata.csv") ##modify .csv name or keep this name. The sample_name in the metadata sheet needs to match the first part (before first underscore), of the zipped raw FASTQ files
#Specify the location of where cutadapt is located
cutadapt_path<-"/opt/homebrew/bin/cutadapt"
```

### Reorganize data tables for downstream steps and obtain primer counts across all sample reads

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
#>  1 rps10       forward     GTTGGTTAGAGYARAAGACT    4001     0  4800     0
#>  2 its         forward     CTTGGTCATTTAGAGGAAGTAA  5439     0  4573     0
#>  3 rps10       f_compt     CAACCAATCTCRTYTTCTGA       0     0     0     0
#>  4 its         f_compt     GAACCAGTAAATCTCCTTCATT     0     0     0     0
#>  5 rps10       f_rev       TCAGAARAYGAGATTGGTTG       0     0     0     0
#>  6 its         f_rev       AATGAAGGAGATTTACTGGTTC     0     0     0     0
#>  7 rps10       f_rc        AGTCTTYTRCTCTAACCAAC       0     0     0     1
#>  8 its         f_rc        TTACTTCCTCTAAATGACCAAG     0   295     0   721
#>  9 rps10       reverse     ATRYYTAGAAAGAYTYGAACT      0  3924     0  4751
#> 10 its         reverse     GCTGCGTTCTTCATCGATGC       0  5318     0  4501
#> 11 rps10       r_compt     TAYRRATCTTTCTRARCTTGA      0     0     0     0
#> 12 its         r_compt     CGACGCAAGAAGTAGCTACG       0     0     0     0
#> 13 rps10       r_rev       TCAAGYTYAGAAAGATYYRTA      0     0     0     0
#> 14 its         r_rev       CGTAGCTACTTCTTGCGTCG       0     0     0     0
#> 15 rps10       r_rc        AGTTCRARTCTTTCTARRYAT      0     0     0     0
#> 16 its         r_rc        GCATCGATGAAGAACGCAGC     392     0   901     0
```

<img src="man/figures/README-create data tables-1.png" width="100%" />

### Remove primers with Cutadapt and check that no primers still remain

If primers still remain, you may want to adjust the parameters or
manually remove any reads that still have primers (TODO-adjust this
instruction later).

***TODO-include a function to remove any ASV sequences that may still
have primers. There likely won’t be many-but it might be worth
incorporating a function like this***

``` r
cut_trim(
  directory_path,
  cutadapt_path,
  verbose = TRUE,
  maxEE = 2,
  truncQ = 5,
  minLen = 200,
  maxLen = 297,
  minCutadaptlength = 50
) 
#> Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
#> of ggplot2 3.3.4.
#> ℹ The deprecated feature was likely used in the dada2 package.
#>   Please report the issue at <https://github.com/benjjneb/dada2/issues>.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
#> Rows: 16 Columns: 11
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (3): primer_name, orientation, sequence
#> dbl (8): S1_R1_rps10, S1_R2_rps10, S1_R1_its, S1_R2_its, S2_R1_rps10, S2_R2_...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

<img src="man/figures/README-Trim function-1.png" width="100%" />

### Generate ASV abundance matrix.

Depending on chosen specifications, the minOverlap and maxMismatch can
be changed. The default values are the DADA2 default values

***To do-There is a message if analysis has already been run and parts
don’t re-run. Change so there is a way to override and rewrite files if
the user chooses***

``` r
asv_abund_matrix <-
  make_asv_abund_matrix(
    directory_path,
    minOverlap = 15,
    maxMismatch = 2,
    verbose = TRUE,
    multithread = TRUE
  )
#> [1] "File already exists"
#> Duplicate sequences detected and merged.
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

### Taxonomic identification steps

After databases are formatted, taxonomic information will be assigned to
the ASVs that were inferred previously. The output will be a csv file
containing each unique ASV, the abundance of each ASV across each
sample, and the taxonomic assignments, the bootstrap support values, as
well as a column listing the percent identity of each sequence to the
corresponding sequence in the reference database.

The user must specify which barcode they are using. Either ‘rps10’, its,
or ‘rps10_its’

**Databases** Todo-clarify instructions

Within the main directory, user must also have retrieved the rps10
oomycetedb (add link) or the Unite fungal database of choice. The names
of the databases must be ‘oomycetedb.fasta’ and ‘fungidb.fasta’, or else
they won’t be recognized. The headers will altered and new versions of
the databse will be saved. This will ensure there is consistency across
databases.

***TODO-save renamed db in temporary db or else they will take up too
much space.***

``` r
summary <- assignTax(
  directory_path,
  data_tables,
  asv_abund_matrix,
  multithread = TRUE,
  barcode = "rps10_its",
)
#> # A tibble: 55 × 7
#>    sequence                  dada2_tax dada2_pid S1_rps10 S2_rps10 S1_its S2_its
#>    <chr>                     <chr>     <chr>     <chr>    <chr>    <chr>  <chr> 
#>  1 GAAAATCTTTGTGTCGGTGGTTCA… Eukaryot… 72.91242… 1434     30       0      0     
#>  2 GAAAATCTTTGTGTCGGTGGTTCA… Eukaryot… 99.77168… 835      337      0      0     
#>  3 GAAAATCTTTGTGTCGGTGGTTCA… Eukaryot… 100       0        831      0      0     
#>  4 GAAAATCTTTGTGTCGGTGGTTCA… Eukaryot… 86.87782… 0        546      0      0     
#>  5 AAAAAGTCGTAACAAGGTTTCCGT… Eukaryot… 90.94488… 0        0        472    0     
#>  6 AAAAAGTCGTAACAAGGTTTCCGT… Eukaryot… 90.94488… 0        0        419    0     
#>  7 AAAAAGTCGTAACAAGGTTTCCGT… Eukaryot… 90.55118… 0        0        0      331   
#>  8 AAGTCGTAACAAGGTTTCCGTAGG… Eukaryot… 76.21145… 0        0        41     205   
#>  9 GAAAATCTTTGTGTCGGTGGTTCA… Eukaryot… 99.54441… 0        240      0      0     
#> 10 AAAAAGTCGTAACAAGGTTTCCGT… Eukaryot… 90.55118… 0        0        0      233   
#> # ℹ 45 more rows
#> Duplicate sequences detected and merged.
#>   sample_nameBarcode input filtered denoisedF denoisedR merged nonchim
#> 1           S1_rps10  4119     2284      2284      2274   2269    2269
#> 2           S2_rps10  5678     1424      1985      1985   1984    1984
#> 3             S1_its  4987     1987      1386      1379   1312    1312
#> 4             S2_its  4827     1394      1346      1341   1294    1294
```

### Convert ASV matrix to Taxmap and Phyloseq objects

While some users will appreciate having a matrix with all ASV sequences,
taxonomic assignments and abundances of each ASV per sample, others may
prefer to reformat matrices into Phyloseq or Taxmap objects.

Two wrapper functions are shared. They are derived from Metacoder
functions ‘parse_tax_data’ and ‘as_phyloseq’.

``` r
obj_dada<-asvmatrix_to_taxmap(asv_abund_matrix, min_read_depth=10, minimum_bootstrap=75)
#> Rows: 55 Columns: 7
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (2): sequence, dada2_tax
#> dbl (5): dada2_pid, S1_rps10, S2_rps10, S1_its, S2_its
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
phylo_obj<-taxmap_to_phyloseq(obj_dada)
```

***To do-paste in other examples of outputs***

## Dependencies

## Citation

## Future development
