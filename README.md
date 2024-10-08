
# *Demulticoder*-An R package to streamline metabarcoding analyses using DADA2

<img src="man/figures/logo.png" align="right" height="139" alt="" />

*This package is actively under development. Until this message has been
removed, use with caution. Additional testing, documentation, and
examples are in progress.*

### Introduction

The ***demulticoder*** package is a Cutadapt and DADA2 wrapper package
for metabarcodng analyses. The main commands and outputs are intuitive
and comprehensive, which helps to account for the complex and iterative
nature of metabarcoding analyses.

Here is a brief schematic of the general workflow:

<img src="man/figures/rps10_fig1_smaller.drawio.png" width="60%" height="60%" style="display: block; margin: auto auto auto 0;" />

### Key features

- The ability to do analysis on either demultiplexed or pooled amplicons
  within samples

- Amplicons from multiple datasets be trimmed of primers, filtered,
  denoised, merged, and given taxonomic assignments in one go (with
  different parameters for each dataset if desired)

- The package handles not just 16S or ITS datasets when using default
  UNITE fungal or Silva 16S databases but also oomycete rps10 analyses
  using [oomycetedb](https://grunwaldlab.github.io/OomyceteDB//), or up to two custom
  databases (provided they are formatted as described here:
  <https://benjjneb.github.io/dada2/training.html>).

### Installation

To install the development version of package:

``` r
devtools::install_github("grunwaldlab/demulticoder")
```

### Quick start

**1. Set-up input directory and files**

After installing the package, make a data directory and add the
following files:  
- PE short read amplicon data. The files must end in either
\*\_R1.fastq.gz\* , or \*\_R2.fastq.gz\* and each sample must have both
R1 and R2 files.

- [**metadata.csv**](https://github.com/grunwaldlab/demulticoder/blob/main/inst/extdata/metadata.csv)
  file (there will be a unique row for each sample, and samples will be
  entered twice if they contain pooled amplicions, as in the example
  template)

- [**primerinfo_params.csv**](https://github.com/grunwaldlab/demulticoder/blob/main/inst/extdata/primerinfo_params.csv)
  file (there will be a new row for each unique barcode and associated
  primer sequences, and there are also optional Cutadapt, DADA2 or
  filtering parameters that can be added or adjusted)

**2. Prepare reads**

``` r
output<-prepare_reads(
  data_directory = "<DATADIR>",
  output_directory = "<OUTDIR>")
```

**3. Cut and trim reads**

``` r
cut_trim(
  output,
  cutadapt_path="<CUTADAPTPATH>")
```

**4. Make ASV abundance matrix**

``` r
make_asv_abund_matrix(
  output)
```

**5. Assign taxonomy**

``` r
assign_tax(
  output,
  asv_abund_matrix)
```

**6. Convert ASV matrix to taxmap and phyloseq objects**

``` r
objs<-convert_asv_matrix_to_objs(output)
```

### Check out the website to view the documentation and see more examples

For more information, key functions, inputs, and example vignettes,
check out the documentation at:
<https://grunwaldlab.github.io/demulticoder>

### Citation

The package was developed by Martha Sudermann, Zachary Foster, Samantha
Dawson, Hung Phan, Jeff Chang, and Niklaus Grunwald.

Stay tuned for associated manuscript.
