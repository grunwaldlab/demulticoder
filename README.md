
<img src="man/figures/demulticoder_logo_iteration1.png" align="right" height="95" alt="" />

## **demulticoder**: An R package for the simultaneous analysis of multiplexed metabarcodes

### Introduction

The **demulticoder** package is a **Cutadapt** and **dada2** wrapper
package for metabarcodng analyses. The main commands and outputs are
intuitive and comprehensive, which helps to account for the complex and
iterative nature of metabarcoding analyses.

Here is a brief schematic of the general workflow:

<img src="man/figures/Figure1.svg" width="35%" height="35%" style="display: block; margin: auto auto auto 0;" />

### Key Features

- It automates the use of **dada2** to analyze data derived from
  multiple metabarcodes.  
- It reduces the number of manual input steps  
- Handles analysis of two metabarcodes multiplexed into the same
  sequencing batch  
- Analyze different types of metabarcodes simultaneously  
- Reproducible workflows for oomycetes
- Supported metabarcodes: 16S rDNA, ITS1, *rps10*, and up to two
  additional metabarcodes

### Installation

**Dependencies**:  
First install **Cutadapt** program following the instructions here:
<https://cutadapt.readthedocs.io/en/stable/installation.html>

Let’s locate where the **Cutadapt** executable is. You must do this from
a **Terminal** window:

``` sh
#If you installed with pip or pipx, or homebrew, run this command from a Terminal window
which cutadapt
cutadapt --version
```

If you followed the **Cutadapt** installation instructions to create a
conda environment called cutadapt (change to whatever you named your
environment), to install it in, open up a **Terminal** window and type
these commands:

``` sh
#Run commands from a Terminal window
conda activate cutadapt
which cutadapt
cutadapt --version
```

Second, make sure the following R packages are installed:

- **dada2** (Latest version is 3.20)
  - To install, follow these instructions:
    <https://www.bioconductor.org/packages/release/bioc/html/dada2.html>  
- **phyloseq**
  - To install:
    <https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html>  
- **metacoder** (Available through CRAN)

To install the development version of package (while submission to CRAN
is in progress):

``` r
#Here we install demulticoder (instructions will be updated once available through CRAN)
devtools::install_github("grunwaldlab/demulticoder")
library("demulticoder")

#Let's make sure other packages are loaded:
library("devtools")
library("dada2")
library("phyloseq")
library("metacoder")
```

### Quick start

**1. Set-up input directory and files**

To demonstrate how to use the package, we have a small test data set
that comes loaded with the package. This data set will be used in the
workflow example below.

Already loaded in the test data set directory are the following files:

- **PE short read amplification data**
  - Files: S1_R1.fastq.gz, S1_R2.fastq.gz, S2_R1.fastq.gz,
    S2_R1.fastq.gz  
  - The files must end in either *R1.fastq.gz* , or *R2.fastq.gz* and
    each sample must have both R1 and R2 files.
- [**metadata.csv**](https://github.com/grunwaldlab/demulticoder/blob/main/inst/extdata/metadata.csv)
  - New row for each unique sample
  - Samples entered twice if samples contain two pooled metabolites, as
    in the test data template
- [**primerinfo_params.csv**](https://github.com/grunwaldlab/demulticoder/blob/main/inst/extdata/primerinfo_params.csv)
  - New row for each unique barcode and associated primer sequence
  - Optional **Cutadapt** and **dada2** parameters
- **Taxonomy databases**
  - UNITE fungal database (abridged version)
  - completed

See
[**Documentation**](https://grunwaldlab.github.io/demulticoder/articles/Documentation.html)
for how to format databases and input files.

For more details on each step, check out the [**Getting
Started**](https://grunwaldlab.github.io/demulticoder/articles/Getting_started.html)
tab on the package website

**2. Prepare reads**

``` r
output<-prepare_reads(
  data_directory = system.file("extdata", package = "demulticoder"), # This allows us to use the test directory located within the package
  output_directory = tempdir(), # Change to you preferred location on your local computer (Example: "~/demulticoder_test")
  overwrite_existing = TRUE)
```

**3. Cut and trim reads** User must install **Cutadapt** on their local
machine and append the path to the executable.

``` r
cut_trim(
  output,
  cutadapt_path="/usr/bin/cutadapt", # Change to the location on your computer. (Example: "/usr/bin/cutadapt")
  overwrite_existing = TRUE) 
```

**4. Make ASV abundance matrix**

``` r
make_asv_abund_matrix(
  output,
  overwrite_existing = TRUE)
```

**5. Assign taxonomy**

``` r
assign_tax(
  output,
  asv_abund_matrix,
  overwrite_existing = TRUE)
```

**6. Convert ASV matrix to taxmap and phyloseq objects**

``` r
objs<-convert_asv_matrix_to_objs(output)
```

### Check out the website to view the documentation and see more examples

For more information on source code, check out the package repository:
<https://grunwaldlab.github.io/demulticoder/>

### For source code:

<https://github.com/grunwaldlab/demulticoder/>

### Citation

The package was developed by Martha Sudermann, Zachary Foster, Samantha
Dawson, Hung Phan, Jeff Chang, and Niklaus Grünwald

An associated manuscript is currently under review.
