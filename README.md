
# Demulticoder R package

**This package is actively under development. Additional testing,
documentation, website, and examples forthcoming. Until this message has
been removed, use with caution.**

### Introduction

The ***demulticoder*** package is a Cutadapt and DADA2 wrapper package
for metabarcodng analyses. The main commands and outputs are intuitive
and comprehensive, which helps to account for the complex and iterative
nature of metabarcoding analyses.

### Key features

- The ability to do analyses on either demultiplexed or pooled amplicons
  within sample  
- Amplicons from multiple datasets are trimmed of primers, filtered,
  denoised, merged, and given taxonomic assignments in one go (with
  different parameters for each dataset if desired)  
- The package handles not just 16S or ITS amplicons when using default
  UNITE fungal or Silva 16S databases for taxonomic assignment steps, but also oomycete rps10 analyses
  using oomycetedb (<https://oomycetedb.org>), or up to two custom
  databases (provided they are formatted as described here:
  <https://benjjneb.github.io/dada2/training.html>).

### Installation

To install the development version of package:

``` r
devtools::install_github("grunwaldlab/demulticoder")
```

### Check out the website to see how to use the package

For more information, key functions, inputs, and example vignettes,
check out the documentation at:
<https://grunwaldlab.github.io/demulticoder>

### Citation

Stay tuned

The package was developed by Martha Sudermann, Zachary Foster, Samantha
Dawson, Hung Phan, Niklaus Grunwald, Jeff Chang.
