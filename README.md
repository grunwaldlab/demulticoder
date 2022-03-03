# rps10_metabarcoding_tool
***

## Project Description

This project is a pipeline to get from raw fastq data to trimmed and filter sequences. The input will be three items; primer information CSV file, a metadata information CSV file, and fastq files from sequencing. The output will be filtered sequences that are ready to move on to other analysis, such as quality analysis or taxonomic identification. Let's just change this an test it out. 

***

## General Info

This project is still in development.

## How to Install and Run the Project 

### Dependencies 

* R packages needed: dada2, ShortRead, BioStrings, dplyr, purrr, furrr, tidyr, readr, ggplot2, gridExtra, sessioninfor (TODO: put versions)

<<<<<<< HEAD
* Cutadapt, the Python trimming program that is used. The version used here is 3.4. Go [here](https://cutadapt.readthedocs.io/en/stable/) to see documentation and installation instructions.
=======
* Cutadapt, the Python trimming program that is used. The version used here is 3.5. Go [here](https://cutadapt.readthedocs.io/en/stable/) to see documentation and installation instructions.
>>>>>>> 3fdb8e114ead4c3e3c5275bf4c8ac094ef90a9de

* Correct .csv files as inputs: In order to be successful, the user needs to provide both a primer data and metadata .csv file. These files need to make the exact formatting of the given metadata and primer info files provided on this project's GitHub. The metadata may have as many columns as needed as long as the first column is the sample names and the second is the primer names used in that sample. 

The primer data .csv needs to have the exact 3 columns shown, with as many rows as needed. 

<center>

![Primer Information CSV Example - names, forward, and reverse primers - all in their own columns. Column names need to match.](https://github.com/grunwaldlab/rps10_metabarcoding_tool/blob/main/screen_shots/primer_example.PNG) 
  
</center>


### How to Get the Development Environment Set Up and Running

## How to Use the Project

This will need to contain examples, screenshots, possible issues and how to fix. 

## Credits 
Martha Sudermann, Samantha Dawson, Zachary Foster, Hung Phan


