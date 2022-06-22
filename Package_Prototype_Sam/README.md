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


* Cutadapt, the Python trimming program that is used. The version used here is 3.7. Go [here](https://cutadapt.readthedocs.io/en/stable/) to see documentation and installation instructions.


* Correct .csv files as inputs: In order to be successful, the user needs to provide both a primer data and metadata .csv file. These files need to make the exact formatting of the given metadata and primer info files provided on this project's GitHub. The metadata may have as many columns as needed as long as the first column is the sample names and the second is the primer names used in that sample. 

The primer data .csv needs to have the exact 3 columns shown, with as many rows as needed. 

<center>

![Primer Information CSV Example - names, forward, and reverse primers - all in their own columns. Column names need to match.](https://github.com/grunwaldlab/rps10_metabarcoding_tool/blob/main/screen_shots/primer_example.PNG) 
  
</center>


### How to Get the Development Environment Set Up and Running

## How to Use the Project

This will need to contain examples, screenshots, possible issues and how to fix. 

## Known Issues
* read_fastq function will not run correctly if the fastq data is not formatted correctly. This can be changed from inside of the function. The correct naming for the direction is 1 or R1.
* For Unix and Unix-like systems (macOS, Linux, etc.), the cutadapt_run function will catch the error if the cutadapt parameter of the function is in quotation marks. This is not true for Windows-based systems. 
* The cutadapt program will not execute correctly if the columns of the metadata CSV file is not in the correct order. The correct order should be: sample (or sample_id), primer, well, organism, etc.
* Before running cutadapt, if the cutadapt_data tibble has any N in it, the program will not run correctly.

## Credits 
Martha Sudermann, Samantha Dawson, Zachary Foster, Hung Phan


