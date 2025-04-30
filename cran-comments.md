## This is a third resubmission. In this version, we:  
* Account for the fact that cutadapt is not always installed an accessible when doing testthat tests. 

## This is a second resubmission. In this version, we:

* Checked the spellings in the DESCRIPTION file. The names in question ('taxmap' and 'phyloseq') are spelled correctly but refer to R packages, so single quotes were added to these names. 

* Title was shortened.

* Two references describing the methods were added to the description field of DESCRIPTION.

* Based on a previous suggestion, examples were wrapped with \donttest{} rather than \dontrun{}. Each of the functions can be run in less than 5 seconds, but the outputs of one are input into the next function, so to run the final functions, all previous steps must also be run, which increases the run time to over 20 seconds, even with reduced data set files.

* Since we cannot wrap examples, we now include tests using the testthat package. We checked that tests pass and now have tests in tests/testthat

* We double checked that functions don't write by default to users file space or package directory. We revised the setup_directories() function so that output directory is tempdir(). This was accidentally missed in the previous review.

* Additional revisions: made a few minor updates to the website to improve readability, highlight dependencies, etc. 

## This is a resubmission. In this version I have:

* Removed redundant words from title and added single quotes to package and software names

* Clarified the meaning of acronyms in the description text

* Replaced \dontrun{} with \donttest{}. The examples can be executed, but they are not executable in < 5 seconds. Read files and reference databases were reduced in size, but the time still exceeds 5 seconds for some steps. To get meaningful outputs at the end, it would be difficult to reduce file sizes further.

* Replaced print()/cat() with message()/warning where these were present in functions.

* No functions write by default to the user's filespace. This was previously the case in example vignettes, by accident. Any files generated now are written to tempdir(). 

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
