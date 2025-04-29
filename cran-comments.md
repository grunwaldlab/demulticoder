## This is a second resubission. In this version, we:

* Checked the spellings in the DESCRIPTION file. Taxmap and phyloseq are spelled correctly but refer to R packages, so single quotes were added to these names. 

* Title was shortened.

* 

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
