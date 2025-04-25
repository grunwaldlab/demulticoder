## This is a resubmission. In this version I have:

* Removed redundant words from title and added single quotes to package and software names

* Clarified the meaning of acronyms in the description text

* Replaced \dontrun{} with \donttest{}. The examples can be executed, but they are not executable in < 5 seconds. Read files and reference databases were reduced in size, but the time still exceeds 5 seconds for some steps. To get meaningful outputs at the end, it would be difficult to reduce file sizes further.

* Replaced print()/cat() with message()/warning where these were present in functions.

* No functions write by default to the user's filespace. This was previously the case in example vignettes, by accident. Any files generated now are written to tempdir(). 

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
