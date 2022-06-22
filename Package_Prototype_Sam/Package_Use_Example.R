# you can now put the ?function to see what this function does!

?prepare_metadata



#run this for package stuff 
seed <- 1
set.seed(seed)

primer_path <- "absolute path to your primer csv"
metadata_path <- "absolute path to your metadata csv"
fastq_path <- "absolute path to your fastq files"
intermediate_path <- "absolute path to your intermediate folder"

create_intermediate()
primer_data <- prepare_primers(primer_path)
metadata <- prepare_metadata(metadata_path, primer_data)
#this may take while
fastq_data <- prepare_fastq(fastq_path, intermediate_path)
#this also could take a hot sec
pre_primer_hit_data <- pre_primer_hit_data(primer_data, fastq_data, intermediate_path)
print(pre_primer_hit_data)
pre_primer_plot <- primer_hit_plot(pre_primer_hit_data, fastq_data, metadata)
print(pre_primer_plot)
#editing ZACHS data with the library (names were mismatched)
#not a robust funtion - need it for testing (hopefully user will not have this issue)

#only use following line if using Zach's Mock Community
#fastq_data <- change_sample_ids(fastq_data)
cutadapt_data <- cutadapt_tibble(fastq_data, metadata)
#running the actual cutadapt program
cutadapt_run("your absolute path to cutadapt program", cutadapt_data)
#checking if cutadapt works (this may take awhile)
post_primer_hit_data <- post_primer_hit_data(pre_primer_hit_data, primer_data, new_fastq_data, intermediate_path)
print(post_primer_hit_data)
post_primer_plot <- primer_hit_plot(post_primer_hit_data, new_fastq_data, metadata)
print(post_primer_plot)

#filter steps
filter_and_trim(intermediate_path, new_fastq_data)




