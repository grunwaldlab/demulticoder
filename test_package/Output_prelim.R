metadata <- returnList$metadata
print(metadata)
abundance <- read_csv(file.path("data/final_asv_abundance_table.csv"))
print(abundance)

abundance$species_boot = as.numeric(str_match(abundance$dada2_tax, pattern = "--([0-9]+)--Species")[,2])
abundance %>%
  filter(species_boot > 80, dada2_pid > 98) %>%
  mutate(dada2_tax = str_match(dada2_tax, pattern = "^(.+)--Species")[,2]) %>%
  mutate(dada2_pid = round(dada2_pid, 2)) %>%
  select(-species_boot) %>%
  select(-sequence, everything(), sequence) %>%
  write_csv(file = "~/species_found.csv")

abundance %>%
  mutate(dada2_tax = str_match(dada2_tax, pattern = "^(.+)--Species")[,2]) %>%
  mutate(dada2_pid = round(dada2_pid, 2), nr_blast_pid = round(nr_blast_pid, 2)) %>%
  select(-is_low_abund, -species_boot, -ref_blast_pid, -ref_blast_tax) %>%
  select(-sequence, everything(), sequence) %>%
  write_csv(file = "results/asvs_found.csv")