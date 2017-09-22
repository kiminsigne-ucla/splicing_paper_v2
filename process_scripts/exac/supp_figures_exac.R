# ExAC control statistics
data_other %>% filter(sub_id == "SKP", v2_index <= 0.5) %>% nrow()
data_other %>% filter(sub_id == "SKP") %>% nrow()

data_other %>% filter(ensembl_id == "RANDOM-EXON", v2_index <= 0.5) %>% nrow()
data_other %>% filter(ensembl_id == "RANDOM-EXON") %>% nrow()

data_all %>% filter(sub_id == "BRK", v2_index <= 0.5) %>% nrow() 
data_all %>% filter(sub_id == "BRK") %>% nrow() 