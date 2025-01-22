## code to prepare `example_data` dataset goes here
scType_scRNAseq_example_data <- readRDS("data-raw/exampleData.RDS")

scType_Visium_example_data <- readRDS("data-raw/frontal_cortex_subset.RDS")

usethis::use_data(scType_scRNAseq_example_data, scType_Visium_example_data, overwrite = TRUE)
