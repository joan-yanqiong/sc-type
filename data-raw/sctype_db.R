## code to prepare `sctype_db` dataset goes here
ScTypeDB_full <- openxlsx::read.xlsx("data-raw/ScTypeDB_full.xlsx")
ScTypeDB_short <- openxlsx::read.xlsx("data-raw/ScTypeDB_short.xlsx")

usethis::use_data(ScTypeDB_full, ScTypeDB_short, overwrite = TRUE)
