# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
require(GBMutils)
# Set working directory
set_wd()


library(SeuratData)
InstallData("pbmc3k")

data("pbmc3k")
