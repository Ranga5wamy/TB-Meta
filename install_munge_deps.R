
# install_munge_deps.R

# install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# core packages you need
cran_pkgs <- c("data.table", "dplyr", "ggplot2", "tibble", "purrr")

to_install_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install_cran) > 0) {
  install.packages(to_install_cran, repos = "https://cloud.r-project.org")
}

# Bioconductor packages
bioc_pkgs <- c("MungeSumstats", "tidyGWAS")

to_install_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) {
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
}

message("All required packages are installed.")


