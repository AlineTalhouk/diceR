# Transpose while keeping sample ID and gene names in dimnames
library(magrittr)
raw <- get(load(file = "data-raw/hgsc-raw.rda"))
hgsc <- raw %>%
  set_rownames(.$UNIQID) %>%
  magrittr::extract(-1) %>%
  t() %>%
  as.data.frame() %>%
  structure(class.true = gsub(".*_", "\\1", rownames(.)))
devtools::use_data(hgsc, overwrite = TRUE, compress = "xz")
