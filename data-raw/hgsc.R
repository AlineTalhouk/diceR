# Transpose while keeping sample ID and gene names in dimnames
raw <- get(load(file = "data-raw/hgsc-raw.rda"))
hgsc <- raw %>% 
  magrittr::set_rownames(.$UNIQID) %>% 
  magrittr::extract(-1) %>% 
  t() %>% 
  as.data.frame()
devtools::use_data(hgsc, compress = "xz")
