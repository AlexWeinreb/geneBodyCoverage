# Because these packages are used in inst/genebody_coverage.R we have them in Imports,
# but we don't use them in the actual package functions so check complains
ignore_unused_imports <- function() {
  getopt::getopt
  dplyr::select
  tidyr::unnest
  readr::read_delim
}
