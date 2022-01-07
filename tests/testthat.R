library(testthat)
library(loosends)
Sys.setenv("DEFAULT_GENOME" = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens", "R_TESTS"="")

test_check("loosends")
