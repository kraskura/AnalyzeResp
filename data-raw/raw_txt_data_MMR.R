## code to prepare `DATASET` dataset goes here

raw_txt_data_MMR<-read.delim("./test_file_TESTER_MMR.txt", skip=19)

# system.file("extdata", "raw_txt_data_MMR.txt", package = "respMR")

usethis::use_data(raw_txt_data_MMR, overwrite = TRUE, compress = "xz", internal=TRUE)
