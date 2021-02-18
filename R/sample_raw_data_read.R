

#' Read text file and then start the data analysis sequence
#'
#' @param path specifies the path to the internally available data
#'
#' @return a \code{tibble}
#' @export
#' @importFrom utils read.delim
#' @examples
#' txt = system.file("extdata", "raw_txt_data_MMR.txt", package = "AnalyzeResp")
#' sample_raw_data_read(txt)
sample_raw_data_read<-function(path){
  utils::read.delim(path, skip=19)
}

# these data will be available for tests, vignettes, etc.
