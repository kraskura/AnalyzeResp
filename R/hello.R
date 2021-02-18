

#' Title
#' @param x the name
#' @return from \code{\link{print}}
#' @export
#' @importFrom tibble as_data_frame
#' @examples
#' hello("Bradley")

#' \dontrun{
#' hello("Krista")
#' }
#'
hello <- function(x) {
  print(paste0("Hello, world!", x, "this is the function!"))
  as_data_frame(matrix(ncol=2, nrow=1))
}
