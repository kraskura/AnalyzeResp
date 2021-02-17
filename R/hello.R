

#' Title
#' @param x the name
#' @return from \code{\link{print}}
#' @export
#'
#' @examples
#' hello("Bradley")

#' \dontrun{
#' hello("Krista")
#' }
#'
hello <- function(x) {
  print(paste0("Hello, world!", x, "this is the function!"))
}
