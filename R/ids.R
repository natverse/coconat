#' Convert any id into a character vector
#'
#' @param x A character or numeric vector of ids
#'
#' @return A character vector
#' @export
#'
#' @examples
#' id2char(1000)
#' id2char("1000")
#' id2char(1e5)
id2char <- function(x) {
  if(is.character(x)) return(x)
  if(is.numeric(x)) {
    cx=as.character(bit64::as.integer64(x))
    return(cx)
  }
  # because class(NA)=="logical"
  if(is.logical(x) && all(is.na(x)))
    return(rep(NA_character_, length(x)))
  stop("I don't know how to handle input of class:",
       paste(class(x), collapse=','))
}
