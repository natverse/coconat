#' Efficient calculation of cosine similarity for sparse or dense matrices
#'
#' @details This is much faster than eg \code{lsa::cosine} with the added
#'   benefit that it works for sparse input matrices
#'
#' @param x A data matrix suitable for clustering
#' @param sparse Whether to return a sparse or dense matrix (default dense)
#' @param transpose When \code{F} (the default) calculates the cosine distance
#'   between columns. When \code{T} calculates the distance between rows.
#'
#' @return A square matrix, dense unless \code{sparse=TRUE} and \code{x} is
#'   sparse.
#' @export
#'
#' @examples
#' da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
#' am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
#' cosine_sim(am)
#'
#' \dontrun{
#' fam_pnkc2=flywire_adjacency_matrix2(
#'   inputids = 'class:ALPN_R', outputids = 'class:Kenyon_Cell_R',
#'   sparse=T , threshold = 2)
#' kckc.cos=cosine_sim(fam_pnkc2)
#' pnpn.cos=cosine_sim(fam_pnkc2, transpose=T)
#' }
cosine_sim <- function(x, sparse=FALSE, transpose=FALSE) {
  cx=class(x)
  if(!is.matrix(x) && !isTRUE(attr(cx, "package") == "Matrix"))
    stop("I don't recognise that as a matrix!")
  cpx <- if(transpose) Matrix::tcrossprod(x) else Matrix::crossprod(x)
  cosx=Matrix::cov2cor(cpx)
  if(sparse) cosx else as.matrix(cosx)
}
