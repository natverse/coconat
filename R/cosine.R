#' Efficient calculation of cosine similarity for sparse or dense matrices
#'
#' @details This is much faster than e.g. \code{lsa::cosine} with the added
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


#' Cosine matrix utility functions
#'
#' @description These functions are intended for use by package authors rather
#'   than end users.
#'
#' @param x A matrix or a named list of input/output matrices
#' @param partners Whether to select input or output matrices when both are
#'   available
#' @param action Whether to zero out or drop any NA values in the cosine matrix
#'   (these may be present when some columns have no entries)
#'
#' @return A matrix. When both inputs and outputs are used these will be
#'   weighted by the total number of input and output synapses.
#' @export
prepare_cosine_matrix <- function(x, partners=c("inputs", "outputs"), action=c("zero", 'drop')) {
  x <- fix_nas(x, action=action)
  if(is.list(x)) {
    x <- if(length(partners)==2) {
      stopifnot(isTRUE(all.equal(dimnames(x$cin), dimnames(x$cout))))
      (x$cin*x$win+x$cout*x$wout)/(x$win+x$wout)
    }
    else if(partners=='outputs')
      x$cout
    else
      x$cin
  }
  x
}

fix_nas <- function(x, action=c("zero", 'drop')) {
  action=match.arg(action)
  if(is.list(x)) {
    # calculate the fixed matrices
    cico=sapply(x[c("cin", "cout")], fix_nas, action=action, simplify = F)
    # replace the original ones (keeping the win/wout elements)
    x[c("cin", "cout")]=cico
    return(x)
  }
  isnax=is.na(x)
  if(!any(isnax))
    return(x)
  if(action=='zero') {
    x[isnax]=0
    x
  } else {
    rs=rowSums(is.na(x))
    cs=colSums(is.na(x))
    # keep rows that have < ncol(x) NAs
    x[rs!=ncol(x), cs!=nrow(x)]
  }
}
