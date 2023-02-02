#' Convert a partner summary table into an adjacency matrix
#'
#' @param x dataframe as produced by \code{flywire_partner_summary},
#'   \code{neuprint_connection_table} or equivalent.
#' @param sparse Whether to return a sparse matrix (default \code{TRUE} in case
#'   you are making a big one)
#' @param inputcol,outputcol Character vector specifying the columns containing
#'   input and output identifiers.
#' @param inputids,outputids Optional vectors of input/output ids to ensure that
#'   these are present in the output matrix.
#'
#' @return A matrix with inputs as rows and outputs (downstream neurons) as
#'   columns
#' @export
#'
#' @examples
#' \dontrun{
#' da2ds=neuprintr::neuprint_connection_table('DA2', details=TRUE, partners='out', threshold=15)
#' am=partner_summary2adjacency_matrix(da2ds, inputcol = 'partner', outputcol = 'bodyid')
#' library(Matrix)
#' image(am, xlab='DA2 PNs', ylab='outputs')
#' }
partner_summary2adjacency_matrix<- function(x, sparse = TRUE,
                                             inputcol=NULL, outputcol=NULL,
                                             inputids=NULL, outputids=NULL) {
  # flywire_partner_summary2 gives columns like pre_pt_root_id
  cx=colnames(x)
  checkmate::assert_choice(inputcol, cx)
  checkmate::assert_choice(outputcol, cx)
  nas=is.na(x[[inputcol]]) | is.na(x[[outputcol]])
  if(any(nas)) {
    warning("Dropping: ", sum(nas), " rows with missing input/output ids!")
    x=x[!nas,,drop=F]
  }
  if(is.null(inputids))
    inputids=unique(x[[inputcol]])
  else {
    inputids=unique(inputids)
    stopifnot(all(unique(x[[inputcol]]) %in% inputids))
  }
  if(is.null(outputids))
    outputids=unique(x[[outputcol]])
  else {
    outputids=unique(outputids)
    stopifnot(all(unique(x[[outputcol]]) %in% outputids))
  }
  i=match(x[[inputcol]], inputids)
  j=match(x[[outputcol]], outputids)
  sm = Matrix::sparseMatrix(
    i = i,
    j = j,
    dims = c(length(inputids), length(outputids)),
    x = x$weight,
    dimnames = list(as.character(inputids), as.character(outputids))
  )
  if (isTRUE(sparse))
    sm
  else as.matrix(sm)
}
