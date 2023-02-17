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
#' @param standardise_input whether to standardise the column names/types in the
#'   input dataframe. The default should work for flywire \code{fafbseg} and
#'   \code{neuprintr} input and ensure that we identify appropriate
#'   \code{pre_id}/\code{post_id} columns.
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
                                            inputcol="pre_id", outputcol="post_id",
                                            inputids=NULL, outputids=NULL,
                                            standardise_input=T) {
  if(standardise_input)
    x=standardise_partner_summary(x)
  cx=colnames(x)
  checkmate::assert_choice(inputcol, cx)
  checkmate::assert_choice(outputcol, cx)
  nas=is.na(x[[inputcol]]) | is.na(x[[outputcol]])
  if(any(nas)) {
    warning("Dropping: ", sum(nas), "/", length(nas), " rows due to missing input/output ids!")
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


standardise_partner_summary <- function(x) {
  # flywire_partner_summary2 gives columns like pre_pt_root_id
  colnames(x)=sub("_pt_root_", "_", colnames(x))
  cx=colnames(x)
  if("bodyid" %in% cx) {
    # looks like a neuprint result
    # but then doesn't have the expected fields ...
    stopifnot(all(c("partner", "prepost") %in% cx))
    # prepost=0 implies input query
    x$pre_id=ifelse(x$prepost==0, x$partner, x$bodyid)
    x$post_id=ifelse(x$prepost==1, x$partner, x$bodyid)
  } else {
    colnames(x)[cx=='query']=ifelse("pre_id" %in% cx, "post_id", "pre_id")
  }
  x$pre_id=bit64::as.integer64(x$pre_id)
  x$post_id=bit64::as.integer64(x$post_id)
  colnames(x)[cx=='cell_type']='type'
  x
}
