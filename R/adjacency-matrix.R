#' Convert a partner summary table into an adjacency matrix
#'
#' @details The \code{inputcol} and \code{outputcol} arguments can name columns
#'   containing other values besides the unique numerical identifiers for
#'   neurons. For example you can refer to a cell type column, thereby
#'   generating a \emph{grouped} connectivity matrix. This is very useful for
#'   bringing together neurons with similar connectivity patterns across brain
#'   hemispheres and individuals.
#'
#'   Passing a function to \code{inputids} and/or \code{outputids} allows
#'   partner neurons to be grouped with maximum flexibility. The input to the
#'   function will be the dataframe \code{x} (after standardisation if this has
#'   been requested). The output must be a single vector which can be
#'   interpreted as a factor to group partner neurons.
#'
#' @param x dataframe as produced by \code{flywire_partner_summary},
#'   \code{neuprint_connection_table} or equivalent.
#' @param sparse Whether to return a sparse matrix (default \code{TRUE} in case
#'   you are making a big one)
#' @param inputcol,outputcol Character vector specifying the columns containing
#'   input and output identifiers. See \bold{details} section for more
#'   information.
#' @param inputids,outputids Optional vectors of input/output ids to ensure that
#'   these are present in the output matrix. Alternatively these may contain a
#'   function that takes the dataframe \code{x} as input and returns a grouping
#'   vector. See \bold{details} section for more information.
#'
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
#' da2ds=neuprintr::neuprint_connection_table('DA2_lPN', details=TRUE, partners='out', threshold=5)
#' am=partner_summary2adjacency_matrix(da2ds, inputcol = 'bodyid', outputcol = 'partner')
#' library(Matrix)
#' image(am, ylab='DA2 PNs', xlab='outputs')
#'
#' amg=partner_summary2adjacency_matrix(da2ds, inputcol = 'bodyid', outputcol = 'type')
#' image(amg, ylab='DA2 PNs', xlab='output cell types')
#' }
partner_summary2adjacency_matrix<- function(x, sparse = TRUE,
                                            inputcol="pre_id", outputcol="post_id",
                                            inputids=NULL, outputids=NULL,
                                            standardise_input=T) {
  if(standardise_input)
    x=standardise_partner_summary(x)
  cx=colnames(x)
  if(is.function(inputids)) {
    inputcol='igroup'
    x[[inputcol]]=inputids(x)
    inputids=NULL
  } else checkmate::assert_choice(inputcol, cx)
  if(is.function(outputids)) {
    outputcol='ogroup'
    x[[outputcol]]=outputids(x)
    outputids=NULL
  } else checkmate::assert_choice(outputcol, cx)

  nas=is.na(x[[inputcol]]) | is.na(x[[outputcol]])
  if(any(nas)) {
    warning("Dropping: ",
            sum(nas), "/", length(nas), " neurons representing ",
            sum(x[['weight']][nas]), "/", sum(x[['weight']]),
    " synapses due to missing ids!")
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
