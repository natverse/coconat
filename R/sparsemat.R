#' Efficient column and row scaling for sparse matrices
#'
#' @details Conventionally in connectivity matrices, rows are upstream input
#'   neurons and columns are downstream output neurons. \code{colScaleM} will
#'   therefore normalise by the total input onto each downstream neuron while
#'   \code{rowScaleM} will normalise by the total output from each upstream
#'   neuron.
#'
#'   These functions will work with dense inputs but are less efficient than
#'   their base R cousins (e.g. \code{scale}). As an example, normalising a 4x4k
#'   subset of VNC connectivity data took \itemize{
#'
#'   \item \code{colScaleM} sparse 3.5ms
#'
#'   \item \code{colScaleM} dense 1460ms
#'
#'   \item \code{scale} dense 365ms
#'
#'   }
#'
#'   In conclusion sparse matrices with \code{colScaleM}/\code{rowScaleM} are
#'   the way to go for real connectivity matrices! For large matrices they will
#'   be \emph{much} faster and for small ones, both methods will be fast anyway.
#'
#'   Note also that these functions were originally called \code{colScale} etc
#'   but then the \code{\link[Matrix]{Matrix}} package added functions of the
#'   same name. The functions in this package differ in two respects from the
#'   \code{Matrix::\link[Matrix:colScale]{colScale,rowScale}} functions. First
#'   they calculate the required row or column sums to perform the scaling.
#'   Second, they handle the inevitable zero sum rows or columns via the default
#'   \code{na.rm} argument.
#'
#' @param A a sparse (or dense) matrix
#' @param na.rm when \code{T}, the default, converts any NA values (usually
#'   resulting from divide by zero errors when a neuron has no partners) to 0
#' @description \code{colScale} normalises a matrix by the sum of each column
#'
#' @return A scaled sparse matrix
#' @export
#'
#' @seealso \code{Matrix::\link[Matrix]{colScale}} on which these are based and,
#'   for the original version,
#'   \url{https://stackoverflow.com/questions/39284774/column-rescaling-for-a-very-large-sparse-matrix-in-r}
#'
#' @examples
#' library(Matrix)
#' set.seed(42)
#' A <- Matrix(rbinom(100,10,0.05), nrow = 10)
#' rownames(A)=letters[1:10]
#' colnames(A)=LETTERS[1:10]
#'
#' colScaleM(A)
#' rowScaleM(A)
#' geomScaleM(A)
colScaleM <- function(A, na.rm=TRUE) {
  scalefac = 1 / Matrix::colSums(A)
  if(na.rm) scalefac[!is.finite(scalefac)]=0
  Matrix::colScale(A, scalefac)
}

#' @export
#' @rdname colScaleM
#' @description \code{rowScale} normalises a matrix by the sum of each row
rowScaleM <- function(A, na.rm=TRUE) {
  scalefac = 1 / Matrix::rowSums(A)
  if(na.rm) scalefac[!is.finite(scalefac)]=0
  Matrix::rowScale(A, scalefac)
}

#' @export
#' @rdname colScaleM
#' @description sqrt of both column and row normalisation
geomScaleM <- function(A, na.rm=TRUE) {
  sf1=sqrt(1/Matrix::rowSums(A))
  if(na.rm) sf1[!is.finite(sf1)]=0
  sf2=sqrt(1/Matrix::colSums(A))
  if(na.rm) sf2[!is.finite(sf2)]=0
  Matrix::dimScale(A, d1 = sf1, d2 = sf2)
}

#' Sparse 0/1 aggregation matrix mapping ids to groups
#'
#' @description Construct a sparse matrix with one row per id and one column per
#'   group, containing a 1 where an id belongs to a group. Right-multiplying a
#'   connectivity matrix by \code{grouping_matrix(colnames(M), group)} sums the
#'   columns of \code{M} within each group (e.g. collapsing partner neurons to
#'   cell types).
#'
#' @param ids Character vector (or coercible) of identifiers, typically the row
#'   or column names of a matrix to be aggregated.
#' @param group Grouping labels. Either a vector parallel to \code{ids} or a
#'   \emph{named} vector that will be looked up by \code{ids}. Ids whose group
#'   is \code{NA} are dropped (contribute no 1s), so they vanish from any
#'   aggregation.
#'
#' @return A sparse \code{\link[Matrix]{Matrix}} of dimension
#'   \code{length(ids)} x \code{number of groups} with ids as row names and
#'   group labels as column names.
#' @export
#' @seealso \code{\link{effective_connectivity}}
#' @examples
#' library(Matrix)
#' M <- Matrix(c(1,0,2, 0,3,0), nrow=2, byrow=TRUE,
#'   dimnames=list(c('a','b'), c('x','y','z')))
#' # collapse columns x,y -> g1 and z -> g2
#' G <- grouping_matrix(colnames(M), c(x='g1', y='g1', z='g2'))
#' M %*% G
grouping_matrix <- function(ids, group) {
  ids <- as.character(ids)
  if(!is.null(names(group)))
    group <- group[ids]
  if(length(group)!=length(ids))
    stop("`group` must be parallel to `ids` or a named vector indexed by `ids`!")
  group <- as.character(group)
  glevels <- unique(group[!is.na(group)])
  j <- match(group, glevels)
  keep <- !is.na(j)
  Matrix::sparseMatrix(
    i = which(keep),
    j = j[keep],
    x = 1,
    dims = c(length(ids), length(glevels)),
    dimnames = list(ids, glevels))
}

# expand/reorder the rows of B so that rownames(B) == ids, filling missing
# ids with zero rows. Used to align consecutive matrices for multiplication.
align_rows <- function(B, ids) {
  ids <- as.character(ids)
  if(is.null(rownames(B)))
    stop("Cannot align a matrix without row names!")
  m <- match(ids, rownames(B))
  ok <- !is.na(m)
  sel <- Matrix::sparseMatrix(
    i = which(ok),
    j = m[ok],
    x = 1,
    dims = c(length(ids), nrow(B)),
    dimnames = list(ids, rownames(B)))
  sel %*% B
}

#' Effective connectivity through multi-step pathways
#'
#' @description Compute an estimate of the effective connectivity between the
#'   inputs of the first matrix and the outputs of the last matrix in a chain,
#'   passing through one or more intermediate neuron layers. This follows the
#'   standard approach of input-normalising each step (so the inputs to every
#'   postsynaptic cell sum to 1) and then multiplying the matrices together.
#'
#' @details Each supplied matrix must be oriented with presynaptic (input)
#'   neurons as rows and postsynaptic (output) neurons as columns (the
#'   convention of \code{\link{partner_summary2adjacency_matrix}}). Consecutive
#'   matrices are chained by name: the columns of matrix \code{k} are matched to
#'   the rows of matrix \code{k+1}, and any output of matrix \code{k} that is
#'   absent from the rows of matrix \code{k+1} contributes a zero (dead-end)
#'   path.
#'
#'   When \code{normalise=TRUE} (the default) every matrix is column-normalised
#'   with \code{\link{colScaleM}} \emph{before} multiplication, implementing the
#'   assumption that an interneuron conveys information about its inputs in
#'   proportion to their synaptic weight. This is the procedure described by
#'   Schlegel et al. (2021) \doi{10.7554/eLife.62576}.
#'
#'   \code{group} collapses the final (output) dimension to groups such as cell
#'   types. Because grouping is applied \emph{after} normalisation and
#'   multiplication, and matrix multiplication is associative, grouping the
#'   final product is identical to grouping the last matrix's columns before the
#'   earlier multiplications - the per-neuron normalisation is always preserved.
#'
#'   \bold{Iterative use}. Passing every hop as one list is convenient but
#'   materialises all the adjacency matrices at once, which can be a memory hog
#'   for deep or highly connected walks. The more scalable (and often more
#'   common) pattern is to advance one hop at a time, keeping only the running
#'   product and the current hop in memory: normalise each new hop yourself with
#'   \code{\link{colScaleM}} and call \code{effective_connectivity(list(running,
#'   step), normalise = FALSE)}. Here \code{normalise = FALSE} is essential -
#'   \code{running} is already a product of normalised matrices and must not be
#'   rescaled, while \code{step} you have normalised in advance. This also lets
#'   you inspect, threshold or prune the intermediate frontier between hops
#'   (e.g. dropping weakly connected partners before fetching the next layer, as
#'   \code{coconatfly}'s multihop clustering does) rather than committing to the
#'   whole chain up front.
#'
#' @param matrices A list of sparse (or dense) adjacency matrices to chain, each
#'   with presynaptic neurons as rows and postsynaptic neurons as columns and
#'   with row/column names to align successive matrices.
#' @param group Optional grouping for the columns (outputs) of the final matrix,
#'   passed to \code{\link{grouping_matrix}} - either parallel to the final
#'   columns or a named vector indexed by them.
#' @param normalise Whether to column-normalise (input-normalise) each matrix
#'   with \code{\link{colScaleM}} before multiplication (default \code{TRUE}).
#'
#' @return A sparse matrix of effective connectivity with the inputs of the
#'   first matrix as rows and the outputs (or output groups) of the last matrix
#'   as columns.
#' @export
#' @seealso \code{\link{colScaleM}}, \code{\link{grouping_matrix}},
#'   \code{\link{partner_summary2adjacency_matrix}}
#' @references Schlegel et al. (2021) \emph{eLife} \doi{10.7554/eLife.62576}
#' @examples
#' library(Matrix)
#' set.seed(42)
#' # query (Q) -> interneurons (I) -> targets (T)
#' A1 <- Matrix(rbinom(6,5,0.4), nrow=2,
#'   dimnames=list(c('q1','q2'), c('i1','i2','i3')))
#' A2 <- Matrix(rbinom(9,5,0.4), nrow=3,
#'   dimnames=list(c('i1','i2','i3'), c('t1','t2','t3')))
#' effective_connectivity(list(A1, A2))
#'
#' # equivalent iterative form: advance one hop at a time, normalising each new
#' # hop yourself and keeping only the running product (memory friendly, and you
#' # can prune the frontier between hops)
#' running <- colScaleM(A1)
#' running <- effective_connectivity(list(running, colScaleM(A2)),
#'                                   normalise = FALSE)
#' running
effective_connectivity <- function(matrices, group=NULL, normalise=TRUE) {
  if(!is.list(matrices) || length(matrices)<1)
    stop("`matrices` must be a list of one or more matrices!")
  if(isTRUE(normalise))
    matrices <- lapply(matrices, colScaleM)
  eff <- matrices[[1]]
  for(k in seq_along(matrices)[-1]) {
    B <- align_rows(matrices[[k]], colnames(eff))
    eff <- eff %*% B
  }
  if(!is.null(group))
    eff <- eff %*% grouping_matrix(colnames(eff), group)
  eff
}
