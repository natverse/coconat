#' Generic connectivity similarity between columns (or rows) of a matrix
#'
#' @description \code{connectivity_similarity} computes pairwise similarity
#'   between columns (or rows) of a matrix using a pluggable metric. This is the
#'   main generic function for connectivity-based clustering; individual metric
#'   functions like \code{\link{cosine_sim}} and \code{\link{jaccard_sim}} can
#'   also be used directly.
#'
#' @param x A (sparse) matrix, typically an adjacency matrix from
#'   \code{\link{partner_summary2adjacency_matrix}}
#' @param metric Character specifying the similarity metric. One of
#'   \code{"cosine"}, \code{"jaccard"}, \code{"weighted_jaccard"}, or
#'   \code{"tanimoto"}.
#' @param sparse Whether to return a sparse matrix (default \code{FALSE})
#' @param transpose When \code{FALSE} (the default) calculates similarity
#'   between columns. When \code{TRUE} calculates similarity between rows.
#'
#' @return A square similarity matrix with values in \code{[0,1]}.
#' @export
#' @seealso \code{\link{cosine_sim}}, \code{\link{jaccard_sim}},
#'   \code{\link{tanimoto_sim}}
#' @examples
#' da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
#' am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
#' connectivity_similarity(am, metric="cosine")
#' connectivity_similarity(am, metric="jaccard")
#' connectivity_similarity(am, metric="weighted_jaccard")
#' connectivity_similarity(am, metric="tanimoto")
connectivity_similarity <- function(x, metric = c("cosine", "jaccard", "weighted_jaccard", "tanimoto"),
                                    sparse = FALSE, transpose = FALSE) {
  metric <- match.arg(metric)
  switch(metric,
    cosine = cosine_sim(x, sparse = sparse, transpose = transpose),
    jaccard = jaccard_sim(x, weighted = FALSE, sparse = sparse, transpose = transpose),
    weighted_jaccard = jaccard_sim(x, weighted = TRUE, sparse = sparse, transpose = transpose),
    tanimoto = tanimoto_sim(x, sparse = sparse, transpose = transpose)
  )
}


#' Jaccard similarity for sparse or dense matrices
#'
#' @description Computes pairwise Jaccard similarity between columns (or rows)
#'   of a matrix. When \code{weighted=FALSE}, uses binary Jaccard (presence/
#'   absence). When \code{weighted=TRUE}, uses the generalised (weighted) Jaccard
#'   index: \code{sum(min(a,b)) / sum(max(a,b))}.
#'
#' @details Both variants are optimised for sparse matrices. The binary variant
#'   uses \code{Matrix::crossprod} on the binarised matrix for efficient
#'   intersection/union computation. The weighted variant uses the identity
#'   \code{min(a,b) = (a + b - |a - b|) / 2} to leverage sparse matrix
#'   arithmetic.
#'
#' @param x A data matrix suitable for clustering (non-negative values expected)
#' @param weighted If \code{FALSE} (the default), compute binary Jaccard
#'   similarity. If \code{TRUE}, compute weighted (generalised) Jaccard
#'   similarity.
#' @param sparse Whether to return a sparse matrix (default \code{FALSE})
#' @param transpose When \code{FALSE} (the default) calculates similarity
#'   between columns. When \code{TRUE} calculates similarity between rows.
#'
#' @return A square similarity matrix with values in \code{[0,1]}.
#' @importFrom methods as
#' @export
#' @seealso \code{\link{cosine_sim}}, \code{\link{connectivity_similarity}}
#' @examples
#' da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
#' am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
#' # Binary Jaccard
#' jaccard_sim(am)
#' # Weighted Jaccard
#' jaccard_sim(am, weighted=TRUE)
jaccard_sim <- function(x, weighted = FALSE, sparse = FALSE, transpose = FALSE) {
  cx <- class(x)
  if (!is.matrix(x) && !isTRUE(attr(cx, "package") == "Matrix"))
    stop("I don't recognise that as a matrix!")
  if (!inherits(x, "dgCMatrix"))
    x <- as(x, "dgCMatrix")
  crossfun <- if(transpose) Matrix::tcrossprod else Matrix::crossprod
  n <- if(transpose) nrow(x) else ncol(x)
  nms <- if(transpose) rownames(x) else colnames(x)

  if (!weighted) {
    # Binary Jaccard: J = I / (s_i + s_j - I)
    # Binarise, crossprod for intersection, transform values in place
    b <- x
    b@x <- rep(1, length(b@x))
    A <- as(crossfun(b), "generalMatrix")
    sizes <- if(transpose) Matrix::rowSums(b) else diff(b@p)
    col_idx <- rep(seq_along(diff(A@p)), diff(A@p))
    row_idx <- A@i + 1L
    isect <- A@x
    A@x <- isect / (sizes[row_idx] + sizes[col_idx] - isect)
    Matrix::diag(A) <- 1
    sim <- A
  } else {
    # Weighted Jaccard: sum(min(a,b)) / sum(max(a,b))
    # For non-negative values: min(a,b) = sum_{t=1}^{max} I(a>=t)*I(b>=t)
    # Each threshold gives a sparse crossprod (same fast op as cosine).
    max_val <- if (length(x@x)) max(x@x) else 0
    if (max_val == 0) {
      sim <- Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n), x = 1,
                                  dims = c(n, n), dimnames = list(nms, nms))
    } else {
      cs <- if(transpose) Matrix::rowSums(x) else Matrix::colSums(x)
      vals <- x@x
      p <- x@p
      ri <- x@i
      col_idx <- rep(seq_len(ncol(x)), diff(p))
      min_sums <- matrix(0, n, n)
      for (t in seq_len(max_val)) {
        keep <- vals >= t
        bt <- Matrix::sparseMatrix(
          i = ri[keep] + 1L, j = col_idx[keep],
          x = 1, dims = x@Dim, dimnames = x@Dimnames)
        min_sums <- min_sums + as.matrix(crossfun(bt))
      }
      max_sums <- outer(cs, cs, "+") - min_sums
      sim <- min_sums / max_sums
      sim[is.nan(sim)] <- 0
      diag(sim) <- 1
      dimnames(sim) <- list(nms, nms)
      sim <- Matrix::Matrix(sim, sparse = TRUE)
    }
  }

  if (sparse) sim else as.matrix(sim)
}


#' Tanimoto (extended Jaccard) similarity for sparse or dense matrices
#'
#' @description Computes pairwise Tanimoto similarity between columns (or rows)
#'   of a matrix. Also known as the extended Jaccard coefficient, defined as:
#'   \code{dot(a,b) / (||a||^2 + ||b||^2 - dot(a,b))}. On binary data this
#'   reduces to the standard Jaccard index. Computed via a single
#'   \code{crossprod} call, so performance is comparable to \code{cosine_sim}.
#'
#' @param x A data matrix suitable for clustering (non-negative values expected)
#' @param sparse Whether to return a sparse matrix (default \code{FALSE})
#' @param transpose When \code{FALSE} (the default) calculates similarity
#'   between columns. When \code{TRUE} calculates similarity between rows.
#'
#' @return A square similarity matrix with values in \code{[0,1]}.
#' @importFrom methods as
#' @export
#' @seealso \code{\link{jaccard_sim}}, \code{\link{cosine_sim}},
#'   \code{\link{connectivity_similarity}}
#' @examples
#' da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
#' am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
#' tanimoto_sim(am)
tanimoto_sim <- function(x, sparse = FALSE, transpose = FALSE) {
  cx <- class(x)
  if (!is.matrix(x) && !isTRUE(attr(cx, "package") == "Matrix"))
    stop("I don't recognise that as a matrix!")
  if (!inherits(x, "dgCMatrix"))
    x <- as(x, "dgCMatrix")
  crossfun <- if (transpose) Matrix::tcrossprod else Matrix::crossprod

  # dot(a,b) for all pairs via crossprod (returns symmetric dsCMatrix)
  A <- crossfun(x)
  norms2 <- Matrix::diag(A)

  # T(a,b) = dot(a,b) / (||a||^2 + ||b||^2 - dot(a,b))
  # Work directly on the symmetric matrix (upper triangle stored)
  if (length(A@x)) {
    col_idx <- rep(seq_along(diff(A@p)), diff(A@p))
    row_idx <- A@i + 1L
    denom <- norms2[row_idx] + norms2[col_idx] - A@x
    nonzero <- denom != 0
    A@x[nonzero] <- A@x[nonzero] / denom[nonzero]
    A@x[!nonzero] <- 0
  }
  Matrix::diag(A) <- 1
  sim <- A

  if (sparse) sim else as.matrix(sim)
}
