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
#' @param triangle If \code{TRUE}, return a \code{\link{dist}} object (lower
#'   triangle only, half memory). Default \code{FALSE}.
#' @param distance If \code{TRUE}, return distance (\code{1 - similarity})
#'   instead of similarity. Default \code{FALSE}.
#'
#' @return A square similarity matrix, or a \code{\link{dist}} object when
#'   \code{triangle = TRUE}.
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
                                    sparse = FALSE, transpose = FALSE,
                                    triangle = FALSE, distance = FALSE) {
  metric <- match.arg(metric)
  switch(metric,
    cosine = cosine_sim(x, sparse = sparse, transpose = transpose,
                        triangle = triangle, distance = distance),
    jaccard = jaccard_sim(x, weighted = FALSE, sparse = sparse, transpose = transpose,
                          triangle = triangle, distance = distance),
    weighted_jaccard = jaccard_sim(x, weighted = TRUE, sparse = sparse, transpose = transpose,
                                   triangle = triangle, distance = distance),
    tanimoto = tanimoto_sim(x, sparse = sparse, transpose = transpose,
                            triangle = triangle, distance = distance)
  )
}


warn_hourly_check <- memoise::memoise(function(msg) TRUE,
                                      ~ memoise::timeout(3600))

warn_hourly <- function(msg) {
  cached <- memoise::has_cache(warn_hourly_check)(msg)
  warn_hourly_check(msg)  # populate cache before warning
  if (!cached)
    warning(msg, call. = FALSE, immediate. = TRUE)
}

warn_natcpp <- function() {
  has <- requireNamespace("natcpp", quietly = TRUE)
  if (has) {
    msg <- paste0("natcpp (>= 0.3.0) is required for fast weighted Jaccard. ",
                  "You have ", utils::packageVersion("natcpp"), ". ",
                  "Please update with: ",
                  "install.packages('natcpp', repos='https://natverse.r-universe.dev')")
  } else {
    msg <- paste0("Install the natcpp package for much faster weighted Jaccard: ",
                  "install.packages('natcpp', repos='https://natverse.r-universe.dev')")
  }
  warn_hourly(msg)
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
#'   intersection/union computation. The weighted variant accumulates
#'   \code{sum(min(a,b))} across shared features for each pair; compiled C++
#'   backends are available for both dense and sparse output.
#'
#' @param x A data matrix suitable for clustering (non-negative values expected)
#' @param weighted If \code{FALSE} (the default), compute binary Jaccard
#'   similarity. If \code{TRUE}, compute weighted (generalised) Jaccard
#'   similarity.
#' @param weighted_method Algorithm to use when \code{weighted=TRUE}.
#'   \code{"auto"} (the default) selects the best available backend: compiled
#'   C++ dense for moderate output sizes, compiled C++ sparse when the output
#'   exceeds \code{10000 x 10000}, or pure R dense as a fallback when compiled
#'   code is unavailable. Explicit choices: \code{"cpp_dense"} and
#'   \code{"cpp_sparse"} for the compiled backends, \code{"dense"} and
#'   \code{"sparse"} for the pure R implementations.
#' @param sparse Whether to return a sparse matrix (default \code{FALSE})
#' @param transpose When \code{FALSE} (the default) calculates similarity
#'   between columns. When \code{TRUE} calculates similarity between rows.
#' @param triangle If \code{TRUE}, return a \code{\link{dist}} object (lower
#'   triangle only, half memory). Default \code{FALSE}.
#' @param distance If \code{TRUE}, return distance (\code{1 - similarity})
#'   instead of similarity. Default \code{FALSE}.
#'
#' @return A square similarity matrix, or a \code{\link{dist}} object when
#'   \code{triangle = TRUE}.
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
jaccard_sim <- function(x, weighted = FALSE, sparse = FALSE, transpose = FALSE,
                        weighted_method = c("auto", "cpp_dense", "cpp_sparse",
                                            "dense", "sparse"),
                        triangle = FALSE, distance = FALSE) {
  cx <- class(x)
  if (!is.matrix(x) && !isTRUE(attr(cx, "package") == "Matrix"))
    stop("I don't recognise that as a matrix!")
  if (!inherits(x, "dgCMatrix"))
    x <- as(x, "dgCMatrix")
  weighted_method <- match.arg(weighted_method)
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
    dimnames(A) <- list(nms, nms)
    return(sim_to_output(A, sparse = sparse, triangle = triangle,
                         distance = distance))
  } else {
    if (weighted_method == "auto") {
      has_cpp <- requireNamespace("natcpp", quietly = TRUE) &&
        utils::packageVersion("natcpp") >= "0.3.0"
      ncomp <- if (transpose) nrow(x) else ncol(x)
      if (has_cpp) {
        weighted_method <- if (ncomp > 10000L) "cpp_sparse" else "cpp_dense"
      } else {
        warn_natcpp()
        weighted_method <- "dense"
      }
    }
    sim <- switch(weighted_method,
      cpp_dense = jaccard_sim_weighted_cpp_dense(x, sparse = sparse, transpose = transpose,
                                                  triangle = triangle, distance = distance),
      cpp_sparse = jaccard_sim_weighted_cpp_sparse(x, sparse = sparse, transpose = transpose,
                                                    triangle = triangle, distance = distance),
      dense = jaccard_sim_weighted_dense_r(x, sparse = sparse, transpose = transpose,
                                           triangle = triangle, distance = distance),
      sparse = jaccard_sim_weighted_sparse_r(x, sparse = TRUE, transpose = transpose,
                                              triangle = triangle, distance = distance)
    )
  }

  if (inherits(sim, "dist")) sim
  else if (sparse) sim else as.matrix(sim)
}

jaccard_weighted_feature_view <- function(x, transpose = FALSE) {
  if (transpose) {
    feat_idx <- rep.int(seq_len(ncol(x)), diff(x@p))
    items <- x@i + 1L
    nfeat <- ncol(x)
    totals <- Matrix::rowSums(x)
  } else {
    feat_idx <- x@i + 1L
    items <- rep.int(seq_len(ncol(x)), diff(x@p))
    nfeat <- nrow(x)
    totals <- Matrix::colSums(x)
  }

  vals <- x@x
  counts <- tabulate(feat_idx, nbins = nfeat)
  offsets <- cumsum(c(0L, counts))
  fill <- offsets[-length(offsets)] + 1L
  feat_items <- integer(length(vals))
  feat_vals <- numeric(length(vals))

  for (k in seq_along(vals)) {
    f <- feat_idx[k]
    dest <- fill[f]
    feat_items[dest] <- items[k]
    feat_vals[dest] <- vals[k]
    fill[f] <- dest + 1L
  }

  list(
    ncomp = if (transpose) nrow(x) else ncol(x),
    nfeat = nfeat,
    totals = totals,
    offsets = offsets,
    items = feat_items,
    vals = feat_vals
  )
}

jaccard_sim_weighted_cpp_sparse <- function(x, sparse = TRUE, transpose = FALSE,
                                             triangle = FALSE, distance = FALSE) {
  n <- if (transpose) nrow(x) else ncol(x)
  nms <- if (transpose) rownames(x) else colnames(x)
  sim <- natcpp::c_weighted_jaccard_sparse(x, transpose = transpose,
                                           triangle = triangle, distance = distance)
  dimnames(sim) <- list(nms, nms)
  if (sparse) sim else as.matrix(sim)
}

jaccard_sim_weighted_cpp_dense <- function(x, sparse = FALSE, transpose = FALSE,
                                            triangle = FALSE, distance = FALSE) {
  n <- if (transpose) nrow(x) else ncol(x)
  nms <- if (transpose) rownames(x) else colnames(x)
  if (length(x@x) == 0L) {
    if (triangle) {
      d <- rep(if (distance) 1 else 0, n * (n - 1L) / 2L)
      return(structure(d, Size = n, Labels = nms, Diag = FALSE, Upper = FALSE,
                       class = "dist"))
    }
    sim <- if (distance) matrix(1, n, n) - diag(n) else diag(n)
    dimnames(sim) <- list(nms, nms)
    if (sparse) return(Matrix::Matrix(sim, sparse = TRUE))
    return(sim)
  }
  res <- natcpp::c_weighted_jaccard_dense(x, transpose = transpose,
                                          triangle = triangle, distance = distance)
  if (triangle) {
    attr(res, "Labels") <- nms
    return(res)
  }
  dimnames(res) <- list(nms, nms)
  if (sparse) Matrix::Matrix(res, sparse = TRUE) else res
}

jaccard_sim_weighted_sparse_r <- function(x, sparse = TRUE, transpose = FALSE,
                                          triangle = FALSE, distance = FALSE) {
  warn_natcpp()
  if (distance)
    warning("distance=TRUE with sparse output produces a mostly-dense matrix; ",
            "consider using the dense backend with triangle=TRUE instead.",
            call. = FALSE)
  n <- if (transpose) nrow(x) else ncol(x)
  nms <- if (transpose) rownames(x) else colnames(x)
  if (length(x@x) == 0L) {
    if (triangle) {
      d <- rep(if (distance) 1 else 0, n * (n - 1L) / 2L)
      return(structure(d, Size = n, Labels = nms, Diag = FALSE, Upper = FALSE,
                       class = "dist"))
    }
    sim_val <- if (distance) 0 else 1
    sim <- Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n), x = sim_val,
                                dims = c(n, n), dimnames = list(nms, nms))
    return(if (sparse) sim else as.matrix(sim))
  }

  fv <- jaccard_weighted_feature_view(x, transpose = transpose)
  counts <- diff(fv$offsets)
  npairs <- sum((counts * pmax.int(counts - 1L, 0L)) %/% 2L)

  ii <- integer(npairs)
  jj <- integer(npairs)
  mins <- numeric(npairs)
  pos <- 1L

  for (f in seq_len(fv$nfeat)) {
    start <- fv$offsets[f] + 1L
    end <- fv$offsets[f + 1L]
    k <- end - start + 1L
    if (k <= 1L) next

    if (k == 2L) {
      va <- fv$vals[start]
      vb <- fv$vals[end]
      ii[pos] <- fv$items[start]
      jj[pos] <- fv$items[end]
      mins[pos] <- if (va < vb) va else vb
      pos <- pos + 1L
      next
    }

    for (a in start:(end - 1L)) {
      ca <- fv$items[a]
      va <- fv$vals[a]
      for (b in (a + 1L):end) {
        vb <- fv$vals[b]
        ii[pos] <- ca
        jj[pos] <- fv$items[b]
        mins[pos] <- if (va < vb) va else vb
        pos <- pos + 1L
      }
    }
  }

  if (npairs == 0L) {
    sim <- Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n), x = 1,
                                dims = c(n, n), dimnames = list(nms, nms))
    return(if (sparse) sim else as.matrix(sim))
  }

  offdiag <- ii != jj
  min_sums <- Matrix::sparseMatrix(
    i = c(ii, jj[offdiag]),
    j = c(jj, ii[offdiag]),
    x = c(mins, mins[offdiag]),
    dims = c(n, n),
    dimnames = list(nms, nms)
  )

  cs <- fv$totals
  col_idx <- rep(seq_len(ncol(min_sums)), diff(min_sums@p))
  row_idx <- min_sums@i + 1L
  denom <- cs[row_idx] + cs[col_idx] - min_sums@x
  nz <- denom != 0
  min_sums@x[nz] <- min_sums@x[nz] / denom[nz]
  min_sums@x[!nz] <- 0
  Matrix::diag(min_sums) <- 1

  if (triangle) {
    dm <- as.matrix(min_sums)
    d <- dm[lower.tri(dm)]
    if (distance) d <- 1 - d
    return(structure(d, Size = n, Labels = nms, Diag = FALSE, Upper = FALSE,
                     class = "dist"))
  }

  if (distance)
    min_sums@x <- 1 - min_sums@x

  if (sparse) min_sums else as.matrix(min_sums)
}

jaccard_sim_weighted_dense_r <- function(x, sparse = FALSE, transpose = FALSE,
                                         triangle = FALSE, distance = FALSE) {
  n <- if (transpose) nrow(x) else ncol(x)
  nms <- if (transpose) rownames(x) else colnames(x)
  if (length(x@x) == 0L) {
    if (triangle) {
      d <- rep(if (distance) 1 else 0, n * (n - 1L) / 2L)
      return(structure(d, Size = n, Labels = nms, Diag = FALSE, Upper = FALSE,
                       class = "dist"))
    }
    sim <- if (distance) matrix(1, n, n) - diag(n) else diag(1, n)
    dimnames(sim) <- list(nms, nms)
    return(if (sparse) Matrix::Matrix(sim, sparse = TRUE) else sim)
  }

  fv <- jaccard_weighted_feature_view(x, transpose = transpose)
  ncomp <- fv$ncomp
  out <- matrix(0, ncomp, ncomp)

  for (f in seq_len(fv$nfeat)) {
    start <- fv$offsets[f] + 1L
    end <- fv$offsets[f + 1L]
    k <- end - start + 1L
    if (k <= 1L) next
    items_f <- fv$items[start:end]
    vals_f <- fv$vals[start:end]
    out[items_f, items_f] <- out[items_f, items_f] + outer(vals_f, vals_f, pmin)
  }

  # Convert min_sums to similarity column by column to avoid n*n temporary
  totals <- fv$totals
  for (j in seq_len(ncomp)) {
    denom <- totals + totals[j] - out[, j]
    nz <- denom != 0
    out[nz, j] <- out[nz, j] / denom[nz]
    out[!nz, j] <- 0
    out[j, j] <- 1
  }

  if (triangle) {
    d <- out[lower.tri(out)]
    if (distance) d <- 1 - d
    return(structure(d, Size = ncomp, Labels = nms, Diag = FALSE, Upper = FALSE,
                     class = "dist"))
  }

  if (distance) out <- 1 - out

  dimnames(out) <- list(nms, nms)
  if (sparse) Matrix::Matrix(out, sparse = TRUE) else out
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
#' @param triangle If \code{TRUE}, return a \code{\link{dist}} object (lower
#'   triangle only, half memory). Default \code{FALSE}.
#' @param distance If \code{TRUE}, return distance (\code{1 - similarity})
#'   instead of similarity. Default \code{FALSE}.
#'
#' @return A square similarity matrix, or a \code{\link{dist}} object when
#'   \code{triangle = TRUE}.
#' @importFrom methods as
#' @export
#' @seealso \code{\link{jaccard_sim}}, \code{\link{cosine_sim}},
#'   \code{\link{connectivity_similarity}}
#' @examples
#' da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
#' am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
#' tanimoto_sim(am)
tanimoto_sim <- function(x, sparse = FALSE, transpose = FALSE,
                         triangle = FALSE, distance = FALSE) {
  cx <- class(x)
  if (!is.matrix(x) && !isTRUE(attr(cx, "package") == "Matrix"))
    stop("I don't recognise that as a matrix!")
  if (!inherits(x, "dgCMatrix"))
    x <- as(x, "dgCMatrix")
  crossfun <- if (transpose) Matrix::tcrossprod else Matrix::crossprod
  nms <- if (transpose) rownames(x) else colnames(x)

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
  dimnames(A) <- list(nms, nms)

  sim_to_output(A, sparse = sparse, triangle = triangle, distance = distance)
}


#' @noRd
sim_to_output <- function(sim, sparse = FALSE, triangle = FALSE, distance = FALSE) {
  if (triangle) {
    if (inherits(sim, "dsCMatrix")) {
      d <- dsCMatrix_to_dist(sim, distance = distance)
    } else {
      dm <- as.matrix(sim)
      d <- dm[lower.tri(dm)]
      if (distance) d <- 1 - d
      nms <- rownames(dm)
      d <- structure(d, Size = nrow(dm), Labels = nms, Diag = FALSE,
                     Upper = FALSE, class = "dist")
    }
    return(d)
  }
  if (distance) {
    sim <- if (is.matrix(sim)) 1 - sim else 1 - as.matrix(sim)
    return(sim)
  }
  if (sparse) sim else as.matrix(sim)
}

#' Extract lower triangle from a dsCMatrix as a dist object without densifying
#' @noRd
dsCMatrix_to_dist <- function(x, distance = FALSE) {
  n <- nrow(x)
  fill_val <- if (distance) 1 else 0
  d <- rep(fill_val, n * (n - 1L) / 2L)

  col0 <- rep(seq_along(diff(x@p)), diff(x@p)) - 1L
  row0 <- x@i
  vals <- x@x

  # Keep only off-diagonal entries (upper triangle: row0 < col0)
  offdiag <- row0 != col0
  row0 <- row0[offdiag]
  col0 <- col0[offdiag]
  vals <- vals[offdiag]

  if (length(vals)) {
    # dsCMatrix stores upper triangle (row < col).
    # dist layout is lower triangle column-major: for pair (i,j), i > j:
    #   pos = (j-1)*n - j*(j-1)/2 + (i-j)   [1-based]
    # Here row0 < col0 (0-based), so i=col0+1, j=row0+1:
    r <- row0 + 1L
    cc <- col0 + 1L
    pos <- (r - 1L) * n - r * (r - 1L) / 2L + (cc - r)
    d[pos] <- if (distance) 1 - vals else vals
  }

  structure(d, Size = n, Labels = rownames(x),
            Diag = FALSE, Upper = FALSE, class = "dist")
}

#' Convert a symmetric sparse similarity matrix to a dist object
#'
#' Efficiently converts a symmetric sparse similarity matrix (\code{dsCMatrix})
#' to a \code{\link{dist}} object without materialising the full dense matrix.
#' Similarity values are converted to distances as \code{1 - similarity}.
#'
#' @param x A symmetric sparse matrix of class \code{dsCMatrix} with values in
#'   \code{[0,1]}, as returned by e.g. \code{\link{jaccard_sim}} with
#'   \code{sparse=TRUE}.
#' @return A \code{\link{dist}} object
#' @export
sim2dist <- function(x) {
  if (!inherits(x, "dsCMatrix"))
    stop("x must be a dsCMatrix (symmetric sparse column-compressed matrix)")
  n <- nrow(x)
  # Pre-fill with distance=1 (similarity=0 for absent pairs)
  d <- rep(1, n * (n - 1L) / 2L)

  # Extract upper triangle entries (excluding diagonal)
  col0 <- rep(seq_along(diff(x@p)), diff(x@p)) - 1L
  row0 <- x@i
  vals <- x@x

  offdiag <- row0 != col0
  row0 <- row0[offdiag]
  col0 <- col0[offdiag]
  vals <- vals[offdiag]

  # Map upper triangle (r,c) r<c to dist position (1-based)
  r <- row0 + 1L
  cc <- col0 + 1L
  pos <- (r - 1L) * n - r * (r - 1L) / 2L + (cc - r)
  d[pos] <- 1 - vals

  structure(d, Size = n, Labels = rownames(x),
            Diag = FALSE, Upper = FALSE, class = "dist")
}
