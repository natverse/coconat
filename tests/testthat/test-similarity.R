# Naive reference implementation for verification
naive_jaccard <- function(x, weighted = FALSE) {
  n <- ncol(x)
  sim <- matrix(0, n, n, dimnames = list(colnames(x), colnames(x)))
  for (i in seq_len(n)) {
    for (j in seq(i, n)) {
      a <- x[, i]
      b <- x[, j]
      if (!weighted) {
        ab <- (a != 0) & (b != 0)
        aub <- (a != 0) | (b != 0)
        sim[i, j] <- sim[j, i] <- if (sum(aub) == 0) 0 else sum(ab) / sum(aub)
      } else {
        min_sum <- sum(pmin(a, b))
        max_sum <- sum(pmax(a, b))
        sim[i, j] <- sim[j, i] <- if (max_sum == 0) 0 else min_sum / max_sum
      }
    }
  }
  diag(sim) <- 1
  sim
}

test_that("jaccard_sim matches naive implementation on sample data", {
  da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
  am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')

  # Binary
  jb <- jaccard_sim(am)
  jb_naive <- naive_jaccard(as.matrix(am), weighted = FALSE)
  expect_equal(jb, jb_naive)

  # Weighted
  jw <- jaccard_sim(am, weighted = TRUE)
  jw_naive <- naive_jaccard(as.matrix(am), weighted = TRUE)
  expect_equal(jw, jw_naive)
})

test_that("jaccard_sim hand-computed values are correct", {
  # 3 rows (partners), 3 columns (neurons) with known structure
  #       n1  n2  n3
  # p1  [ 4   2   0 ]
  # p2  [ 0   3   3 ]
  # p3  [ 1   1   0 ]
  m <- matrix(c(4,0,1, 2,3,1, 0,3,0), nrow = 3, ncol = 3)
  colnames(m) <- c("n1", "n2", "n3")

  # Binary Jaccard between columns (nonzero patterns):
  # n1: {p1, p3}        n2: {p1, p2, p3}    n3: {p2}
  # J(n1,n2) = |{p1,p3}| / |{p1,p2,p3}| = 2/3
  # J(n1,n3) = |{}| / |{p1,p2,p3}| = 0/3 = 0
  # J(n2,n3) = |{p2}| / |{p1,p2,p3}| = 1/3
  jb <- jaccard_sim(m)
  expect_equal(jb["n1","n2"], 2/3)
  expect_equal(jb["n1","n3"], 0)
  expect_equal(jb["n2","n3"], 1/3)

  # Weighted Jaccard:
  # J_w(n1,n2) = (min(4,2)+min(0,3)+min(1,1)) / (max(4,2)+max(0,3)+max(1,1))
  #            = (2+0+1) / (4+3+1) = 3/8
  # J_w(n1,n3) = (min(4,0)+min(0,3)+min(1,0)) / (max(4,0)+max(0,3)+max(1,0))
  #            = (0+0+0) / (4+3+1) = 0/8 = 0
  # J_w(n2,n3) = (min(2,0)+min(3,3)+min(1,0)) / (max(2,0)+max(3,3)+max(1,0))
  #            = (0+3+0) / (2+3+1) = 3/6 = 1/2
  jw <- jaccard_sim(m, weighted = TRUE)
  expect_equal(jw["n1","n2"], 3/8)
  expect_equal(jw["n1","n3"], 0)
  expect_equal(jw["n2","n3"], 1/2)
})

test_that("jaccard_sim: identical columns give 1, disjoint columns give 0", {
  #       n1  n2  n3
  # p1  [ 5   5   0 ]
  # p2  [ 3   3   0 ]
  # p3  [ 0   0   7 ]
  m <- Matrix::Matrix(c(5,3,0, 5,3,0, 0,0,7), nrow = 3, ncol = 3, sparse = TRUE)
  colnames(m) <- c("n1", "n2", "n3")

  # n1 and n2 are identical → both Jaccard variants should be 1
  # n1/n2 and n3 are completely disjoint → both should be 0
  jb <- jaccard_sim(m)
  expect_equal(jb["n1","n2"], 1)
  expect_equal(jb["n1","n3"], 0)
  expect_equal(jb["n2","n3"], 0)

  jw <- jaccard_sim(m, weighted = TRUE)
  expect_equal(jw["n1","n2"], 1)
  expect_equal(jw["n1","n3"], 0)
  expect_equal(jw["n2","n3"], 0)
})

test_that("jaccard_sim: proportional columns differ from cosine", {
  # a = 2*b: cosine similarity = 1, but weighted Jaccard = 0.5
  m <- matrix(c(2,4,6, 1,2,3), nrow = 3, ncol = 2)
  colnames(m) <- c("a", "b")

  # Cosine should be 1 (proportional vectors have cosine sim = 1)
  cs <- cosine_sim(m)
  expect_equal(cs["a","b"], 1)

  # Binary Jaccard should be 1 (same nonzero pattern)
  jb <- jaccard_sim(m)
  expect_equal(jb["a","b"], 1)

  # Weighted Jaccard: sum(min)/sum(max) = (1+2+3)/(2+4+6) = 6/12 = 0.5
  jw <- jaccard_sim(m, weighted = TRUE)
  expect_equal(jw["a","b"], 0.5)
})

test_that("jaccard_sim: binary matrix gives same result for binary and weighted", {
  m <- Matrix::Matrix(c(1,0,1,0, 0,1,1,0, 1,1,0,1), nrow = 4, ncol = 3, sparse = TRUE)
  colnames(m) <- c("a", "b", "c")

  jb <- jaccard_sim(m)
  jw <- jaccard_sim(m, weighted = TRUE)
  expect_equal(jb, jw)
})

test_that("jaccard_sim works with dense matrix input", {
  m_sparse <- Matrix::Matrix(c(4,0,1, 2,3,1, 0,3,0), nrow = 3, ncol = 3, sparse = TRUE)
  colnames(m_sparse) <- c("n1", "n2", "n3")
  m_dense <- as.matrix(m_sparse)

  expect_equal(jaccard_sim(m_dense), jaccard_sim(m_sparse))
  expect_equal(jaccard_sim(m_dense, weighted = TRUE),
               jaccard_sim(m_sparse, weighted = TRUE))
})

test_that("jaccard_sim sparse output parameter works", {
  m <- Matrix::Matrix(c(4,0,1, 2,3,1, 0,3,0), nrow = 3, ncol = 3, sparse = TRUE)
  colnames(m) <- c("n1", "n2", "n3")

  js <- jaccard_sim(m, sparse = TRUE)
  expect_true(inherits(js, "Matrix"))
  expect_equal(as.matrix(js), jaccard_sim(m, sparse = FALSE))

  jws <- jaccard_sim(m, weighted = TRUE, sparse = TRUE)
  expect_true(inherits(jws, "Matrix"))
  expect_equal(as.matrix(jws), jaccard_sim(m, weighted = TRUE, sparse = FALSE))
})

test_that("all weighted jaccard methods agree", {
  m_sparse <- Matrix::Matrix(c(4,0,1, 2,3,1, 0,3,0), nrow = 3, ncol = 3, sparse = TRUE)
  colnames(m_sparse) <- c("n1", "n2", "n3")
  rownames(m_sparse) <- c("p1", "p2", "p3")
  m_dense <- as.matrix(m_sparse)

  ref <- jaccard_sim(m_sparse, weighted = TRUE, weighted_method = "dense")
  ref_t <- jaccard_sim(m_sparse, weighted = TRUE, transpose = TRUE,
                       weighted_method = "dense")

  for (method in c("sparse", "dense", "cpp_dense", "cpp_sparse")) {
    suppressWarnings({
      # Default (dense return, transpose=FALSE)
      expect_equal(
        jaccard_sim(m_sparse, weighted = TRUE, weighted_method = method),
        ref, info = paste(method, "sparse input")
      )
      # Dense matrix input
      expect_equal(
        jaccard_sim(m_dense, weighted = TRUE, weighted_method = method),
        ref, info = paste(method, "dense input")
      )
      # transpose=TRUE
      expect_equal(
        jaccard_sim(m_sparse, weighted = TRUE, transpose = TRUE,
                    weighted_method = method),
        ref_t, info = paste(method, "transpose")
      )
      # sparse return
      expect_equal(
        as.matrix(jaccard_sim(m_sparse, weighted = TRUE, sparse = TRUE,
                              weighted_method = method)),
        ref, info = paste(method, "sparse return")
      )
    })
  }
})

test_that("jaccard_sim transpose parameter works", {
  m <- Matrix::Matrix(c(4,0,1, 2,3,1, 0,3,0), nrow = 3, ncol = 3, sparse = TRUE)
  colnames(m) <- c("n1", "n2", "n3")
  rownames(m) <- c("p1", "p2", "p3")

  # Default: similarity between columns (3x3)
  jc <- jaccard_sim(m)
  expect_equal(dim(jc), c(3L, 3L))
  expect_equal(rownames(jc), c("n1", "n2", "n3"))

  # Transpose: similarity between rows (3x3 here too, but different values)
  jt <- jaccard_sim(m, transpose = TRUE)
  expect_equal(dim(jt), c(3L, 3L))
  expect_equal(rownames(jt), c("p1", "p2", "p3"))

  # Transpose should give same result as manually transposing
  expect_equal(jt, jaccard_sim(Matrix::t(m)))
})

test_that("jaccard_sim snapshot on sample data", {
  da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
  am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')

  expect_snapshot(jaccard_sim(am))
  expect_snapshot(jaccard_sim(am, weighted = TRUE))
})

test_that("jaccard_sim handles edge cases", {
  # Single column matrix
  m1 <- Matrix::Matrix(c(1, 0, 3), ncol = 1, sparse = TRUE)
  colnames(m1) <- "a"
  j1 <- jaccard_sim(m1)
  expect_equal(j1, matrix(1, 1, 1, dimnames = list("a", "a")))

  j1w <- jaccard_sim(m1, weighted = TRUE)
  expect_equal(j1w, matrix(1, 1, 1, dimnames = list("a", "a")))

  # All-zero columns: diagonal should still be 1, off-diagonal 0
  m2 <- Matrix::Matrix(0, nrow = 3, ncol = 2, sparse = TRUE)
  colnames(m2) <- c("a", "b")
  j2 <- jaccard_sim(m2)
  expect_true(all(diag(j2) == 1))
  expect_equal(j2["a", "b"], 0)

  j2w <- jaccard_sim(m2, weighted = TRUE)
  expect_true(all(diag(j2w) == 1))
  expect_equal(j2w["a", "b"], 0)

  # Rejects non-matrix input
  expect_error(jaccard_sim(1:10), "matrix")
  expect_error(jaccard_sim(data.frame(a = 1:3)), "matrix")
})

test_that("tanimoto_sim matches naive implementation", {
  naive_tanimoto <- function(x) {
    n <- ncol(x)
    sim <- matrix(0, n, n, dimnames = list(colnames(x), colnames(x)))
    for (i in seq_len(n)) {
      for (j in seq(i, n)) {
        a <- x[, i]; b <- x[, j]
        d <- sum(a * b)
        denom <- sum(a^2) + sum(b^2) - d
        sim[i, j] <- sim[j, i] <- if (denom == 0) 0 else d / denom
      }
    }
    diag(sim) <- 1
    sim
  }

  da2ds15 <- readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
  am <- partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
  expect_equal(tanimoto_sim(am), naive_tanimoto(as.matrix(am)))
})

test_that("tanimoto_sim hand-computed values are correct", {
  m <- matrix(c(4,0,1, 2,3,1, 0,3,0), nrow = 3, ncol = 3)
  colnames(m) <- c("n1", "n2", "n3")

  # T(n1,n2) = dot(n1,n2) / (||n1||^2 + ||n2||^2 - dot(n1,n2))
  # dot = 4*2 + 0*3 + 1*1 = 9
  # ||n1||^2 = 16+0+1 = 17, ||n2||^2 = 4+9+1 = 14
  # T = 9 / (17 + 14 - 9) = 9/22
  ts <- tanimoto_sim(m)
  expect_equal(ts["n1","n2"], 9/22)
  expect_equal(ts["n1","n3"], 0)  # disjoint
  # dot(n2,n3) = 0+9+0 = 9, ||n3||^2 = 9
  # T = 9 / (14 + 9 - 9) = 9/14
  expect_equal(ts["n2","n3"], 9/14)
})

test_that("tanimoto_sim equals binary jaccard on 0/1 data", {
  m <- Matrix::Matrix(c(1,0,1,0, 0,1,1,0, 1,1,0,1), nrow = 4, ncol = 3, sparse = TRUE)
  colnames(m) <- c("a", "b", "c")
  expect_equal(tanimoto_sim(m), jaccard_sim(m))
})

test_that("tanimoto_sim handles edge cases", {
  # Single column
  m1 <- Matrix::Matrix(c(1, 0, 3), ncol = 1, sparse = TRUE)
  colnames(m1) <- "a"
  expect_equal(tanimoto_sim(m1), matrix(1, 1, 1, dimnames = list("a", "a")))

  # All-zero columns
  m2 <- Matrix::Matrix(0, nrow = 3, ncol = 2, sparse = TRUE)
  colnames(m2) <- c("a", "b")
  t2 <- tanimoto_sim(m2)
  expect_true(all(diag(t2) == 1))
  expect_equal(t2["a", "b"], 0)

  # Rejects non-matrix input
  expect_error(tanimoto_sim(1:10), "matrix")
})

test_that("tanimoto_sim sparse and transpose parameters work", {
  m <- Matrix::Matrix(c(4,0,1, 2,3,1, 0,3,0), nrow = 3, ncol = 3, sparse = TRUE)
  colnames(m) <- c("n1", "n2", "n3")
  rownames(m) <- c("p1", "p2", "p3")

  ts <- tanimoto_sim(m, sparse = TRUE)
  expect_true(inherits(ts, "Matrix"))
  expect_equal(as.matrix(ts), tanimoto_sim(m))

  tt <- tanimoto_sim(m, transpose = TRUE)
  expect_equal(dim(tt), c(3L, 3L))
  expect_equal(rownames(tt), c("p1", "p2", "p3"))
  expect_equal(tt, tanimoto_sim(Matrix::t(m)))
})

test_that("connectivity_similarity dispatches correctly", {
  da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
  am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')

  expect_equal(connectivity_similarity(am, metric = "cosine"), cosine_sim(am))
  expect_equal(connectivity_similarity(am, metric = "jaccard"), jaccard_sim(am))
  expect_equal(connectivity_similarity(am, metric = "weighted_jaccard"),
               jaccard_sim(am, weighted = TRUE))
  expect_equal(connectivity_similarity(am, metric = "tanimoto"), tanimoto_sim(am))
})

test_that("prepare_similarity_matrix is an alias for prepare_cosine_matrix", {
  da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
  am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
  cm <- cosine_sim(am, transpose = TRUE)
  expect_equal(prepare_similarity_matrix(cm), prepare_cosine_matrix(cm))
})

test_that("connectivity_heatmap works with distfun", {
  da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
  am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
  jw <- jaccard_sim(am, weighted = TRUE, transpose = TRUE)
  pm <- prepare_similarity_matrix(jw)

  # Default distfun (1-x)
  expect_silent(cl <- connectivity_heatmap(pm, heatmap = FALSE))
  expect_true(inherits(cl, "hclust"))

  # Custom distfun
  expect_silent(cl2 <- connectivity_heatmap(pm, heatmap = FALSE,
    distfun = function(x) as.dist(1 - x)))
  expect_equal(cl, cl2)
})
