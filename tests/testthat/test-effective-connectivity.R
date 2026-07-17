test_that("grouping_matrix aggregates columns", {
  M <- Matrix::Matrix(c(1,0,2, 0,3,0), nrow=2, byrow=TRUE,
    dimnames=list(c('a','b'), c('x','y','z')))
  # parallel grouping vector
  G <- grouping_matrix(colnames(M), c('g1','g1','g2'))
  expect_equal(dim(G), c(3L, 2L))
  expect_equal(colnames(G), c('g1','g2'))
  # x,y -> g1 ; z -> g2  (a: 1+0=1, 2 ; b: 0+3=3, 0)
  expect_equal(as.matrix(M %*% G),
    matrix(c(1,2, 3,0), nrow=2, byrow=TRUE,
      dimnames=list(c('a','b'), c('g1','g2'))))

  # named grouping vector is looked up by ids (order independent)
  Gn <- grouping_matrix(colnames(M), c(z='g2', x='g1', y='g1'))
  expect_equal(as.matrix(Gn), as.matrix(G))

  # NA group is dropped entirely
  Gna <- grouping_matrix(colnames(M), c('g1', NA, 'g2'))
  expect_equal(colnames(Gna), c('g1','g2'))
  # y (the NA-grouped column) drops out: g1 gets only x, g2 gets only z
  expect_equal(as.matrix(M %*% Gna),
    matrix(c(1,2, 0,0), nrow=2, byrow=TRUE,
      dimnames=list(c('a','b'), c('g1','g2'))))

  expect_error(grouping_matrix(colnames(M), c('g1','g2')))
})

test_that("effective_connectivity matches a hand-computed pathway", {
  # Q -> I -> T
  A1 <- Matrix::Matrix(c(2,0, 1,3), nrow=2, byrow=TRUE,
    dimnames=list(c('q1','q2'), c('i1','i2')))
  A2 <- Matrix::Matrix(c(4,0, 2,2), nrow=2, byrow=TRUE,
    dimnames=list(c('i1','i2'), c('t1','t2')))
  # colScaleM(A1): cols sum to 1 -> i1=(2/3,1/3), i2=(0,1)
  # colScaleM(A2): t1=(2/3,1/3), t2=(0,1)
  # effective q1: t1=4/9, t2=0 ; q2: t1=5/9, t2=1
  eff <- effective_connectivity(list(A1, A2))
  expect_equal(as.matrix(eff),
    matrix(c(4/9, 0, 5/9, 1), nrow=2, byrow=TRUE,
      dimnames=list(c('q1','q2'), c('t1','t2'))))

  # grouping the outputs to a single type sums the columns
  effg <- effective_connectivity(list(A1, A2), group=c(t1='T', t2='T'))
  expect_equal(as.matrix(effg),
    matrix(c(4/9, 14/9), ncol=1,
      dimnames=list(c('q1','q2'), 'T')))
})

test_that("effective_connectivity is associative wrt grouping", {
  A1 <- Matrix::Matrix(c(2,0, 1,3), nrow=2, byrow=TRUE,
    dimnames=list(c('q1','q2'), c('i1','i2')))
  A2 <- Matrix::Matrix(c(4,0, 2,2), nrow=2, byrow=TRUE,
    dimnames=list(c('i1','i2'), c('t1','t2')))
  g <- c(t1='T', t2='T')
  # grouping inside == grouping the ungrouped product afterwards
  eff <- effective_connectivity(list(A1, A2))
  manual <- eff %*% grouping_matrix(colnames(eff), g)
  expect_equal(effective_connectivity(list(A1, A2), group=g), manual)
})

test_that("effective_connectivity aligns by name and treats dead ends as zero", {
  A1 <- Matrix::Matrix(c(2,0, 1,3), nrow=2, byrow=TRUE,
    dimnames=list(c('q1','q2'), c('i1','i2')))
  # i2 is absent from the rows of A2 -> its paths are dead ends (zero)
  A2 <- Matrix::Matrix(c(4,0), nrow=1,
    dimnames=list('i1', c('t1','t2')))
  eff <- effective_connectivity(list(A1, A2))
  # only i1 contributes: q1 via i1=2/3, q2 via i1=1/3 ; A2n col t1 = i1 -> 1
  expect_equal(as.matrix(eff),
    matrix(c(2/3, 0, 1/3, 0), nrow=2, byrow=TRUE,
      dimnames=list(c('q1','q2'), c('t1','t2'))))

  # rows of A2 reordered relative to cols of A1 -> still correct
  A2b <- Matrix::Matrix(c(2,2, 4,0), nrow=2, byrow=TRUE,
    dimnames=list(c('i2','i1'), c('t1','t2')))
  effb <- effective_connectivity(list(A1, A2b))
  expect_equal(as.matrix(effb),
    matrix(c(4/9, 0, 5/9, 1), nrow=2, byrow=TRUE,
      dimnames=list(c('q1','q2'), c('t1','t2'))))
})

test_that("effective_connectivity single matrix just normalises", {
  A1 <- Matrix::Matrix(c(2,0, 1,3), nrow=2, byrow=TRUE,
    dimnames=list(c('q1','q2'), c('i1','i2')))
  expect_equal(effective_connectivity(list(A1)), colScaleM(A1))
  expect_equal(effective_connectivity(list(A1), normalise=FALSE), A1)
  expect_error(effective_connectivity(A1))
})
