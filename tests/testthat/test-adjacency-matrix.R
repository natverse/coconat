test_that("partner_summary2adjacency_matrix works", {
  # hbconn=neuprintr::neuprint_login(dataset='hemibrain:v1.2.1')
  # da2ds15=neuprintr::neuprint_connection_table('DA2_lPN', details=TRUE, partners='out', threshold=15, conn = hbconn)
  # saveRDS(da2ds15, file = 'inst/sampledata/da2ds15.rds', version=2)

  da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))

  am <- partner_summary2adjacency_matrix(da2ds15, inputcol = 'bodyid', outputcol = 'partner')
  expect_true(inherits(am,'dgCMatrix'))
  expect_snapshot(am2 <- partner_summary2adjacency_matrix(da2ds15, inputcol = 'bodyid', outputcol = 'partner', sparse = F))
  expect_snapshot(am3 <- partner_summary2adjacency_matrix(da2ds15, inputcol = 'pre_id', outputcol = 'post_id', sparse = F))

  bl=c(`1796817841` = 380, `1797505019` = 207, `1796818119` = 489,
       `1827516355` = 197, `818983130` = 380)

  expect_equal(Matrix::rowSums(am), bl)
  expect_equal(rowSums(am2), Matrix::rowSums(am))
})
