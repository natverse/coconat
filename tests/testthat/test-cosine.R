test_that("cosine_sim works", {
  da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
  am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
  expect_snapshot(cosine_sim(am))
})
