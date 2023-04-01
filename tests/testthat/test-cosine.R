test_that("cosine_sim works", {
  da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
  am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
  expect_snapshot(cosine_sim(am))

  expect_snapshot(pm <- prepare_cosine_matrix(cosine_sim(am, transpose = T)))

  expect_silent(hm <- cosine_heatmap(pm))
  expect_true(inherits(hm$Colv, "dendrogram"))

  expect_equal(hm, cosine_heatmap(pm, heatmap=stats::heatmap))

  expect_silent(cl <- cosine_heatmap(pm, heatmap=FALSE))
  expect_true(inherits(cl, "hclust"))
})

