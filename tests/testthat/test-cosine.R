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
  # neuprintr::neuprint_login(dataset='hemibrain:v1.2.1')
  # clm=neuprintr::neuprint_get_meta(cl$labels)
  # saveRDS(clm, file = 'inst/sampledata/hemibrain-metadata.rds')
  clm=readRDS(system.file('sampledata/hemibrain-metadata.rds', package = 'coconat'))
  clm2=add_cluster_info(clm, cl, h=1)
  expect_equal(clm2$group_h1,
               c(2L, 2L, 5L, 1L, 4L, 1L, 3L, 1L, 3L, 3L, 3L, 2L, 1L, 1L, 1L,
                 1L, 1L, 1L, 4L, 2L, 2L, 1L, 5L, 4L, 3L, 3L, 3L, 1L, 1L, 1L, 3L,
                 2L, 1L, 1L, 3L, 1L, 4L, 4L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 6L,
                 1L, 5L, 5L, 5L, 6L, 4L, 3L, 3L, 3L, 3L, 2L, 6L, 2L, 2L, 1L, 1L,
                 1L, 1L, 5L, 5L, 4L, 4L, 3L, 2L, 2L, 2L, 1L, 1L, 1L, 5L, 5L, 4L,
                 4L, 4L, 3L))
})

