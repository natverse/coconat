test_that("dataset support", {
  expect_silent(register_dataset('flywire', shortname = 'fw', species = 'Drosophila melanogaster', sex = 'F', namespace = 'testthat'))

  expect_equal(dataset_names(namespace = 'testthat'), "flywire")
  expect_silent(register_dataset('hemibrain', shortname = 'hb', species = 'Drosophila melanogaster', sex = 'F', namespace = 'testthat'))
  expect_equal(dataset_names(namespace = 'testthat'), c("flywire", "hemibrain"))
  expect_equal(dataset_names("h", namespace = 'testthat'), "hemibrain")

  expect_equal(dataset_names("h", namespace = 'testthat', return.short = T), "hb")

  expect_silent(register_dataset('flywirex', shortname = 'fx', namespace = 'testthat', inherits = 'flywire'))
  expect_equal(dataset_details('flywirex', namespace = 'testthat')$sex, 'F')

  remove_namespace('testthat')
})
