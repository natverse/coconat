test_that("id2har works", {
  expect_equal(id2char(1000), "1000")
  expect_equal(id2char(NA), NA_character_)
  expect_equal(id2char("1000"), "1000")
  expect_equal(id2char(1e5), "100000")
  expect_error(id2char(TRUE))
})
