context("Validate count_watershed_data output")


test_that("count_watershed_data valid output", {

  data_dir <- '/null'
  cities <- c("Portland | OR")

  expect_error(
    df <- sup(count_watershed_data(data_dir = data_dir, cities = cities))
  )

})
