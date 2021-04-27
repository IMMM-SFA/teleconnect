context("Validate IO functionality")


test_that("get_cities functionality", {

  expected_cols <- c("city", "state", "city_uid", "intake", "DVSN_ID", "DVTR_ID", "city_state", "key_watershed")
  expected_dim <- c(888, 8)

  # generate output
  df <- get_cities()

  # check equality
  expect_equal(expected_cols, colnames(df))
  expect_equal(expected_dim, dim(df))

})


test_that("get_crop_mapping functionality", {

  expected_cols <- c("item", "GTAP_crop", "GCAM_commodity", "CROSIT_crop", "CROSIT_cropID",
                     "IFA2002_crop", "IFA_commodity", "MIRCA_crop_name", "MIRCA_crop",
                     "MIRCA_Crop24and26", "MH_crop",   "MH2011_crop", "MH2014_proxy",
                     "GTAP_use",  "LPJmL_crop", "GEPIC_crop",  "Pegasus_crop",  "C3avg_include")
  expected_dim <- c(182, 18)

  df <- get_crop_mapping()

  # check equality
  expect_equal(expected_cols, colnames(df))
  expect_equal(expected_dim, dim(df))

})
