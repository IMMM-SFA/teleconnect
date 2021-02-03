context("Verify constants")


test_that("verify constants values", {

  expect_equal("+proj=longlat +datum=WGS84 +no_defs", proj4_string)

  expect_equal(c(0, 81, 83, 88, 111), non_land_cdl_classes)

  expect_equal(c(82,121,122,123,124), developed_values)

  expect_equal(c(0,63:65,81,83:92,111,112,131:195), non_devcrop_class)

  expect_equal("Normal_Storage", hydro_variable)

  expect_equal(1233.48e-9, AF_to_BCM)

  expect_equal(0.00378541 * 1e-3, Mgallons_to_BCM) # Million gallons to billion cubic meters

  m2expect_equal(1000000, m2_to_km2)

  expect_equal(1000, m_to_km)

  expect_equal(0.001, mm_to_m)

  expect_equal(86400, day_to_sec)

  expect_equal(c(3589, 3850, 3678, 3582, 3615, 4226, 4227, 2336, 3448, 2162, 3474), missing_watersheds) # for runoff

  expect_equal(200, avg_wateruse_ltr_per_day)

  expect_equal(0.0438, MGD_to_m3sec)

  expect_equal(31.71, BCMyr_to_m3sec)

  expect_equal(373.36, BCMmonth_to_m3sec)

  expect_equal(1.15741e-8, ltrday_to_m3sec)

  expect_equal(0.0283168, cfs_to_m3sec)

  expect_equal("Q0001C", nhdplus_flow_metric)

})
