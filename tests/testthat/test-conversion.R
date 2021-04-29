context("Validate conversion functions")


test_that("raster_to_polygon valid output", {

  # build raster object
  rast_obj <- raster()

  # load comparison data of expected output
  comp_ply <- readRDS("./data/raster_to_polygon.rds")

  # create output from function
  sim_ply <- raster_to_polygon(rast_obj)

  # validate match
  expect_equal(comp_ply, sim_ply)

})
