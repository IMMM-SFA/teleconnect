context("Analysis functions")


test_that("init_bbox() functionality", {

  # load sample shapefile from "sf"
  ply <- sf::st_read(system.file("shape/nc.shp", package="sf"))

  # build what to expect as an output
  comp_bbox <- list()
  comp_bbox[[ "xmin" ]] <- -84.59462
  comp_bbox[[ "ymin" ]] <- 33.61123
  comp_bbox[[ "ymax" ]] <- 36.86041
  comp_bbox[[ "xmin_source" ]] <- -84.32385
  comp_bbox[[ "xmax_source" ]] <- -75.45698

  # run the function and return the output list
  sim_bbox <- init_bbox(ply)

  # ensure that the function produces what is expected
  expect_equal(comp_bbox, sim_bbox, tolerance=1e-5)

})


test_that("poly_intersect() functionality", {

  # build sample polygon sf object
  ply <- sample_polygon()

  # create a bounding box attribute list
  sim_bbox <- init_bbox(ply)




})





