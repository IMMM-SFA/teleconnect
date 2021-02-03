context("Analysis functions")


#' Create a sample polygon sf object to use for tests
#'
#' @return sf polygon object
#' @importFrom sf st_polygon
#'
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
sample_polygon <- function() {

  # create a sample polygon to use for tests
  outer = matrix(c(0,0,10,0,10,10,0,10,0,0),ncol=2, byrow=TRUE)
  hole1 = matrix(c(1,1,1,2,2,2,2,1,1,1),ncol=2, byrow=TRUE)
  hole2 = matrix(c(5,5,5,6,6,6,6,5,5,5),ncol=2, byrow=TRUE)
  pts = list(outer, hole1, hole2)

  return(pl1 = st_polygon(pts))

}


test_that("init_bbox() functionality", {

  # build sample polygon sf object
  ply <- sample_polygon()

  # build what to expect as an output
  comp_bbox <- list()
  comp_bbox[[ "xmin" ]] <- -1
  comp_bbox[[ "ymin" ]] <- -1
  comp_bbox[[ "ymax" ]] <- 11
  comp_bbox[[ "xmin_source" ]] <- 0
  comp_bbox[[ "xmax_source" ]] <- 10

  # run the function and return the output list
  sim_bbox <- init_bbox(ply)

  # ensure that the function produces what is expected
  expect_equal(comp_bbox, sim_bbox)

})


test_that("poly_intersect() functionality", {

  # build sample polygon sf object
  ply <- sample_polygon()

  # create a bounding box attribute list
  sim_bbox <- init_bbox(ply)

  # load expected output data
  comp_ply <- readRDS("data/comp_poly_intersect.rds")

  # generate polygon intersection data
  sim_ply <- poly_intersect(sim_bbox$xmin, sim_bbox$xmax, sim_bbox$ymin, sim_bbox$ymax, ply)

  # ensure that the function produces what is expected
  expect_equal(comp_ply, sim_ply)

})
