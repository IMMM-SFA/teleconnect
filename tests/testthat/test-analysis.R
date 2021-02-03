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


# build sample polygon sf object
test_ply <- sample_polygon()

# create a bounding box attribute list
test_bbox <- init_bbox(test_ply)


test_that("init_bbox() functionality", {

  # build what to expect as an output
  comp_bbox <- list()
  comp_bbox[[ "xmin" ]] <- -1
  comp_bbox[[ "ymin" ]] <- -1
  comp_bbox[[ "ymax" ]] <- 11
  comp_bbox[[ "xmin_source" ]] <- 0
  comp_bbox[[ "xmax_source" ]] <- 10

  # run the function and return the output list
  sim_bbox <- init_bbox(test_ply)

  # ensure that the function produces what is expected
  expect_equal(comp_bbox, sim_bbox)

})


test_that("poly_intersect() functionality", {

  # load expected output data
  comp_ply <- readRDS("data/comp_poly_intersect.rds")

  # generate polygon intersection data
  sim_ply <- poly_intersect(test_bbox$xmin, test_bbox$xmax, test_bbox$ymin, test_bbox$ymax, test_ply)

  # ensure that the function produces what is expected
  expect_equal(comp_ply, sim_ply)

})



test_that("chopped_area() functionality", {

  # load expected output data
  comp_ply <- readRDS("data/comp_chopped_area.rds")

  # generate polygon intersection data
  sim_ply <- chopped_area(test_bbox$xmax, test_bbox, test_ply)

  # ensure that the function produces what is expected
  expect_equal(comp_ply, sim_ply)

})



test_that("target_function() functionality", {

  # load expected output data
  comp_ply <- readRDS("data/target_function.rds")

  # generate function output
  sim_ply <- target_function(fraction=0.4, test_bbox$xmax, test_bbox, test_ply, total_area=sf::st_area(test_ply))

  # ensure that the function produces what is expected
  expect_equal(comp_ply, sim_ply)

})






