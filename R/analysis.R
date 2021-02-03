#' Analysis tools for teleconnect


#' Create an altered bounding box from an input sf polygon object
#'
#' @param ply An sf polygon object
#' @return A bounding box keyed list with altered coordinates
#' @importFrom sf st_bbox
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
init_bbox <- function(ply) {

  bb <- st_bbox(ply)
  delta <- (bb$ymax[[1]] - bb$ymin[[1]]) / 10
  xmin <- bb$xmin[[1]] - delta
  ymin <- bb$ymin[[1]] - delta
  ymax <- bb$ymax[[1]] + delta
  xmin_source <- bb$xmin[[1]]
  xmax_source <- bb$xmax[[1]]

  poly_box <- list()
  poly_box[[ "xmin" ]] <- xmin
  poly_box[[ "ymin" ]] <- ymin
  poly_box[[ "ymax" ]] <- ymax
  poly_box[[ "xmin_source" ]] <- xmin_source
  poly_box[[ "xmax_source" ]] <- xmax_source

  return(poly_box)
}


#' Create a slice of the original polygon object from a smaller bounding box
#'
#' @param xmin A float or integer value for the x (longitude) coordinate minimum
#' @param xmax A float or integer value for the x (longitude) coordinate maximum
#' @param ymin A float or integer value for the y (latitude) coordinate minimum
#' @param ymax A float or integer value for the y (latitude) coordinate maximum
#' @param ply An sf polygon object
#' @return An sf polygon object
#' @importFrom sf st_polygon st_as_sf st_set_crs st_crs st_intersection st_geometry
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
poly_intersect <- function(xmin, xmax, ymin, ymax, ply) {

  return (matrix(c(xmin, ymin, xmin, ymax, xmax, ymax, xmax, ymin, xmin, ymin), ncol = 2, byrow = TRUE) %>%
          list() %>%
          st_polygon() %>%
          as('Spatial') %>%
          st_as_sf() %>%
          st_set_crs(st_crs(ply)) %>%
          st_intersection(st_geometry(ply)))
}


#' Calculate the area of a slice of a polygon as numeric
#'
#' @param xmax A float or integer value for the x (longitude) coordinate maximum
#' @param coords A list containing xmin, ymin, and ymax
#' @param ply An sf polygon object
#' @return numeric area of a polygon slice
#' @importFrom sf st_area
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
chopped_area <- function(xmax, coords, ply) {

  # create polygon bounds list for slice
  poly_slice <- poly_intersect(coords$xmin, xmax, coords$ymin, coords$ymax, ply)

  # get slice area from the interescted source polygon
  part_area <- st_area(poly_slice) %>% as.numeric()

  return(part_area)
}


#' Calculate the area of the intersected polygon portion
#'
#' @param fraction The fraction of area being processed represented from 0.0 to 1.0
#' @param xmax A float or integer value for the x (longitude) coordinate maximum
#' @param coords A list containing xmin, ymin, and ymax
#' @param ply An sf polygon object
#' @param total_area total area of the whole polygon
#' @return the area of the polygon portion
#' @importFrom sf st_area
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
target_function <- function(fraction, xmax, coords, ply, total_area) {

  target = total_area * fraction

  alt_area <- chopped_area(xmax, coords, ply) - target

  return(alt_area)
}


#' Create next polygon slice in sequence
#'
#' @param ply An sf polygon object
#' @param xmax A float or integer value for the x (longitude) coordinate maximum
#' @param xmin A float or integer value for the x (longitude) coordinate minimum
#' @return new polygon slice that is next in sequence
#' @importFrom sf st_bbox
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
slicer <- function(ply, xmin, xmax){

  bb = st_bbox(ply)
  delta <- (bb$ymax[[1]] - bb$ymin[[1]]) / 10
  ymin <- bb$ymin[[1]] - delta
  ymax <- bb$ymax[[1]] + delta

  r = poly_intersect(xmin[[1]], xmax[[1]], ymin, ymax, ply)

  return(r)
}


#' Create a spatial polygon bounding box sf object
#'
#' Creates a spatial polygon bounding box from a user-provided extent
#' and coordinate reference system.
#'
#' @param x_min A float or integer value for the x (longitude) coordinate minimum
#' @param x_max A float or integer value for the x (longitude) coordinate maximum
#' @param y_min A float or integer value for the y (latitude) coordinate minimum
#' @param y_max A float or integer value for the y (latitude) coordinate maximum
#' @param my_crs An integer for the EPSG number of the desired output coordinate
#' reference system.
#' @return A bounding box polygon as an sf object
#' @importFrom sf st_polygon st_sfc
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
polygon_bounding_box <- function(x_min, x_max, y_min, y_max, my_crs) {

  bbox <- matrix(c(x_min, y_max,
                   x_max, y_max,
                   x_max, y_min,
                   x_min, y_min,
                   x_min, y_max), byrow = TRUE, ncol = 2) %>%
    list() %>%
    st_polygon() %>%
    st_sfc(., crs = my_crs)

  return(bbox)
}


#' Create a spatial polygon bounding box sf object
#'
#' Creates a spatial polygon bounding box from a user-provided extent
#' and coordinate reference system.
#'
#' @param x_min A float or integer value for the x (longitude) coordinate minimum
#' @param x_max A float or integer value for the x (longitude) coordinate maximum
#' @param y_min A float or integer value for the y (latitude) coordinate minimum
#' @param y_max A float or integer value for the y (latitude) coordinate maximum
#' @param my_crs An integer for the EPSG number of the desired output coordinate
#' reference system.
#' @return A bounding box polygon as an sf object
#' @importFrom sf st_polygon st_sfc
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
build_polygon <- function(xmin, xmax, ymin, ymax) {

  return(matrix(c(xmin, ymin, xmin, ymax, xmax, ymax, xmax, ymin, xmin, ymin), ncol = 2, byrow = TRUE))
}
