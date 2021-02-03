# Conversion tools for gamut


#' Create a polygon sf object from a raster
#'
#' Create a polygon sf object from a raster, stack, or brick.
#'
#' @param raster A raster, stack, or brick object
#' @param na.rm boolean. If TRUE will polygonize only non-NA cells. Defualt is FALSE.
#' @importFrom spex qm_rasterToPolygons
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
raster_to_polygon <- function(raster, na.rm = FALSE) {

  return(qm_rasterToPolygons(raster, na.rm))
}
