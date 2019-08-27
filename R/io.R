# Input/Output IO tools for teleconnect


#' Create an sf object from a shapefile
#'
#' Create an sf object from a full path to shapefile with file name and extension
#'
#' @param shp_path character. A full path to the input shapefile with file name and extension
#' @param quiet boolean.
#' @param method tool used to convert to spatial data; either 'sf' or 'rgdal'
#' @return sf object
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @importFrom sf st_read
#' @importFrom rgdal readOGR
#' @export
import_shapefile <- function(shp_path, quiet = TRUE, method = "sf") {

  if(method == "sf") return(st_read(shp_path, quiet = quiet))
  if(method == "rgdal") return(readOGR(shp_path, verbose = !quiet))
}


#' Create a raster object from a file path
#'
#' Create a raster object from a full path to raster file name and extension
#'
#' @param raster_path character. A full path to the input raster file with file name and extension
#' @return raster object
#' @importFrom fasterize fasterize
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
import_raster <- function(raster_path) {

  return(raster(raster_path))
}


#' Import NetCDF to brick raster
#'
#' Import NetCDF to brick raster
#'
#' @param ncdf_file character. A full path to the input NetCDF file with the file name and extension
#' @importFrom raster brick
#' @return raster object
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
import_ncdf_to_raster <- function(ncdf_file) {

  return(brick(ncdf_file))
}

#' Import point data from a CSV file
#'
#' Import point data that contains a value to be spatially joined to the fishnet containing
#' fractional area.  May either be a shapefile or a CSV file containing a latitude and longitude
#' for each record.
#'
#' @param f character. The full path with filename and extension to the points dataset.
#' @param pts_lat_field character. The field name for latitude
#' @param pts_lon_field character. The field name for longitude
#' @param pts_crs int. The native EPSG number for the coordinate reference system used in the
#' creation of the input points data. The default is 4326 (WGS 1984).
#' @param my_crs int. The EPSG number of the desired coordinate reference system. The
#' default is EPSG:3857 the WGS 84 / Pseudo-Mercator -- Used by all modern web
#' mapping applications.
#' @importFrom sf st_as_sf st_transform
#' @return A simple features (sf) spatial data frame object.
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
import_points_from_csv <- function(f, pts_lat_field, pts_lon_field, pts_crs = 4326, my_crs = 3857) {

  pts <- read.csv(file = f, header = TRUE, sep = ',')

  # change latitude, longitude columns to numeric
  cols.num <- c(pts_lat_field, pts_lon_field)
  pts[cols.num] <- sapply(pts[cols.num], as.numeric)

  # convert to sf spatial data frame object and transform to target CRS
  pts.SP <- st_as_sf(pts, coords = c(pts_lat_field, pts_lon_field), crs = pts_crs) %>%
    st_transform(crs = my_crs)

  return(pts)
}

#' get_cities
#'
#' Read internal data file that specifies all 236 US cities from UWB database.
#' Provides consistent mapping across cities, intakes and watersheds.
#' @import vroom
#' @author Sean Turner (sean.turner@pnnl.gov)
get_cities <- function(){
  vroom(paste0(system.file("extdata/", package = "teleconnect"),
               "city_to_intake_mapping.csv"),
        col_types = cols(city = col_character(),
                         state = col_character(),
                         city_uid = col_integer(),
                         intake = col_character(),
                         DVSN_ID = col_integer(),
                         city_state = col_character())
  )
}

#' get_ucs_power_plants
#'
#' Read the UCS database
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param ucs_file_path path of UCS xlsx file within data_dir
#' @param method tool used to convert to spatial data; either 'sp' or 'sf'
#' @importFrom readxl read_xlsx
#' @importFrom sf st_as_sf st_transform
#' @importFrom sp CRS SpatialPointsDataFrame
#' @importFrom dplyr select
#' @author Sean Turner (sean.turner@pnnl.gov)
#' @export
get_ucs_power_plants <- function(data_dir,
                                 ucs_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
                                 method = "sp"){
  read_xlsx(paste0(data_dir,
                   ucs_file_path),
            sheet = "MAIN DATA", skip = 4) %>%
    select(cooling = `Requires cooling?`,
           cooling_tech = `Cooling Technology`,
           `Power Plant Type` = Fuel,
           `Nameplate Capacity (MW)`,
           lat = Latitude, lon = Longitude) ->
    ucs_plants

  if (method == "sf") return(st_as_sf(ucs_plants,
                                      coords = c("lon", "lat"),
                                      crs = CRS(proj4_string)))

  if (method == "sp") return(SpatialPointsDataFrame(data = ucs_plants,
                                                    coords = ucs_plants[c("lon", "lat")],
                                                    proj4string = CRS(proj4_string)))
}

#' Get raster value count from polygon input areas
#'
#' Get the count of raster values represented in the input raster dataset
#' when restricted to the input watershed polygons for a target city.
#'
#' @param raster_file character. A full path to the input raster file with file name and extension
#' @param polygon_shpfile character. A full path to the input polygon shapefile with file name and extension
#' @param city_name character. The target city name in the input polygon shapefile.
#' @return count of unique raster values in target polygons
#' @importFrom sf st_crs st_transform
#' @importFrom raster crop projection mask unique
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
get_raster_nvals <- function(raster_file, polygon_shpfile, city_name) {

  r <- import_raster(raster_file)

  # get the coordinate system of the input raster
  r_crs <- st_crs(projection(r))

  # read in shapefile and transform projection to the raster CRS
  polys <- import_shapefile(polygon_shpfile) %>%
    subset(City_Name == city_name) %>%
    st_transform(crs = r_crs)

  # calculate the number of unique land classes from the input raster that are in the target polygons
  n_lcs <- crop(r, polys) %>%
    mask(polys) %>%
    unique() %>%
    length()

  return(n_lcs)
}
