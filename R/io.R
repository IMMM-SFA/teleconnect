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
#' @importFrom raster raster
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
  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/city_to_intake_mapping.csv"),
        col_types = cols(city = col_character(),
                         state = col_character(),
                         city_uid = col_integer(),
                         intake = col_character(),
                         DVSN_ID = col_integer(),
                         city_state = col_character())
  )
}

#' get_crop_mapping
#'
#' Read internal data file that specifies the GCAM classification for certain crop types.
#' @import vroom
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export
get_crop_mapping <- function(){
  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/FAO_ag_items_PRODSTAT.csv"),
        col_types = cols(item = col_character(),
                         GTAP_crop = col_character(),
                         GCAM_commodity = col_character(),
                         CROSIT_crop = col_character(),
                         CROSIT_cropID = col_double(),
                         IFA2002_crop = col_character(),
                         IFA_commodity = col_character(),
                         MIRCA_crop_name = col_character(),
                         MIRCA_crop = col_character(),
                         MIRCA_Crop24and26 = col_character(),
                         MH_crop = col_character(),
                         MH2011_crop = col_character(),
                         MH2014_proxy = col_character(),
                         GTAP_use = col_character(),
                         LPJmL_crop = col_character(),
                         GEPIC_crop = col_character(),
                         Pegasus_crop = col_character(),
                         C3avg_include = col_double())
  )
}

#' get_ucs_power_plants
#'
#' Read the UCS database
#' @param ucs_file_path full path of UCS xlsx file within data_dir
#' @param method tool used to convert to spatial data; either 'sp' or 'sf'
#' @importFrom readxl read_xlsx
#' @importFrom sf st_as_sf st_transform
#' @importFrom sp CRS SpatialPointsDataFrame
#' @importFrom dplyr select left_join bind_rows
#' @importFrom stringr str_to_title
#' @import vroom
#' @author Sean Turner (sean.turner@pnnl.gov)
#' @export
get_ucs_power_plants <- function(ucs_file_path,
                                 method = "sp"){

  # read in EIA Plant and utility data (2010)
  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/EIA_Plant_2010.csv"),
        col_select = c('UTILITY_ID',
                         'PLANT_CODE'),
        col_types = cols(UTILITY_ID = col_integer(),
                         PLANT_CODE = col_integer()),
        delim = ",") -> utility_data_2010

  # read in EIA Plant and Balancing Authorities (2015)
  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/EIA_Plant_2015.csv"),
        col_select = c('Plant Code',
                       'BA_NAME'),
        col_types = cols(`Plant Code` = col_integer(),
                         BA_NAME = col_character()),
        skip = 1, delim = ",") %>% unique() -> utility_data_2015

  # read in Utility ID -> Balancing Authority mapping
  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/Electric_Retail_Service_Territories.csv"),
        col_select = c('UTILITY_ID',
                       'CNTRL_AREA'),
        col_types = cols(UTILITY_ID = col_integer(),
                         CNTRL_AREA = col_character()),
        delim = ",") -> ba_data

  # read UCS generator database
  read_xlsx(ucs_file_path,
            sheet = "MAIN DATA", skip = 4) %>%
    select(cooling = `Requires cooling?`,
           cooling_tech = `Cooling Technology`,
           PLANT_CODE = `Plant Code`,
           `Power Plant Type` = Fuel,
           lat = Latitude, lon = Longitude) %>%
    # aggregate to plant level by removing generator variables (e.g., nameplate) ...
    # and taking unique columns...
    unique() -> ucs

  # add utilities to UCS table...
  left_join(ucs, utility_data_2010, by = "PLANT_CODE") -> ucs_util
  # ... then use those utilities to map in Balancing Authorities
  left_join(ucs_util, ba_data, by = "UTILITY_ID") -> ucs_util_and_BA_partial

  # for matches not made above, attempt to join BA from EIA tables
  ucs_util_and_BA_partial %>% filter(is.na(CNTRL_AREA)) %>%
    left_join(utility_data_2015, by = c("PLANT_CODE" = "Plant Code")) %>%
    mutate(CNTRL_AREA = BA_NAME) %>% select(-BA_NAME) -> ucs_missing_BA_from_EIA

  #combine
  ucs_util_and_BA_partial %>% filter(!is.na(CNTRL_AREA)) %>%
    dplyr::bind_rows(ucs_missing_BA_from_EIA) %>%
    mutate(CNTRL_AREA = stringr::str_to_title(CNTRL_AREA)) -> ucs_with_utilities_and_BAs


  if (method == "sf") return(st_as_sf(ucs_with_utilities_and_BAs,
                                      coords = c("lon", "lat"),
                                      crs = CRS(proj4_string)))

  if (method == "sp") return(SpatialPointsDataFrame(data = ucs_with_utilities_and_BAs,
                                                    coords = ucs_with_utilities_and_BAs[c("lon", "lat")],
                                                    proj4string = CRS(proj4_string)))
}

#' Get raster frequency per class from polygon input areas
#'
#' Get the frequency per class of raster values represented in the input raster dataset
#' when restricted to the input watershed polygons for a target city.  This ignores NA.
#'
#' @param raster_object character. An object of class RasterLayer.
#' @param polygon character. A polygon to define spatial boundary of raster value counts (e.g. a given city's watersheds)
#' @return matrix of raster value and frequency in target polygons
#' @importFrom sf st_crs st_transform
#' @importFrom raster crop projection mask unique freq
#' @author Chris R. Vernon (chris.vernon@pnnl.gov)
#' @export
get_raster_val_classes <- function(raster_object, polygon) {

  # transform polygon to sf object if not already
  if(class(polygon)[[1]] != "sf") polygon <- st_as_sf(polygon)

  # get the coordinate system of the input raster
  r_crs <- st_crs(projection(raster_object))

  # read in shapefile and transform projection to the raster CRS
  polys <- polygon %>%
    st_transform(crs = r_crs)

  # calculate the frequency of unique land classes from the input raster that are in the target polygons
  n_lcs <- crop(raster_object, polys) %>%
    mask(polys) %>%
    unique() # freq(useNA = "no")

  return(n_lcs)
}


#' Mask raster to polygon
#'
#' @details masks a raster file against a chosen polygon.
#' @param raster_object character. An object of class RasterLayer.
#' @param polygon character. A polygon to define spatial boundary of raster value counts (e.g. a given city's watersheds)
#' @importFrom sf st_crs st_transform st_as_sf
#' @importFrom raster crop projection mask unique
#' @author Sean Turner (sean.turner@pnnl.gov)
#' @export
mask_raster_to_polygon <- function(raster_object, polygon) {

  # transform polygon to sf object if not already
  if(class(polygon)[[1]] != "sf") polygon <- st_as_sf(polygon)

  # get the coordinate system of the input raster
  r_crs <- st_crs(projection(raster_object))

  # read in shapefile and transform projection to the raster CRS
  polys <- polygon %>%
    st_transform(crs = r_crs)

  # crop and mask
  n_lcs <- crop(raster_object, polys) %>%
    mask(polys)

  return(n_lcs)
}

#' reclassify_raster
#' @param crop_cover_levels levels of the crop cover raster file.
#' @details Read internal data file that specifies the GCAM classification for certain crop types. Reclassify CDL based on GCAM.
#' @importFrom dplyr group_indices left_join filter rename
#' @importFrom car recode
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export
reclassify_raster <- function(crop_cover_levels){

  # Load in the GCAM classes CSV.
  gcam_csv <- get_crop_mapping()

  # Extract only the GCAM Class column
  gcam_csv %>%
    dplyr::select(GCAM_commodity) %>%
    filter(!is.na(GCAM_commodity)) %>%
    mutate(GCAM_ID = dplyr::group_indices(., GCAM_commodity)) %>%
    unique() ->
    gcam_classes

  # Obtain the attribute table for the raster.
  crop_cover_levels %>%
    filter(Class_Names != "") %>%
    dplyr::select(ID, CDL_Class = Class_Names) %>%
    mutate(CDL_Class = as.character(CDL_Class))->
    cdl_classes

  # Assign the GCAM IDs to matching CDL Class.
        # All Double Crops are under MiscCrop
        # CDL crops were assigned to where they fit best if there was no overlap with GCAM classes.
        # Fallow/Idle Cropland assigned to MiscCrop because it was not present in GCAM classes.
        # Any crop that was potentially used for oil is assigned to OilCrop. Example: Camelina
        # Cantelopes is not present in GCAM classes and is counted as MiscCrop because Watermelon is MiscCrop. Both are melons.
        # Any grain that was not wheat was assigned to OtherGrain
        # 1 = Corn; 2 = FiberCrop; 3 = FodderGrass; 4 = FodderHerb; 5 = MiscCrop; 6 = OilCrop;
        # 7 = OtherGrain; 8 = PalmFruit; 9 = Rice; 10 = Root_Tuber; 11 = SugarCrop; 12 = Wheat
  cdl_classes$GCAM_ID <- car::recode(cdl_classes$ID, "c(1,12,13) = 1;
                                                      c(2,32,36) = 2;
                                                      c(59,60) = 3;
                                                      c(37,58,224) = 4; c(10,11,14,42,44,47,48,49,50,51,52,53,54,55,56,57,61,66,67,68,69,70,71,72,74,75,
                                                      76,77,204,206,207,208,209,210,212,213,214,215,216,217,218,219,220,221,222,223,227,229,230,231,232,233,
                                                      234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,225,226) = 5;
                                                      c(5,6,26,31,33,34,35,38,211) = 6;
                                                      c(4,21,25,27,28,29,39,205) = 7;
                                                      3 = 9;
                                                      c(43,46) = 10;
                                                      c(41,45) = 11;
                                                      c(22,23,24,30) = 12;
                                                      c(82,121,122,123,124,63,64,65,81,83,87,88,92,
                                                      111,112,131,141,142,143,152,176,190,195,0) = NA", as.factor = FALSE)

  # Join tables by GCAM_ID and add column to seperate crops and land classes. Rename columns for simplicity.
  reclass_table <- dplyr::left_join(cdl_classes, gcam_classes, by = "GCAM_ID") %>%
    mutate(is_crop = if_else(is.na(GCAM_commodity), FALSE, TRUE))

  reclass_table_df <- dplyr::rename(reclass_table, CDL_ID = ID,
                                            GCAM_Class = GCAM_commodity)

}
