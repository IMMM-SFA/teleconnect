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
#' @return table of crop types present and their frequency of occurrence
#' @importFrom sf st_crs st_transform
#' @importFrom raster crop projection mask unique freq
#' @importFrom tibble as_tibble
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
    freq(useNA = "no") %>%
    as_tibble() %>%
    rename(Group.1 = value, x = count)

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
#' @importFrom dplyr group_indices left_join filter rename case_when
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

  # First link uses the GCAM crop table and links the CDL table by similar crop names. 40 IDs get filled by this operation.
  gcam_csv %>%
    dplyr::select(item, GTAP_crop, GCAM_commodity, IFA2002_crop) %>%
    left_join(cdl_classes, by = c("GTAP_crop" = "CDL_Class")) %>%
    left_join(cdl_classes, by = c("IFA2002_crop" = "CDL_Class")) %>%
    mutate(ID = case_when(
      is.na(ID.x) & is.na(ID.y) ~ NA_integer_,
      is.na(ID.x) & !is.na(ID.y) ~ ID.y,
      !is.na(ID.x) & is.na(ID.y) ~ ID.x,
      !is.na(ID.x) & !is.na(ID.y) ~ ID.x)) %>%
    left_join(gcam_classes, by = "GCAM_commodity") %>%
    dplyr::select(ID, GCAM_ID, item) %>%
    left_join(cdl_classes, gcam_csv, by = "ID") -> gcam_link

  # This is then merged with the CDL IDs. All of the IDs that are still NA are made into a new table to do the second link.
  crop_cover_levels %>%
    left_join(gcam_link, by = "ID") %>%
    left_join(gcam_classes, by = "GCAM_ID") %>%
    dplyr::select(ID, Class_Names, GCAM_ID, CDL_Class, GCAM_commodity) -> crop_level_join

  crop_level_join[is.na(crop_level_join$GCAM_commodity), ] -> na_crops

  # This second link uses grepl to assign crops that are still NA to GCAM crops.
  na_crops %>%
    mutate(GCAM_Class = case_when(
      Class_Names %in% "Oats" ~ "OtherGrain",
      grepl("Grains", Class_Names) | grepl("Other Small Grains", Class_Names) ~ "OtherGrain",
      Class_Names %in% c("Corn", "Sweet Corn", "Pop or Orn Corn") ~ "Corn",
      Class_Names %in% c("Sod/Grass Seed", "Switchgrass") ~ "FodderGrass",
      Class_Names %in% c("Canola", "Rape Seed", "Camelina") ~ "OilCrop",
      Class_Names %in% c("Other Hay/Non Alfalfa", "Clover/Wildflowers", "Vetch") ~ "FodderHerb",
      grepl("Wheat", Class_Names) ~ "Wheat",
      grepl("Speltz", Class_Names) ~ "Wheat",
      grepl("Sweet Potatoes", Class_Names) ~ "Root_Tuber",
      grepl("Sugar", Class_Names) ~ "SugarCrop",
      grepl("Flaxseed", Class_Names) ~ "FiberCrop",
      grepl("Dbl Crop Barley", Class_Names) ~ "OtherGrain",
      grepl("Dbl Crop Corn", Class_Names) ~ "Corn",
      grepl("Wht", Class_Names) ~ "Wheat",
      grepl("Dbl Crop Lettuce/", Class_Names) ~ "MiscCrop",
      grepl("Dbl Crop Soybeans/", Class_Names) ~ "OilCrop",
      grepl("Developed", Class_Names) | grepl("Water", Class_Names) | grepl("Perennial", Class_Names) |
      grepl("Barren", Class_Names) | grepl("Forest", Class_Names) | grepl("land", Class_Names) | grepl("Background", Class_Names) |
      grepl("Clouds", Class_Names) | grepl("Undefined", Class_Names) | grepl("Aqua", Class_Names) ~ NA_character_,
      TRUE ~ "MiscCrop")) %>%
    dplyr::select(ID, Class_Names, GCAM_Class) -> na_crops_class

  # This joins both tables together so that every CDL crop is assigned a GCAM crop
  left_join(crop_level_join, na_crops_class, by = "ID") -> full_crop_table

  # FIlls in NA of GCAM_commodity
  full_crop_table %>% mutate(GCAM_commodity = if_else(is.na(GCAM_commodity), GCAM_Class, GCAM_commodity)) %>%
    dplyr::select(ID, GCAM_commodity) -> linked_crop

  # Brings in GCAM Commodity and GCAM IDs. Then adds true false column for land/crop covers.
  cdl_classes %>%
    left_join(linked_crop, by = "ID") %>%
    left_join(gcam_classes, by = "GCAM_commodity") %>%
    mutate(is_crop = if_else(is.na(GCAM_ID), FALSE, TRUE)) -> gcam_cdl_link

  # Rename so the column names match those within the code.
  rename(gcam_cdl_link, CDL_ID = ID, GCAM_Class = GCAM_commodity) -> reclass_table

}

#' get_land_category
#' @details Classify watershed condition based on percent development and cultivation.
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export
get_land_category <- function(percent_area){

  if(percent_area <= 1){
    watershed_condtion <- "Very Low"
  }else if(percent_area > 1 & percent_area <= 5){
    watershed_condtion <-"Low"
  }else if(percent_area > 5 & percent_area <= 15){
    watershed_condtion <-"Average"
  }else if(percent_area > 15 & percent_area <= 30){
    watershed_condtion <-"High"
  }else if(percent_area > 30){
    watershed_condtion <-"Very High"
  }

}


#' get_city_utility_mapping
#' @details Load plant code file to map cities to Utility IDs
#' @import vroom
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export
get_city_utility_mapping <- function(){

  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/EIA_Plant_2018.csv"),
        col_select = c('UTILITY_ID',
                       'City',
                       'State'),
        col_types = cols(UTILITY_ID = col_factor(),
                         City = col_character(),
                         State = col_character()),
        delim = ",",
        skip = 1) -> util_data_2018

  util_data_2018$city_state <- with(util_data_2018, paste(City, "|", State))
  city_to_utility <- util_data_2018[,c("UTILITY_ID", "city_state")]

}
