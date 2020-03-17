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

  return(raster(raster_path, RAT = TRUE))
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
#' @export
get_cities <- function(){
  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/city_to_intake_mapping.csv"),
        col_types = cols(city = col_character(),
                         state = col_character(),
                         city_uid = col_integer(),
                         intake = col_character(),
                         DVSN_ID = col_integer(),
                         DVTR_ID = col_integer(),
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
           lat = Latitude, lon = Longitude,
           MWh = `Estimated Generation (MWh)`,
           PLANT_NAME = `Plant Name`,
           consumption = `Calculated Consumption (million gallons/yr)`,
           withdrawal = `Calculated Withdrawal (million gallons/yr)`) %>%
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
#' Get the coverage per class of raster values represented in the input raster dataset
#' when restricted to the input watershed polygons for a target city.
#'
#' @param raster_object character. An object of class RasterLayer.
#' @param polygon character. A polygon to define spatial boundary of raster value counts (e.g. a given city's watersheds)
#' @return table of crop types present and their coverage within the polygon area
#' @importFrom sf st_transform st_union st_as_sf
#' @importFrom raster crs
#' @importFrom exactextractr exact_extract
#' @importFrom dplyr rename n
#' @importFrom tmaptools aggregate_map
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export
get_zonal_data <- function(raster_object, polygon, city) {
  # extract raster projection
  raster_crs <- crs(raster_object)
  # transform projection and union all polygons into one
  polygon %>%
    st_as_sf() %>%
    st_transform(crs = raster_crs) %>%
    st_union() -> ply_union
  # use exactextractr to return classes and coverage fractions of each cell
  exactextractr::exact_extract(raster_object, ply_union) %>%
    .[[1]] %>% .[["value"]] %>%
    tabulate() %>% tibble(count = .) %>%
    mutate(value = 1:n()) %>% filter(count!=0) -> id_list

  rename(id_list, "Group.1" = "value",
                  "x" = "count") -> class_freq
  return(class_freq)
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
    filter(Class_Name != "") %>%
    dplyr::select(ID = Value, CDL_Class = Class_Name) %>%
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
    rename(ID = Value) %>%
    left_join(gcam_link, by = "ID") %>%
    left_join(gcam_classes, by = "GCAM_ID") %>%
    dplyr::select(ID, Class_Name, GCAM_ID, CDL_Class, GCAM_commodity) -> crop_level_join

  crop_level_join[is.na(crop_level_join$GCAM_commodity), ] -> na_crops

  # This second link uses grepl to assign crops that are still NA to GCAM crops.
  na_crops %>%
    mutate(GCAM_Class = case_when(
      Class_Name %in% "Oats" ~ "OtherGrain",
      grepl("Grains", Class_Name) | grepl("Other Small Grains", Class_Name) ~ "OtherGrain",
      Class_Name %in% c("Corn", "Sweet Corn", "Pop or Orn Corn") ~ "Corn",
      Class_Name %in% c("Sod/Grass Seed", "Switchgrass","Other Hay/Non Alfalfa") ~ "FodderGrass",
      Class_Name %in% c("Canola", "Rape Seed", "Camelina") ~ "OilCrop",
      Class_Name %in% c("Clover/Wildflowers", "Vetch") ~ "FodderHerb",
      grepl("Wheat", Class_Name) ~ "Wheat",
      grepl("Speltz", Class_Name) ~ "Wheat",
      grepl("Sweet Potatoes", Class_Name) ~ "Root_Tuber",
      grepl("Sugar", Class_Name) ~ "SugarCrop",
      grepl("Flaxseed", Class_Name) ~ "FiberCrop",
      grepl("Dbl Crop Barley", Class_Name) ~ "OtherGrain",
      grepl("Dbl Crop Corn", Class_Name) ~ "Corn",
      grepl("Wht", Class_Name) ~ "Wheat",
      grepl("Dbl Crop Lettuce/", Class_Name) ~ "MiscCrop",
      grepl("Dbl Crop Soybeans/", Class_Name) ~ "OilCrop",
      grepl("Developed", Class_Name) | grepl("Water", Class_Name) | grepl("Perennial", Class_Name) |
      grepl("Barren", Class_Name) | grepl("Forest", Class_Name) | grepl("land", Class_Name) | grepl("Background", Class_Name) |
      grepl("Clouds", Class_Name) | grepl("Undefined", Class_Name) | grepl("Aqua", Class_Name) ~ NA_character_,
      TRUE ~ "MiscCrop")) %>%
    dplyr::select(ID, Class_Name, GCAM_Class) -> na_crops_class

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
#' @param percent_area percent area of developed and cultivated areas
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

#' get_demeter_file
#' @param irrigation_file_path FUll path to the demeter data file
#' @details Load demeter file for irrigation and rainfed crop data
#' @import vroom
#' @importFrom sp SpatialPointsDataFrame
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export

get_demeter_file <- function(irrigation_file_path){

  vroom(irrigation_file_path,
        col_types = cols(
          latitude = col_double(),
          longitude = col_double(),
          corn_irr = col_double(),
          fibercrop_irr = col_double(),
          foddergrass_irr = col_double(),
          fodderherb_irr = col_double(),
          misccrop_irr = col_double(),
          oilcrop_irr = col_double(),
          othergrain_irr = col_double(),
          palmfruit_irr = col_double(),
          rice_irr = col_double(),
          root_tuber_irr = col_double(),
          sugarcrop_irr = col_double(),
          wheat_irr = col_double(),
          corn_rfd = col_double(),
          fibercrop_rfd = col_double(),
          foddergrass_rfd = col_double(),
          fodderherb_rfd = col_double(),
          misccrop_rfd = col_double(),
          oilcrop_rfd = col_double(),
          othergrain_rfd = col_double(),
          palmfruit_rfd = col_double(),
          rice_rfd = col_double(),
          root_tuber_rfd = col_double(),
          sugarcrop_rfd = col_double(),
          wheat_rfd = col_double(),
          area_sqkm = col_double())) -> demeter
  # convert to spatial points so that it can be masked by watershed.
  SpatialPointsDataFrame(coords = demeter[ ,c(2,1)], demeter, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")) -> usa_irrigation

  return(usa_irrigation)
}

#' get_irrigation_count
#' @param irrigation_city irrigation data subsetted by city watersheds
#' @details Count irrigated crops vs rainfed crops and determine irrigation status
#' @importFrom tibble enframe
#' @importFrom dplyr mutate case_when filter if_else
#' @importFrom tidyr separate spread
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export

get_irrigation_count <- function(irrigation_city){

  # add values of all columns to get sum (in meters) of area irrigated and rainfed.
    irrigation_city[,c(3:26)] * irrigation_city$area_sqkm -> crop_areas
    colSums(crop_areas) %>% enframe() -> dem_sums

  # create column for irrigation status
  dem_sums %>%
    mutate(name = if_else(grepl("root_tuber", name), gsub("t_t", "tt", name), name)) %>%
    separate(name, into = c("crop", "water"), sep = "_") %>%
    spread(water, value) %>%
    mutate(irr_count = if_else(irr > rfd, TRUE, FALSE)) -> demeter_sep

  sum(demeter_sep$irr) -> irrigation_km2

  # create column that will allow matching with the crop cover raster
  demeter_sep %>% mutate(GCAM_commodity = dplyr::case_when(
    grepl("corn", crop) ~ "Corn",
    grepl("fibercrop", crop) ~ "FiberCrop",
    grepl("foddergrass",crop) ~ "FodderGrass",
    grepl("fodderherb", crop) ~ "FodderHerb",
    grepl("misccrop", crop) ~ "MiscCrop",
    grepl("oilcrop", crop) ~ "OilCrop",
    grepl("othergrain", crop) ~ "OtherGrain",
    grepl("palmfruit", crop) ~ "PalmFruit",
    grepl("rice", crop) ~ "Rice",
    grepl("roottuber", crop) ~ "Root_Tuber",
    grepl("sugarcrop", crop) ~ "SugarCrop",
    grepl("wheat", crop) ~ "Wheat")) %>% select(c("GCAM_commodity","irr")) -> irrigation_table
  #  filter(demeter_sep$irr_count != FALSE) -> irrigated_crops

  return(irrigation_table)
}

#' get_nlud_names
#' @param economic_ids dataframe of nlud sector ids
#' @details Load in nlud_id_table.csv and merge with the economic_ids to get class names
#' @import vroom
#' @importFrom dplyr left_join
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export

get_nlud_names <- function(economic_ids){
  # Load id table
  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/nlud_id_table.csv"),
        col_select = c('Group.1',
                       'NLUD_Class',
                       'Reclass'),
        col_types = cols(Group.1 = col_integer(),
                         NLUD_Class = col_character(),
                         Reclass = col_character()),
        delim = ",", skip = 1) -> nlud_id_table
  # Join ids and class names
  left_join(economic_ids, nlud_id_table, by = "Group.1") -> nlud_table

  return(nlud_table)
}


#' get_hydro_dataset
#' @param economic_ids dataframe of nlud sector ids
#' @details Load in nlud_id_table.csv and merge with the economic_ids to get class names
#' @importFrom readxl read_xlsx
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export

get_hydro_dataset <- function(data_dir, hydro_file_path){
  # Load id table
  read_xlsx(paste0(data_dir, hydro_file_path),
            sheet = "Operational") %>%
    select(NID_ID = NID_ID,
           lat = Lat, lon = Lon,
           PLANT_NAME = PtName,
           generation = CH_MWh) -> plants

  return(plants)
}


#' get_population
#' @param city_state name of city_state
#' @details Census populations
#' @importFrom vroom vroom cols
#' @author Sean Turner (sean.turner@pnnl.gov)
get_population <- function(city_state){
  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/city_population.csv"),
        skip = 1, col_types = cols()) %>%
    filter(city_state == !! city_state) %>%
    .[["Population"]]
}

#' get_watershed_ts
#' @param watershed select watershed
#' @details Load in runoff time series and select the watershed being analyzed and its time series
#' @importFrom vroom vroom cols
#' @importFrom dplyr pull
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export
get_watershed_ts <- function(watershed){

  vroom(paste0(system.file("extdata", package = "teleconnect"),
                                "/teleconnect_runoff_bcm.csv"),
                         delim = ",", skip = 2, col_types = cols()) %>%
    pull(as.character(watershed))
}

#' get_irrigation_bcm
#' @param irrbcm_file_path irrigation bcm file path
#' @details load in irrigation bcm file path
#' @importFrom vroom vroom cols
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
get_irrigation_bcm <- function(){
  vroom(paste0(system.file("extdata", package = "teleconnect"),
               "/HUC2_Irrigation_Data.csv"),
        skip = 2, col_types = cols())

}
