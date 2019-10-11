#' count_watershed_teleconnections
#'
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param watersheds_file_path path of watersheds shapefile within data_dir
#' @param powerplants_file_path path of power plants data file
#' @param crop_file_path path of crop cover raster
#' @param dams_file_path path of National Inventory of Dams "NID" point file
#' @param cities a vector of cities to be included in the count. If omitted, all cities will be included.
#' @details counts teleconnections assoicated with water supply catchments associated with each city
#' @importFrom purrr map_dfr
#' @importFrom dplyr filter group_indices left_join
#' @importFrom tibble tibble
#' @importFrom sf as_Spatial
#' @export
#'
count_watershed_teleconnections <- function(data_dir,
                                            watersheds_file_path = "water/CWM_v2_2/World_Watershed8.shp",
                                            powerplants_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
                                            crop_file_path = "land/2016_30m_cdls/2016_30m_cdls.img",
                                            dams_file_path = "water/nabd_fish_barriers_2012/nabd_fish_barriers_2012.shp",
                                            cities = NULL){

  all_cities <- get_cities()[["city_state"]]

  # use all cities if "cities" argument is omitted
  if(is.null(cities)) cities <- all_cities %>% unique()

  # throw error if any chosen city lies outside available cities
  if(any(!cities %in% all_cities)) {
    cities[!cities %in% all_cities] -> bad_cities
    stop(paste0(paste(bad_cities), ": not part of '
                teleconnect'!"))
  }

  get_cities() %>%
    subset(city_state %in% cities) %>%
    subset(key_watershed == TRUE) ->
    city_watershed_mapping

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, watersheds_file_path),
                   method = "rgdal") %>%
    # subset for desired watersheds
    subset(DVSN_ID %in% city_watershed_mapping$DVSN_ID) ->
    # remove City_Name; it's inconsistent with other shapefiles
    #dplyr::select(-City_Name) ->
    watersheds

  # read ucs plant data
  get_ucs_power_plants(paste0(data_dir, powerplants_file_path)) -> power_plants_USA

  # read croptype raster for US
  import_raster(paste0(data_dir, crop_file_path)) ->
    cropcover_USA

  # read reclassified crop table
  reclassify_raster(crop_cover_levels = levels(cropcover_USA)[[1]]) -> crop_reclass_table

  # read NID point file and select only Flood Control Dams (C = Flood Control)
  import_shapefile(paste0(data_dir, dams_file_path)) %>%
    subset(grepl("C", Purposes)) %>%
    sf::as_Spatial()-> flood_control_dams

  # map through all cities, computing teleconnections
  cities %>%
    map_dfr(function(city){
      filter(city_watershed_mapping, city_state == !!city) %>%
        .[["DVSN_ID"]] -> city_intake_ids

      # subset watersheds shapefile for target city
      watersheds %>%
        subset(DVSN_ID %in% city_intake_ids) ->
        watersheds_city

      # catch cases with only groundwater points (i.e., no watershed polygons)
      if(nrow(watersheds_city) == 0){
        done(city)
        return(
          tibble(city = !! city,
                 n_watersheds = 0,
                 n_hydro = 0,
                 n_thermal = 0,
                 n_landclasses = 0,
                 n_cropcover = 0,
                 n_floodcontroldams = 0,
                 n_utilities = 0,
                 n_balancauth = 0)
        )
      }else{

        # TELECONNECTION - NUMBER OF WATERSHEDS
        tc_n_watersheds <- length(city_intake_ids)
        # NOTE: CURRENTLY COUNTS NESTED WATERSHEDS; ADDITIONAL...
        # ... ALGORITHM NEEDED TO AVOID DOUBLE COUNTING

        # subset power plants for target city watersheds
        power_plants_USA[watersheds_city, ] -> power_plants_city

        # TELECONNECTION - NUMBER OF HYDRO PLANTS
        power_plants_city %>%
          subset(`Power Plant Type` == "Hydropower") %>%
          nrow() ->
          tc_n_hydroplants

        # TELECONNECTION - NUMBER OF THERMAL PLANTS
        power_plants_city %>% subset(cooling == "Yes") %>%
          nrow() ->
          tc_n_thermalplants

        # TELECONNECTION - NUMBER OF UTILITIES
        power_plants_city %>%
          .[["Utility_ID"]] %>% unique() %>%
          length() -> tc_n_utility

        # TELECONNECTION - NUMBER OF BALANCING AUTHORITIES
        power_plants_city %>%
          .[["Balancing_Authority"]] %>% unique() %>%
         length() -> tc_n_ba

        # TELECONNECTION - NUMBER OF CROP TYPES BASED ON GCAM CLASSES. NUMBER OF LAND COVERS.

        # get raster values of crops within the watershed.
        get_raster_val_classes(cropcover_USA, watersheds_city) -> cropcover_ids

            # filter reclass table by IDs that match raster IDs.
            crop_reclass_table %>%
              filter(CDL_ID %in% cropcover_ids) ->
              crop_and_landcover_types

            # filter out where "is_crop" is true and only count crop types.
            crop_and_landcover_types %>%
              filter(is_crop == TRUE)%>%
              .[["GCAM_ID"]] %>% unique() %>%
              length() -> tc_n_cropcover

            # filter out where "is crop" is false and count land classes, adding 1 class to account for agriculture.
            crop_and_landcover_types %>%
              filter(is_crop == FALSE) %>%
              .[["CDL_Class"]] %>% unique() %>%
              length() -> number_landclasses_ex_ag

            if (tc_n_cropcover > 0){
              tc_n_landclasses = number_landclasses_ex_ag + 1
            }else{
              tc_n_landclasses = number_landclasses_ex_ag
            }
        # TELECONNECTION - NUMBER OF FLOOD CONTROL DAMS WITHIN WATERSHED.
        flood_control_dams[watersheds_city, ] %>%
          length() -> tc_fcdam

        done(city)

        # output
        return(
          tibble(city = !! city,
                 n_watersheds = tc_n_watersheds,
                 n_hydro = tc_n_hydroplants,
                 n_thermal = tc_n_thermalplants,
                 n_landclasses = tc_n_landclasses,
                 n_cropcover = tc_n_cropcover,
                 n_floodcontroldams = tc_fcdam,
                 n_utilities = tc_n_utility,
                 n_balancauth = tc_n_ba)
        )
      }

    })
}
