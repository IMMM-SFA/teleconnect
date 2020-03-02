#' count_watershed_teleconnections
#'
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param watersheds_file_path path of watersheds shapefile within data_dir
#' @param powerplants_file_path path of power plants data file
#' @param crop_file_path path of crop cover raster
#' @param dams_file_path path of National Inventory of Dams "NID" point file
#' @param irrigation_file_path path of edited demeter irrigation file.
#' @param nuld_file_path path of land use raster file.
#' @param cities a vector of cities to be included in the count. If omitted, all cities will be included.
#' @param poly_slices integer for how may parts to split the watersheds polygons into to enable faster zonal stats
#' @param n_cores integer for the number of machine cores used to run the polygon slicing function. 2 is default for users with 16GB of RAM.
#' @details counts teleconnections assoicated with water supply catchments associated with each city
#' @importFrom purrr map_dfr
#' @importFrom dplyr filter group_indices left_join
#' @importFrom tibble tibble
#' @importFrom sf as_Spatial
#' @importFrom foreign read.dbf
#' @export
count_watershed_teleconnections <- function(data_dir,
                                            watersheds_file_path = "water/CWM_v2_2/World_Watershed8.shp",
                                            powerplants_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
                                            crop_file_path = "land/2016_90m_cdls/cdl_lowres_usa.img",
                                            crop_attribute_path = "land/2016_90m_cdls/cdl_lowres_usa.img.vat.dbf",
                                            dams_file_path = "water/nabd_fish_barriers_2012/nabd_fish_barriers_2012.shp",
                                            irrigation_file_path = "land/usa_demeter.csv",
                                            nlud_file_path = "land/usa_nlud_LR.tif",
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
  import_raster(paste0(data_dir, crop_file_path)) -> cropcover_USA

  read.dbf(paste0(data_dir,crop_attribute_path)) -> crop_cover_levels

  # read reclassified crop table
  reclassify_raster(crop_cover_levels) -> crop_reclass_table

  # read NLUD economic raster
  import_raster( paste0(data_dir, nlud_file_path)) -> economic_USA

  # read NID point file and select only Flood Control Dams (C = Flood Control)
  import_shapefile(paste0(data_dir, dams_file_path)) %>%
    subset(grepl("C", Purposes)) %>%
    as_Spatial() -> flood_control_dams

  #read demeter irrigation file
  get_demeter_file(paste0(data_dir, irrigation_file_path)) -> usa_irrigation

  # map through all cities, computing teleconnections
  cities %>%
    map_dfr(function(city){
      filter(city_watershed_mapping, city_state == !!city) %>%
        .[["DVSN_ID"]] -> city_intake_ids

      # subset watersheds shapefile for target city
      watersheds %>%
        subset(DVSN_ID %in% city_intake_ids) -> watersheds_city

      # catch cases with only groundwater points (i.e., no watershed polygons)
      if(nrow(watersheds_city) == 0){
        done(city)
        return(
          tibble(city = !! city,
                 n_watersheds = 0,
                 n_hydro = 0,
                 n_thermal = 0,
                 n_floodcontroldams = 0,
                 wtrshd_impact = NA_character_,
                 n_cropcover = 0,
                 n_utilities = 0,
                 n_balancauth = 0,
                 n_cities = 0,
                 n_irrigatedcrops = 0,
                 n_economicsectors = 0)
        )
      }else{

        # TELECONNECTION - NUMBER OF WATERSHEDS
        tc_n_watersheds <- length(city_intake_ids)
        # NOTE: CURRENTLY COUNTS NESTED WATERSHEDS; ADDITIONAL...
        # ... ALGORITHM NEEDED TO AVOID DOUBLE COUNTING

        # TELECONNECTION - NUMBER OF CITIES USING SAME WATERSHEDS
        get_cities() %>% filter(DVSN_ID %in% city_intake_ids, city_state != !!city) %>%
          select(city_state) %>% unique() %>% nrow() -> tc_n_cities

        # subset power plants for target city watersheds
        power_plants_USA[watersheds_city, ] -> power_plants_city

        # TELECONNECTION - NUMBER OF HYDRO PLANTS
        power_plants_city %>%
          subset(`Power Plant Type` == "Hydropower") %>%
          .[["PLANT_CODE"]] %>% unique() %>%
          length() ->
          tc_n_hydroplants

        # TELECONNECTION - NUMBER OF THERMAL PLANTS
        power_plants_city %>% subset(cooling == "Yes") %>%
          .[["PLANT_CODE"]] %>% unique() %>%
          length() ->
          tc_n_thermalplants

        # TELECONNECTION - NUMBER OF UTILITIES
        power_plants_city %>%
          subset(cooling == "Yes" | `Power Plant Type` == "Hydro") %>%
          .[["UTILITY_ID"]] %>% unique() %>%
          length() -> tc_n_utility

        # TELECONNECTION - NUMBER OF BALANCING AUTHORITIES
        power_plants_city %>%
          subset(cooling == "Yes" | `Power Plant Type` == "Hydro") %>%
          .@data %>% dplyr::select(PLANT_CODE, CNTRL_AREA) %>% unique() %>%
          .[["CNTRL_AREA"]] -> tc_ba_na

        sum(is.na(tc_ba_na)) -> n_missing_ba

        if(n_missing_ba > 0) message(paste0("For ", city,", ", n_missing_ba,
                                            " plant(s) not assigned a Balancing Authority"))

        tc_ba_na %>% .[!is.na(.)] %>% unique() %>% length() -> tc_n_ba

        # TELECONNECTION - NUMBER OF CROP TYPES BASED ON GCAM CLASSES. NUMBER OF LAND COVERS.

        # get raster values of crops within the watershed.
        get_zonal_data(cropcover_USA, watersheds_city) -> cropcover_ids

        # filter reclass table by IDs that match raster IDs.
        crop_reclass_table %>%
          filter(CDL_ID %in% cropcover_ids$Group.1) ->
          crop_and_landcover_types

        # filter out where "is_crop" is true and only count crop types.
        crop_and_landcover_types %>%
          filter(is_crop == TRUE)%>%
          .[["GCAM_ID"]] %>% unique() %>%
          length() -> tc_n_cropcover

        # TELECONNECTION - NUMBER OF FLOOD CONTROL DAMS WITHIN WATERSHED.
        flood_control_dams[watersheds_city, ] %>%
          length() -> tc_fcdam

        # TELECONNECTION - CLASSIFY WATERSHED BASED ON % OF DEVELOPED/CULTIVATED AREA.
        # Remove NA and all categories that are not land cover/use(water/background).
        cropcover_ids %>% filter(!Group.1 %in% non_land_cdl_classes) -> all_land
        # New df with only crops and developement categories
        cropcover_ids %>% filter(!Group.1 %in% non_devcrop_class) -> dev_and_crop
        # Add cell count for all the land to get total land coverage.
        totcells <- sum(all_land$x)
        # Add cell count for all development and crop counts.
        totaldevcrop <- sum(dev_and_crop$x)
        # Find percent area for development and crops.
        percent_area <- 100*totaldevcrop/totcells
        # Assign to category based on percent area.
        get_land_category(percent_area) -> watershed_condition

        # TELECONNECTION - Count number of irrigated and rainfed crops in watershed.
        # Get irrigation data points within city's watersheds
        usa_irrigation[watersheds_city, ] %>%
          as_tibble() -> irrigation_city
        # count the number of crop types that are irrigated
        get_irrigation_count(irrigation_city) %>%
          filter(GCAM_Class %in% crop_and_landcover_types$GCAM_Class) -> irr_crops
          length(irr_crops$irr_count) -> tc_n_irrigated_crops

        # TELECONNECTION - Count # of economic sectors within watershed
        # get raster values of land use types within the watershed.
        get_zonal_data(economic_USA, watersheds_city) -> economic_ids

        # merge ids with id table to attach class names
        get_nlud_names(economic_ids) -> nlud_table
        # Filter out water and count the unique class variables within the watershed
        nlud_table %>%
          filter(Reclass != "Water") %>%
          .[["Reclass"]] %>% unique() %>%
          length() -> tc_n_economicsectors

        done(city)

        # output
        return(
          tibble(city = !! city,
                 n_watersheds = tc_n_watersheds,
                 n_hydro = tc_n_hydroplants,
                 n_thermal = tc_n_thermalplants,
                 n_floodcontroldams = tc_fcdam,
                 wtrshd_impact = watershed_condition,
                 n_cropcover = tc_n_cropcover,
                 n_utilities = tc_n_utility,
                 n_balancauth = tc_n_ba,
                 n_cities = tc_n_cities,
                 n_irrigatedcrops = tc_n_irrigated_crops,
                 n_economicsectors = tc_n_economicsectors)
        )
      }

    })
}

#' count_utility_teleconnections
#'
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param powerplants_file_path path of power plants data file
#' @param utility_file_path path of electric retail service areas file
#' @param citypoint_file_path path of city centroid point file
#' @param cities a vector of cities to be included in the count. If omitted, all cities will be included.
#' @details counts teleconnections for service areas associated with each city
#' @importFrom purrr map_dfr
#' @importFrom sf st_intersection st_as_sf st_agr
#' @importFrom tibble tibble as_tibble
#' @export
count_utility_teleconnections <- function(data_dir,
                                          powerplants_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
                                          utility_file_path = "energy/Electric_Retail_service_Territories/Electric_Retail_service_Territories.shp",
                                          citypoint_file_path = "water/CWM_v2_2/City_Centroid.shp",
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
  # Load city mapping file
  get_cities() %>%
    subset(city_state %in% cities) -> city_mapping

  # read shapefiles for utility areas
  import_shapefile(paste0(data_dir, utility_file_path),
                   method = "rgdal") %>% st_as_sf() -> utilities

  # read ucs plant data
  get_ucs_power_plants(paste0(data_dir, powerplants_file_path)) %>%
    as_tibble() -> power_plants_USA

  # Load city centroid file and merge with city mapping file
  import_shapefile(paste0(data_dir, citypoint_file_path)) %>% st_as_sf() %>%
    rename(.,"city_uid" = "City_ID") %>%
    left_join(., city_mapping, by = "city_uid") %>%
    select(c("city_state","geometry")) %>% unique() -> city_points

  # map through all cities, computing utility teleconnections
  cities %>%
    map_dfr(function(city){
      filter(city_points, city_state == !!city) -> utility_city

        # intersect city points and utility polygons to find the service areas city belongs to.
        sf::st_agr(utility_city) = "constant"
        sf::st_agr(utilities) = "constant"
        suppress_intersect(utility_city, utilities) -> city_in_utility

        # count number of utilities that serve the city
        length(city_in_utility$NAME)-> tc_n_utilities

        # subset power plants for target city utilities
        power_plants_USA %>%
          subset(UTILITY_ID %in% city_in_utility$ID) -> power_plants_utility

        # count number of water dependent plants
        power_plants_utility %>%
          subset(`Power.Plant.Type` == "Hydropower" | cooling == "Yes") %>%
          .[["PLANT_CODE"]] %>% unique() %>%
          length() -> tc_n_water_dependent

      # output
      return(
        tibble(city = !! city,
               n_utilites = tc_n_utilities,
               n_waterdependentplants = tc_n_water_dependent)
      )
    })
}
