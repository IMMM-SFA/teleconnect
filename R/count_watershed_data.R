#' count_watershed_data
#'
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param watersheds_file_path path of watersheds shapefile within data_dir
#' @param powerplants_file_path path of power plants data file
#' @param crop_file_path path of crop cover raster
#' @param dams_file_path path of National Inventory of Dams "NID" point file
#' @param irrigation_file_path path of edited demeter irrigation file.
#' @param nuld_file_path path of land use raster file.
#' @param watersheds a vector of cities to be included in the count. If omitted, all cities will be included.
#' @details counts teleconnections associated with water supply catchments associated with each city
#' @importFrom purrr map_dfr map
#' @importFrom dplyr filter group_indices left_join if_else
#' @importFrom tibble tibble
#' @importFrom sf as_Spatial st_as_sf st_cast st_within
#' @importFrom foreign read.dbf
#' @importFrom exactextractr exact_extract
#' @importFrom geosphere areaPolygon distGeo
#' @importFrom tmaptools set_projection
#' @importFrom lwgeom st_startpoint st_endpoint
#' @importFrom circular rad
#' @import dams
#' @export
count_watershed_data <- function(data_dir, watersheds = NULL){

  files <- teleconnect_files

  all_watersheds <- get_cities() %>%
    select(c("DVSN_ID", "key_watershed")) %>%
    filter(key_watershed == TRUE) %>%
    .[["DVSN_ID"]] %>% unique()

  # use all cities if "cities" argument is omitted
  if(is.null(watersheds)) watersheds <- all_watersheds

  # throw error if any chosen city lies outside available cities
  if(any(!watersheds %in% all_watersheds)) {
    watersheds[!watersheds %in% all_watersheds] -> bad_watersheds
    stop(paste0(paste(bad_watersheds), ": not part of '
                teleconnect'!"))
  }

  get_cities() %>%
    subset(DVSN_ID %in% watersheds) %>%
    subset(key_watershed == TRUE) ->
    watershed_mapping

  watershed_mapping[c(7,5,3,8)] %>%
    mutate(DVSN_ID = as.integer(DVSN_ID)) %>% unique() -> water_mapping_select

  watershed_mapping[c(5,4)] %>% mutate(DVSN_ID = as.integer(DVSN_ID)) %>% unique() -> intake_names

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, files["watersheds_file_path"]),
                   method = "rgdal") %>%
    # subset for desired watersheds
    subset(DVSN_ID %in% watershed_mapping$DVSN_ID) -> catchment_shapes

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, files["withdrawal_file_path"]),
                   method = "rgdal") %>% st_as_sf() %>%
    subset(DVSN_ID %in% watershed_mapping$DVSN_ID) -> withdrawal_points

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, files["citypoint_file_path"]),
                   method = "rgdal") %>% st_as_sf() %>%
    rename("city_uid" = "City_ID") %>%
    subset(city_uid %in% watershed_mapping$city_uid) -> city_points

  # read ucs plant data
  get_ucs_power_plants(paste0(data_dir, files["powerplants_file_path"])) -> power_plants_USA

  # read croptype raster for US
  import_raster(paste0(data_dir, files["crop_file_path"])) -> cropcover_USA

  read.dbf(paste0(data_dir, files["crop_attribute_path"])) -> crop_cover_levels

  # read reclassified crop table
  reclassify_raster(crop_cover_levels) -> crop_reclass_table

  # read NLUD economic raster
  import_raster( paste0(data_dir, files["nlud_file_path"])) -> economic_USA

  # read NID point file
  dams::nid_cleaned -> nid_dataset

  # read hydrosource hydropower dataset
  get_hydro_dataset(data_dir, files["hydro_file_path"]) %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = "+proj=longlat +datum=WGS84 +no_defs") %>% as_Spatial() -> hydro_points

  #read demeter irrigation file
  get_demeter_file(paste0(data_dir, files["irrigation_file_path"])) -> usa_irrigation

  sup(import_shapefile(paste0(data_dir, files["transfers_file_path"]), method = "rgdal")) %>%
    st_as_sf() %>% st_transform("+proj=longlat +datum=WGS84 +no_defs") -> interbasin_transfers

  import_raster(paste0(data_dir, files["climate_file_path"])) -> us_climate

  # map through all cities, computing teleconnections
  watersheds %>%
    map_dfr(function(watershed){

      # subset watersheds shapefile for target city
        sup(subset(catchment_shapes, DVSN_ID %in% watershed)) -> watersheds_select

        if(nrow(watersheds_select) == 0){
          done(watershed)
          return(
            tibble(DVSN_ID = watershed,
                 hydro_gen = 0,
                 hydro_stor = 0,
                 therm_consum = 0,
                 therm_gen = 0,
                 therm_withd = 0,
                 n_transfers_out = 0,
                 n_transfers_into = 0,
                 n_transfers_within = 0,
                 wtrshd_impact = NA_character_,
                 n_cropcover = 0,
                 crop_per_km2 = 0,
                 n_irrigatedcrops = 0,
                 n_economicsectors = 0,
                 n_climatezones = 0))

        }else{

          #-------------------------------------------------------
          # TELECONNECTION - NUMBER OF CITIES USING WATERSHED
          watershed_mapping %>%
          filter(DVSN_ID == watershed) %>%
          .[["city_state"]] %>% unique() %>% length() -> tc_n_cities

          #-------------------------------------------------------
          # TELECONNECTION - HYDRO ENERGY
          # subset power plants for target watersheds
          sup(power_plants_USA[watersheds_select, ] %>% as_tibble()) -> watershed_power_plants

          subset(watershed_power_plants, Power.Plant.Type == "Hydropower") -> hydro_plants

          # subset hydro dams for target watershed
          hydro_points[watersheds_select, ] %>% as.data.frame() %>%
            left_join(nid_dataset, by = "NID_ID") %>%
            select(hydro_variable) %>% na.omit() -> hydro_dams

          as.numeric(hydro_dams$Normal_Storage) -> hydro_dams_num

          # get area of watershed
          sup(areaPolygon(watersheds_select)) -> polygon_area

          if_else(nrow(hydro_plants) > 0, (sum(hydro_plants$MWh)) / polygon_area, 0) -> hydro_generation

          if_else(nrow(hydro_dams) > 0, (sum(hydro_dams_num[1])*acrefeet_conv) / polygon_area, 0) -> hydro_storage

          #-------------------------------------------------------
          # TELECONNECTION - THERMAL ENERGY
          # subset power plants for target watersheds
          subset(watershed_power_plants, cooling == "Yes") %>%
            mutate(consumption = as.integer(consumption),
                   withdrawal = as.integer(withdrawal)) -> thermal_plants

          # Generation per unit watershed area MWh/km2
          if_else(nrow(thermal_plants) > 0, (sum(thermal_plants$MWh)) / polygon_area, 0) -> thermal_generation
          # Consumption per unit area Mg/yr/km2
          if_else(nrow(thermal_plants) > 0, (sum(thermal_plants$consumption)) / polygon_area, 0) -> thermal_consumption
          # Withdrawal per unit area Mg/yr/km2
          if_else(nrow(thermal_plants) > 0, (sum(thermal_plants$withdrawal)) / polygon_area, 0) -> thermal_withdrawal

          #-------------------------------------------------------
          # TELECONNECTION - NUMBER OF CROP TYPES BASED ON GCAM CLASSES. NUMBER OF LAND COVERS.

          # get raster values of crops within the watershed.
          get_zonal_data(cropcover_USA, watersheds_select) -> cropcover_ids

          # filter reclass table by IDs that match raster IDs.
          crop_reclass_table %>%
            filter(CDL_ID %in% cropcover_ids$Group.1) ->
            crop_and_landcover_types

          # filter out where "is_crop" is true and only count crop types.
          crop_and_landcover_types %>%
            filter(is_crop == TRUE)%>%
            .[["GCAM_ID"]] %>% unique() %>%
            length() -> tc_n_cropcover

          #-------------------------------------------------------
          # TELECONNECTION - CROP AREA PER UNIT WATERSHED
          cropcover_ids %>% rename(CDL_ID = "Group.1") %>%
          left_join(crop_and_landcover_types, by = "CDL_ID") %>%
            filter(is_crop == TRUE) -> crops

          # Crop area per unit watershed km2 per km2
          sum(crops$x) / sum(cropcover_ids$x) -> crop_per_km2

          #-------------------------------------------------------
          # TELECONNECTION - CLASSIFY WATERSHED BASED ON % OF DEVELOPED/CULTIVATED AREA.
          # Remove NA and all categories that are not land cover/use(water/background).
          cropcover_ids %>% filter(!Group.1 %in% non_land_cdl_classes) -> all_land
          # New df with only crops and developement categories
          cropcover_ids %>% filter(!Group.1 %in% non_devcrop_class) -> dev_and_crop
          # Find percent area for development and crops.
          percent_area <- (100*(sum(dev_and_crop$x))) / (sum(all_land$x))
          # Assign to category based on percent area.
          get_land_category(percent_area) -> watershed_condition

          #--------------------------------------------------------
          # TELECONNECTION - Count number of irrigated and rainfed crops in watershed.
          # Get irrigation data points within city's watersheds
          usa_irrigation[watersheds_select, ] %>%
            as_tibble()  %>%
          # count the number of crop types that are irrigated
            get_irrigation_count() %>%
            filter(GCAM_Class %in% crop_and_landcover_types$GCAM_Class) -> irr_crops
          length(irr_crops$irr_count) -> tc_n_irrigated_crops

          #--------------------------------------------------------
          # TELECONNECTION - Count # of economic sectors within watershed
          # get raster values of land use types within the watershed.
          get_zonal_data(economic_USA, watersheds_select) %>%
            # merge ids with id table to attach class names
            get_nlud_names() -> nlud_table

          # Filter out water and count the unique class variables within the watershed
          nlud_table %>%
            filter(Reclass != "Water") %>%
            .[["Reclass"]] %>% unique() %>%
            length() -> tc_n_economicsectors

          # -------------------------------------------------------
          # TELECONNECTION - Interbasin Transfers in Watershed
          # subset interbasin transfers for the select watershed
          st_as_sf(watersheds_select) -> sf_watershed
          sup(interbasin_transfers[sf_watershed, ]) -> watershed_transfers

          # add a row.id and +1 so that it equals the the row number starting from 1
          watershed_transfers$row.id <- as.numeric(rownames(watershed_transfers))
          watershed_transfers$row.id <- watershed_transfers$row.id + 1

          # find all the transfers that are fully within the watershed shape, then join by row.id
          sup(st_within(interbasin_transfers, sf_watershed)) %>% as_tibble() %>%
            left_join(watershed_transfers, by = "row.id") -> transfers_within

          # subset the interbasin transfers by those that are with in and those that are outward
          subset(interbasin_transfers, OBJECTID %in% transfers_within$OBJECTID) -> inner_transfers

          subset(watershed_transfers, !(OBJECTID %in% inner_transfers$OBJECTID)) -> outer_transfers

          # get start and endpoints of the outer transfers to determine inflow/outflow
          get_startpoints(outer_transfers) -> transfer_start

          get_endpoints(outer_transfers) -> transfer_end

          # find transfers that start within the watershed and transfers that start outside of the watershed
          sup(transfer_start[sf_watershed, ]) %>%
            nrow() -> n_transfers_out
          sup(transfer_end[sf_watershed, ]) %>%
            nrow() -> n_transfers_into
          nrow(inner_transfers) -> n_transfers_within

          #---------------------------------------------------------
          # TELECONNTECTION - Number of climate zones watershed crosses
          get_zonal_data(us_climate, watersheds_select) -> climate_ids

          climate_ids[!(climate_ids$Group.1=="128"),] %>% nrow() -> n_climate_zones

          #---------------------------------------------------------
          # TELECONNECTION - Distance between withdrawal point and city
          # Get cities that watershed is connected to
          watershed_mapping %>%
            filter(DVSN_ID == watershed) -> watershed_cities

          watershed_cities[c("city_uid", "DVSN_ID","DVTR_ID")] -> city_connect

          # Subset city points by previous result, then add row.id to connect to distances
          city_points %>%
            filter(city_uid %in% watershed_cities$city_uid)  -> citypoints_df
          citypoints_df$row.id <- as.numeric(rownames(citypoints_df))

          # Create into spatial object
          citypoints_df %>% as_Spatial() -> wtrsh_city_points

          # Subset withdrawal points by watershed DVSN_ID and create to spatial object
          withdrawal_points %>%
            filter(DVSN_ID == watershed) %>% as_Spatial() -> wtrshd_withdrawal

          # Find distance bewtween withdrawal point and all cities associated
          distGeo(wtrshd_withdrawal, wtrsh_city_points) %>% as.data.frame() -> point_distance
          point_distance$row.id <- as.numeric(rownames(point_distance))

          # Find max distance for watershed
          max(point_distance$.) -> max_wtrshd_distance

          # Join up with cities
          left_join(citypoints_df, point_distance, by = "row.id") %>%
            rename(Distance_to_City = ".") %>%
            left_join(city_connect, by = "city_uid") %>%
            mutate(DVSN_ID = as.integer(DVSN_ID)) %>% as.data.frame() %>%
            select(2,6,7,8) -> city_distance


          done(watershed)
                return(
                tibble(DVSN_ID = as.integer(!! watershed),
                       hydro_gen_MWh_km2 = hydro_generation,
                       hydro_stor_Mg_yr_km2 = hydro_storage,
                       therm_consum_Mg_yr_km2 = thermal_consumption,
                       therm_gen_MWh_km2 = thermal_generation,
                       therm_withd_Mg_yr_km2 = thermal_withdrawal,
                       n_transfers_out,
                       n_transfers_into,
                       n_transfers_within,
                       wtrshd_impact = watershed_condition,
                       n_cropcover = tc_n_cropcover,
                       crop_per_km2,
                       n_irrigatedcrops = tc_n_irrigated_crops,
                       n_economicsectors = tc_n_economicsectors,
                       n_climatezones = n_climate_zones) %>%
                left_join(city_distance, by = "DVSN_ID"))

        }

      }) %>% left_join(water_mapping_select, ., by = "city_uid") %>%
                .[, !names(.) %in% c("DVSN_ID.x", "DVTR_ID")] %>%
                unique() %>% rename(DVSN_ID = "DVSN_ID.y") %>%
                mutate(DVSN_ID = as.integer(DVSN_ID)) %>%
                left_join(intake_names, by = "DVSN_ID") %>%
                .[c(1,2,3,20,4:19)]
  }

