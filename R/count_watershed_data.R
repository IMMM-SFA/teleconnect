#' count_watershed_data
#'
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param file_paths paths to data files
#' @param cities a vector of cities to be included in the count. If omitted, all cities will be included.
#' @details counts teleconnections associated with water supply catchments associated with each city
#' @importFrom purrr map_dfr map
#' @importFrom dplyr filter group_indices left_join if_else tribble
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
count_watershed_data <- function(data_dir,
                                 cities = NULL,
                                 file_paths = c(
                                   watersheds = "water/CWM_v2_2/World_Watershed8.shp",
                                   withdrawal = "water/CWM_v2_2/Snapped_Withdrawal_Points.shp",
                                   citypoint = "water/CWM_v2_2/City_Centroid.shp",
                                   powerplants = "water/UCS-EW3-Energy-Water-Database.xlsx",
                                   crop = "land/2016_90m_cdls/cdl_lowres_usa.img",
                                   crop_attributes = "land/2016_90m_cdls/cdl_lowres_usa.img.vat.dbf",
                                   irrigation = "land/usa_demeter.csv",
                                   nlud = "land/usa_nlud_LR.tif",
                                   hydro = "energy/EHA_Public_PlantFY2019_GIS_6/ORNL_EHAHydroPlant_PublicFY2019final.xlsx",
                                   transfers = "water/transfers/USIBTsHUC6_Dickson.shp",
                                   climate = "land/kop_climate_classes.tif"
                                   )){

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

  watersheds <- city_watershed_mapping[["DVSN_ID"]] %>%
    unique()

  # watershed_mapping %>%
  #   dplyr::select(city_state, DVSN_ID, city_uid, intake) %>%
  #   unique() -> water_mapping_select

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, file_paths["watersheds"]),
                   method = "rgdal") %>%
    # subset for desired watersheds
    subset(DVSN_ID %in% watersheds) -> catchment_shapes

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, file_paths["withdrawal"]),
                   method = "rgdal") %>% st_as_sf() %>%
    subset(DVSN_ID %in% watersheds) -> withdrawal_points

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, file_paths["citypoint"]),
                   method = "rgdal") %>% st_as_sf() %>%
    rename("city_uid" = "City_ID") %>%
    subset(city_uid %in% city_watershed_mapping$city_uid) -> city_points

  # read ucs plant data
  get_ucs_power_plants(paste0(data_dir, file_paths["powerplants"])) -> power_plants_USA

  # read croptype raster for US
  import_raster(paste0(data_dir, file_paths["crop"])) -> cropcover_USA

  read.dbf(paste0(data_dir, file_paths["crop_attributes"])) -> crop_cover_levels

  # read reclassified crop table
  reclassify_raster(crop_cover_levels) -> crop_reclass_table

  # read NLUD economic raster
  import_raster(paste0(data_dir, file_paths["nlud"])) -> economic_USA

  # read NID point file
  dams::nid_cleaned -> nid_dataset

  # read hydrosource hydropower dataset
  get_hydro_dataset(data_dir, file_paths["hydro"]) %>%
    st_as_sf(coords = c("lon", "lat"),
             crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
    as_Spatial() ->
    hydro_points

  #read demeter irrigation file
  get_demeter_file(paste0(data_dir, file_paths["irrigation"])) -> usa_irrigation

  sup(import_shapefile(paste0(data_dir, file_paths["transfers"]), method = "rgdal")) %>%
    st_as_sf() %>% st_transform("+proj=longlat +datum=WGS84 +no_defs") -> interbasin_transfers

  import_raster(paste0(data_dir, file_paths["climate"])) -> us_climate

  # map through all cities, computing teleconnections
  watersheds %>%
    map(function(watershed){

      # subset watersheds shapefile for target watershed
        sup(subset(catchment_shapes, DVSN_ID %in% watershed)) -> watersheds_select

      # deal with groundwater cases by returning blank result
      if(nrow(watersheds_select) == 0){
          done(watershed)
          return(list(counts = tibble(item = as.character(), counts = as.integer()),
                 metrics = tibble(),
                 land_table = tibble(),
                 climate_zones = vector()))
        }else{


          # get area of watershed
          sup(areaPolygon(watersheds_select)) -> polygon_area

          #-------------------------------------------------------
          # TELECONNECTION - NUMBER OF CITIES USING WATERSHED
          get_cities() %>%
            filter(DVSN_ID == watershed) %>%
            .[["city_state"]] %>% unique() %>% length() - 1 ->
            n_other_cities

          #-------------------------------------------------------
          # TELECONNECTION - HYDRO ENERGY
          # subset power plants for target watersheds
          sup(power_plants_USA[watersheds_select, ]) %>%
            as_tibble() ->
            watershed_power_plants

          # subset(watershed_power_plants, Power.Plant.Type == "Hydropower") ->
          #   hydro_plants

          # subset hydro plants for target watershed
          hydro_points[watersheds_select, ] %>% as_tibble() %>%
            select(PLANT_NAME, NID_ID, generation) ->
            hydro_plants
            # left_join(nid_dataset, by = "NID_ID") %>%
            # select(NID_ID, PLANT_NAME, generation,
            #        Max_Storage, Normal_Storage, NID_Storage) ->
            # hydro_plants

          nrow(hydro_plants) -> n_hydro_plants

          if_else(n_hydro_plants == 0, 0,
                  sum(hydro_plants$generation, na.rm = T)) ->
            hydro_gen_MWh

          #-------------------------------------------------------
          # TELECONNECTION - THERMAL ENERGY
          # subset power plants for target watersheds
          subset(watershed_power_plants, cooling == "Yes") %>%
            mutate(consumption = as.integer(consumption),
                   withdrawal = as.integer(withdrawal)) %>%
            select(PLANT_NAME, MWh, consumption, withdrawal) %>%
            group_by(PLANT_NAME) %>%
            summarise(MWh = sum(MWh, na.rm = T),
                      consum_BCM = sum(consumption, na.rm = T) * Mgallons_to_BCM,
                      withdr_BCM = sum(withdrawal, na.rm = T) * Mgallons_to_BCM) ->
            thermal_plants

          nrow(thermal_plants) -> n_thermal_plants

          # Generation
          if_else(nrow(thermal_plants) == 0, 0,
                  sum(thermal_plants$MWh)) ->
            thermal_gen_MWh

          # Consumption
          if_else(nrow(thermal_plants) == 0, 0,
                  sum(thermal_plants$consum_BCM)) ->
            thermal_consum_BCM

          # Withdrawal
          if_else(nrow(thermal_plants) == 0, 0,
                  sum(thermal_plants$withdr_BCM)) ->
            thermal_withdr_BCM

          #-------------------------------------------------------
          # TELECONNECTION - NUMBER OF CROP TYPES BASED ON GCAM CLASSES. NUMBER OF LAND COVERS.

          # get raster values of crops within the watershed.
          get_zonal_data(cropcover_USA, watersheds_select) ->
            cropcover_ids

          # filter reclass table by IDs that match raster IDs.
          crop_reclass_table %>%
            left_join(cropcover_ids, by = c("CDL_ID" = "Group.1")) %>%
            as_tibble() %>% rename(cell_freq = x) %>%
            arrange(CDL_ID) -> land_table

          if (sum(land_table$cell_freq, na.rm = T) != sum(cropcover_ids$x, na.rm = T)){
            warning("land table area is consistent with cropcover_ids area!!")
          }


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
          # usa_irrigation[watersheds_select, ] %>%
          #   as_tibble()  %>%
          # # count the number of crop types that are irrigated
          #   get_irrigation_count() %>%
          #   filter(GCAM_Class %in% crop_and_landcover_types$GCAM_Class) -> irr_crops
          # length(irr_crops$irr_count) -> tc_n_irrigated_crops

          #--------------------------------------------------------
          # TELECONNECTION - Count # of economic sectors within watershed
          # get raster values of land use types within the watershed.
          get_zonal_data(economic_USA, watersheds_select) %>%
            # merge ids with id table to attach class names
            get_nlud_names() -> nlud_table


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
          get_zonal_data(us_climate, watersheds_select) %>%
            filter(Group.1 != 128) %>% .[["Group.1"]] -> climate_ids

          return(
            list(
              counts = tribble(
                ~item,                           ~count,
                "other cities using watershed",  n_other_cities,
                "hydro plants",                  n_hydro_plants,
                "thermal plants",                n_thermal_plants,
               #"reservoirs",                    n_reservoirs,
                "transfers in",                  n_transfers_into,
                "transfers out",                 n_transfers_out,
                "transfers within",              n_transfers_within
              ),
              metrics = tribble(
                ~metric,                           ~unit,   ~value,
                "watershed area",                  "sq_km",  polygon_area,
                "total hydro generation",          "MWh",    hydro_gen_MWh,
                "total thermal generation",        "MWh",    thermal_gen_MWh,
                "total thermal water consumption", "BCM/yr", thermal_consum_BCM,
                "total thermal water withdrawals", "BCM/yr", thermal_withdr_BCM
              # "total reservoir storage",         "BCM",    xxx
                ),
              land = land_table,
              economic_sectors = nlud_table,
              #irrigation = ????
              climate_zones = climate_ids
              )
            )
        }

      }) -> watershed_output

  names(watershed_output) <- watersheds

  return(watershed_output)

}







#---------------------------------------------------------
# # TELECONNECTION - Distance between withdrawal point and city
#
# # Get cities that watershed is connected to
# city_watershed_mapping %>%
#   filter(DVSN_ID == watershed) -> watershed_cities
#
# watershed_cities[c("city_uid", "DVSN_ID","DVTR_ID")] -> city_connect
#
# # Subset city points by previous result, then add row.id to connect to distances
# city_points %>%
#   filter(city_uid %in% watershed_cities$city_uid)  -> citypoints_df
# citypoints_df$row.id <- as.numeric(rownames(citypoints_df))
#
# # Create into spatial object
# citypoints_df %>% as_Spatial() -> wtrsh_city_points
#
# # Subset withdrawal points by watershed DVSN_ID and create to spatial object
# withdrawal_points %>%
#   filter(DVSN_ID == watershed) %>% as_Spatial() -> wtrshd_withdrawal
#
# # Find distance bewtween withdrawal point and all cities associated
# distGeo(wtrshd_withdrawal, wtrsh_city_points) %>% as.data.frame() -> point_distance
# point_distance$row.id <- as.numeric(rownames(point_distance))
#
# # Find max distance for watershed
# max(point_distance$.) -> max_wtrshd_distance
#
# # Join up with cities
# left_join(citypoints_df, point_distance, by = "row.id") %>%
#   rename(Distance_to_City = ".") %>%
#   left_join(city_connect, by = "city_uid") %>%
#   mutate(DVSN_ID = as.integer(DVSN_ID)) %>% as.data.frame() %>%
#   dplyr::select(city_uid, Distance_to_City, DVSN_ID, DVTR_ID) ->
#   city_distance




