#' count_watershed_data
#'
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param file_paths paths to data files
#' @param cities a vector of cities to be included in the count. If omitted, all cities will be included.
#' @details counts teleconnections associated with water supply catchments associated with each city
#' @importFrom purrr map_dfr map
#' @importFrom dplyr filter group_indices left_join if_else tribble group_by summarise arrange
#' @importFrom tibble tibble
#' @importFrom sf as_Spatial st_as_sf st_cast st_within st_make_valid
#' @importFrom foreign read.dbf
#' @importFrom exactextractr exact_extract
#' @importFrom geosphere areaPolygon distGeo
#' @importFrom tmaptools set_projection
#' @importFrom lwgeom st_startpoint st_endpoint
#' @importFrom reservoir yield
#' @importFrom raster intersect extent
#' @importFrom stringr str_remove
#' @import rgeos
#' @import rgdal
#' @import dams
#' @export
count_watershed_data <- function(data_dir,
                                 cities = NULL,
                                 run_all,
                                 file_paths = c(
                                   watersheds = "water/CWM_v2_2/World_Watershed8.shp",
                                   withdrawal = "water/CWM_v2_2/Snapped_Withdrawal_Points.shp",
                                   citypoint = "water/CWM_v2_2/City_Centroid.shp",
                                   powerplants = "water/UCS-EW3-Energy-Water-Database.xlsx",
                                   crop = "land/2016_90m_cdls/cdl_lowres_usa.img",
                                   crop_attributes = "land/2016_90m_cdls/cdl_lowres_usa.img.vat.dbf",
                                   irrigation = "land/Version2_USA_Demeter.csv",
                                   nlud = "land/usa_nlud_LR.tif",
                                   hydro = "energy/EHA_Public_PlantFY2019_GIS_6/ORNL_EHAHydroPlant_PublicFY2019final.xlsx",
                                   transfers = "water/transfers/USIBTsHUC6_Dickson.shp",
                                   climate = "land/kop_climate_classes.tif",
                                   HUC4 = "water/USA_HUC4/huc4_to_huc2.shp",
                                   population = "land/pden2010_block/pden2010_60m.tif",
                                   runoff = "water/Historical_Mean_Runoff/USA_Mean_Runoff.tif"
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

  message(paste0("Processing ", length(watersheds), " watershed(s). This may take several minutes..."))

  # watershed_mapping %>%
  #   dplyr::select(city_state, DVSN_ID, city_uid, intake) %>%
  #   unique() -> water_mapping_select

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, file_paths["watersheds"]),
                   method = "rgdal") %>%
    # subset for desired watersheds
    subset(DVSN_ID %in% watersheds) -> catchment_shapes

  # read ucs plant data
  get_ucs_power_plants(paste0(data_dir, file_paths["powerplants"])) -> power_plants_USA

  # read croptype raster for US
  sup(import_raster(paste0(data_dir, file_paths["crop"]))) -> cropcover_USA

  read.dbf(paste0(data_dir, file_paths["crop_attributes"])) -> crop_cover_levels

  # read reclassified crop table
  reclassify_raster(crop_cover_levels) -> crop_reclass_table

  # read NLUD economic raster
  import_raster(paste0(data_dir, file_paths["nlud"])) -> economic_USA

  # read NID point file
  # dams::nid_cleaned %>%
  #   as.data.frame() -> nid_dataset

  dams::nid_subset -> nid_dataset

  nid_dataset %>%
    filter(!is.na(longitude),
           !is.na(latitude)) %>%
    # filter(!is.na(nid_dataset$Longitude),
    #        !is.na(nid_dataset$Latitude))  %>%
    as_tibble() ->
    nid_no_NA

  nid_spatial <- SpatialPointsDataFrame(coords = nid_no_NA %>%
                                          select(longitude, latitude),
                                          #select(Longitude, Latitude),
                                        data = nid_no_NA,
                                        proj4string = CRS(proj4_string))

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
  # read in irrigation bcm file
  get_irrigation_bcm() -> irrigation_bcm

  # read in US climate raster
  import_raster(paste0(data_dir, file_paths["climate"])) -> us_climate

  # read in HUC4 shapes
  import_shapefile(paste0(data_dir, file_paths["HUC4"]),
                   method = "rgdal") %>% st_as_sf() %>% select(c("HUC4")) -> HUC4_connect
  stringr::str_sub(HUC4_connect$HUC4, 1, 2) -> HUC4_connect$HUC2
  HUC4_connect %>% as_Spatial() %>% sp::spTransform(CRS("+proj=longlat +datum=WGS84 +no_defs")) -> HUC2_sp

  # read in runoff raster
  import_raster(paste0(data_dir, file_paths["runoff"])) -> runoff_raster

  # read in population raster
  import_raster(paste0(data_dir, file_paths["population"])) -> population_raster

  # read in waste flow points
  get_wasteflow_points() -> wasteflow_points

  # read in watershed time series
  sup(get_watershed_ts()) -> flow

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


        # get area of watershed IN SQUARE KILOMETERS
        raster::area(watersheds_select) / m2_to_km2 -> polygon_area

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
        # TELECONNECTION - WATERSHED RUNOFF VALUES
        as.character(watershed) -> name
        flow[c("Monthly_Date",name)] -> select_ts

        colnames(select_ts)[2] <- "runoff_vals"

        select_ts %>%
          filter(runoff_vals == min(runoff_vals)) %>% .[["Monthly_Date"]] -> driest_month
        min(select_ts[c(2)]) * BCMmonth_to_m3sec -> driest_runoff_m3sec

        select_ts %>% separate(Monthly_Date, into = c("Year", "Month"), sep = "/") %>%
          group_by(Year) %>%
          dplyr::summarize_at(vars(runoff_vals), sum) -> select_ts_yearly

        (mean(select_ts_yearly$runoff_vals) * BCMyr_to_m3sec) -> historical_runoff_m3sec

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

        if(run_all == TRUE){

        #--------------------------------------------------------
        # TELECONNECTION - FIND RUNOFF VALUES FOR DEVELOPED AND CULTIVATED AREAS
        # Developed Runoff Calculation
        cropcover_USA %>%
          crop(watersheds_select) ->
          cropcover_agg_nonproj

        runoff_raster %>%
          crop(watersheds_select) ->
          runoff_agg

        cropcover_agg_crop <- crop(cropcover_agg_nonproj, runoff_agg)

        raster::extent(cropcover_agg_crop) <- raster::extent(runoff_agg)

        cropcover_agg <- resample(cropcover_agg_crop, runoff_agg)

        suppressMessages(get_runoff_values(cropcover_agg,
                                           runoff_agg,
                                           developed_values,
                                           polygon_area,
                                           land_table)) -> dev_runoff_m3persec

        # Crop Runoff Calculation
        crop_reclass_table %>% filter(!is.na(GCAM_Class)) %>%
          .[["CDL_ID"]] -> cultivated_values

        suppressMessages(get_runoff_values(cropcover_agg,
                                           runoff_agg,
                                           cultivated_values,
                                           polygon_area,
                                           land_table)) -> cultivated_runoff_m3persec

        # Rest of land runoff calculation

        append(developed_values, cultivated_values) -> crop_dev_vals

        crop_reclass_table %>%
          filter(!CDL_ID %in% crop_dev_vals) %>%
          #filter(!CDL_ID %in% non_land_cdl_classes) %>%
          .[["CDL_ID"]] -> other_vals

        suppressMessages(get_runoff_values(cropcover_agg,
                                           runoff_agg,
                                           other_vals,
                                           polygon_area,
                                           land_table)) -> other_runoff_m3persec

        # Total runoff calculation
        append(crop_dev_vals, other_vals) -> all_vals
        #get_runoff_values(cropcover_agg, runoff_agg, all_vals) -> total_runoff_m3persec

        # Calculate fractions
        dev_runoff_m3persec + cultivated_runoff_m3persec + other_runoff_m3persec -> expec_total_runoff

        }else{
          expec_total_runoff <- 0
          dev_runoff_m3persec <- 0
          cultivated_runoff_m3persec <- 0
        }

        #--------------------------------------------------------
        # TELECONNECTION - Count number of irrigated and rainfed crops in watershed.
        # Get irrigation data points within city's watersheds
        # usa_irrigation[watersheds_select, ] %>%
        #   as_tibble()  %>%
        # # count the number of crop types that are irrigated
        #   get_irrigation_count() %>%
        #   filter(GCAM_Class %in% crop_and_landcover_types$GCAM_Class) -> irr_crops
        # length(irr_crops$irr_count) -> tc_n_irrigated_crops

        # TELECONNECTION - Irrigation Consumption
        usa_irrigation[watersheds_select, ] %>%
          as_tibble()  %>%
          get_irrigation_count() -> irrigation_km2

        watersheds_select %>% st_as_sf() %>% sf::st_make_valid() -> wtrshd_sf
        HUC2_sp %>% st_as_sf() %>% sf::st_make_valid() -> HUC_sf

        sup(sf::st_intersection(wtrshd_sf, HUC_sf)) %>% as_tibble() -> watershed_HUC2

        if(watershed_HUC2[["HUC2"]] %>% unique() %>% length() == 1){
          HUC2_majority <- watershed_HUC2[["HUC2"]] %>% unique()
        }else{
          watershed_HUC2 %>%
            group_by(HUC2) %>% summarise(area = sum(AreaKM2)) %>%
            arrange(-area) %>% .[["HUC2"]] %>% .[1] -> HUC2_majority
        }

        HUC2_majority <- stringr::str_remove(HUC2_majority, "^0+")

        irrigation_bcm %>% filter(HUC2 %in% HUC2_majority) %>%
          select(GCAM_commodity, Irr_BCM_KM2) %>%
          right_join(irrigation_km2, by = "GCAM_commodity") %>%
          mutate(consumption_BCM = Irr_BCM_KM2 * irr) %>%
          .[["consumption_BCM"]] %>% sum(na.rm = T) ->
          total_irr_bcm

        #--------------------------------------------------------
        # TELECONNECTION - Count # of economic sectors within watershed
        # get raster values of land use types within the watershed.
        get_zonal_data(economic_USA, watersheds_select) %>%
          # merge ids with id table to attach class names
          get_nlud_names() -> nlud_table


        if(run_all == TRUE){
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

        }else{
          n_transfers_into <- 0
          n_transfers_out <- 0
          n_transfers_within <- 0
        }
        #---------------------------------------------------------
        # TELECONNTECTION - Number of climate zones watershed crosses
        get_zonal_data(us_climate, watersheds_select) %>%
          filter(Group.1 != 128) %>% .[["Group.1"]] -> climate_ids


        if(run_all == TRUE){
        #---------------------------------------------------------
        # TELECONNECTION - WATERSHED STORAGE
        sup(nid_spatial[watersheds_select, ]) %>%
          as_tibble() -> watershed_dams

        watershed_dams %>% nrow() -> n_dams

        sum(watershed_dams$normal_storage, na.rm = TRUE) * AF_to_BCM ->
          watershed_storage_BCM

        flow %>%
        pull(as.character(watershed)) -> flow_ts

        sup(yield(Q = flow_ts,
                  capacity = watershed_storage_BCM,
                  reliability = 1,
                  plot = FALSE,
                  double_cycle = TRUE)) %>% .[["Yield"]] * 12 ->
          watershed_yield_BCM
        }else{
          n_dams <- 0
          watershed_storage_BCM <- 0
          watershed_yield_BCM <- 0
        }
        #---------------------------------------------------------
        # TELECONNECTION - POPULATION WITHIN WATERSHED
        watersheds_select %>%
          sf::st_as_sf() %>%
          sf::st_transform(crs = crs(population_raster)) %>%
          sf::st_union() -> ply_union
        exactextractr::exact_extract(population_raster, ply_union) %>%
          .[[1]] %>% subset(coverage_fraction == 1) %>%
          .[["value"]] %>% mean() -> mean_pop_per_sqkm

        mean_pop_per_sqkm * polygon_area -> population_total

        population_total * avg_wateruse_ltr_per_day -> ltr_per_day
        #---------------------------------------------------------
        # Use waste flow points as a check
        wasteflow_points[watersheds_select, ] %>% as_tibble() %>%
          subset(PRES_FACILITY_OVERALL_TYPE == "Wastewater") %>%
          subset(DISCHARGE_METHOD == "Outfall To Surface Waters") %>%
          nrow() -> n_treatment_plants

        #---------------------------------------------------------

        return(
          list(
            counts = tribble(
              ~item,                           ~count,
              "other cities using watershed",  n_other_cities,
              "hydro plants",                  n_hydro_plants,
              "thermal plants",                n_thermal_plants,
              "dams",                          n_dams,
              "transfers in",                  n_transfers_into,
              "transfers out",                 n_transfers_out,
              "transfers within",              n_transfers_within,
              "outlow treatment plants",       n_treatment_plants
            ),
            metrics = tribble(
              ~metric,                           ~unit,     ~value,
              "watershed area",                  "sq_km",    polygon_area,
              "total hydro generation",          "MWh",      hydro_gen_MWh,
              "total thermal generation",        "MWh",      thermal_gen_MWh,
              "total thermal water consumption", "BCM/yr",   thermal_consum_BCM,
              "total thermal water withdrawals", "BCM/yr",   thermal_withdr_BCM,
              "total reservoir storage",         "BCM",      watershed_storage_BCM,
              "total watershed yield",           "BCM/yr",   watershed_yield_BCM,
              "total runoff from watershed",     "m3/s",     expec_total_runoff,
              "development runoff",              "m3/s",     dev_runoff_m3persec,
              "cultivated runoff",               "m3/s",     cultivated_runoff_m3persec,
              "total irrigation consumption",    "BCM",      total_irr_bcm,
              "watershed population",            "people",   population_total,
              "population water consumption",    "ltr/day",  ltr_per_day,
              "average historical runoff",       "m3/sec",   historical_runoff_m3sec,
              "driest_runoff_m3sec",             "m3/sec",   driest_runoff_m3sec,
            ),
            land = land_table,
            economic_sectors = nlud_table,
            climate_zones = climate_ids
          )
        )
      }
      done(watershed)
    }) -> watershed_output

  names(watershed_output) <- watersheds

  return(watershed_output)

}
