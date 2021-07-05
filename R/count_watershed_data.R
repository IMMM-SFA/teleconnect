#' count_watershed_data
#'
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param cities a vector of cities to be included in the count. If omitted, all cities will be included.
#' @param file_paths file paths to all geospatial input datasets
#' @param run_all to be depreciated.  Runs current configuration.
#' @details counts teleconnections associated with water supply catchments associated with each city
#' @importFrom purrr map_dfr map
#' @importFrom dplyr filter group_indices left_join if_else tribble group_by summarise arrange right_join
#' @importFrom tibble tibble
#' @importFrom foreign read.dbf
#' @importFrom exactextractr exact_extract
#' @importFrom geosphere areaPolygon distGeo
#' @importFrom lwgeom st_startpoint st_endpoint
#' @importFrom sf as_Spatial st_as_sf st_cast st_within st_make_valid
#' @importFrom reservoir yield
#' @importFrom raster intersect extent
#' @importFrom stringr str_remove
#' @importFrom tidyr separate
#' @import rgeos
#' @import rgdal
#' @import dams
#' @export
count_watershed_data <- function(data_dir,
                                 cities = NULL,
                                 run_all = TRUE,
                                 file_paths = c(
                                   watersheds = "water/CWM_v2_2/World_Watershed8.shp",
                                   withdrawal = "water/CWM_v2_2/Snapped_Withdrawal_Points.shp",
                                   citypoint = "water/CWM_v2_2/City_Centroid.shp",
                                   powerplants = "water/UCS-EW3-Energy-Water-Database.xlsx",
                                   crop = "land/2016_90m_cdls/cdl_lowres_usa.img",
                                   crop_attributes = "land/2016_90m_cdls/cdl_lowres_usa.img.vat.dbf",
                                   irrigation = "land/Version2_USA_Demeter.csv",
                                   nlud = "land/usa_nlud_LR.tif",
                                   hydro = "energy/ORNL_EHAHydroPlant_FY2020revised.xlsx",
                                   transfers = "water/transfers/USIBTsHUC6_Dickson.shp",
                                   climate = "land/kop_climate_classes.tif",
                                   HUC4 = "water/USA_HUC4/huc4_to_huc2.shp",
                                   population = "land/pden2010_60m.tif",
                                   runoff = "water/UWSCatCH/USA_Mean_Runoff.tif",
                                   nhd_flow = "water/UWSCatCH/UWSCatCH_Intake_Flows.shp",
                                   contributions = "water/UWSCatCH/Watershed_Contributions.csv"
                                 )){

  all_cities <- get_cities()[["city_state"]]

  # use all cities if "cities" argument is omitted
  if(is.null(cities)) cities <- all_cities %>% unique()

  # throw error if any chosen city lies outside available cities
  if(any(!cities %in% all_cities)) {
    cities[!cities %in% all_cities] -> bad_cities
    stop(paste0(paste(bad_cities), ": not part of '
                gamut'!"))
  }

  get_cities() %>%
    subset(city_state %in% cities) %>%
    subset(key_watershed == TRUE) ->
    city_watershed_mapping

  watersheds <- city_watershed_mapping[["DVSN_ID"]] %>%
    unique()

  message(paste0("Processing ", length(watersheds), " watershed(s). This may take several minutes..."))

  # read shapefiles for watersheds
  import_shapefile(file.path(data_dir, file_paths["watersheds"]),
                   method = "rgdal") %>%
    # subset for desired watersheds
    subset(DVSN_ID %in% watersheds) -> catchment_shapes

  # read ucs plant data
  sup(get_ucs_power_plants(file.path(data_dir, file_paths["powerplants"]))) -> power_plants_USA

  # read croptype raster for US
  sup(import_raster(file.path(data_dir, file_paths["crop"]))) -> cropcover_USA

  read.dbf(file.path(data_dir, file_paths["crop_attributes"])) -> crop_cover_levels

  # read reclassified crop table
  reclassify_raster(crop_cover_levels) -> crop_reclass_table

  # read NLUD economic raster
  import_raster(file.path(data_dir, file_paths["nlud"])) -> economic_USA

  # read NID point file
  dams::nid_subset -> nid_dataset

  nid_dataset %>%
    filter(!is.na(longitude),
           !is.na(latitude)) %>%
    as_tibble() ->
    nid_no_NA

  nid_spatial <- SpatialPointsDataFrame(coords = nid_no_NA %>%
                                        select(longitude, latitude),
                                        data = nid_no_NA,
                                        proj4string = CRS(proj4_string))

  # read hydrosource hydropower dataset
  get_hydro_dataset(data_dir, file_paths["hydro"]) %>%
    st_as_sf(coords = c("lon", "lat"),
             crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
    as_Spatial() ->
    hydro_points

  #read demeter irrigation file
  get_demeter_file(file.path(data_dir, file_paths["irrigation"])) -> usa_irrigation

  sup(import_shapefile(file.path(data_dir, file_paths["transfers"]), method = "rgdal")) %>%
    st_as_sf() %>% st_transform("+proj=longlat +datum=WGS84 +no_defs") -> interbasin_transfers
  # read in irrigation bcm file
  get_irrigation_bcm() -> irrigation_bcm

  # read in US climate raster
  import_raster(file.path(data_dir, file_paths["climate"])) -> us_climate

  # read in HUC4 shapes
  import_shapefile(file.path(data_dir, file_paths["HUC4"]),
                   method = "rgdal") %>% st_as_sf() %>% select(c("HUC4")) -> HUC4_connect
  stringr::str_sub(HUC4_connect$HUC4, 1, 2) -> HUC4_connect$HUC2
  HUC4_connect %>% as_Spatial() %>% sp::spTransform(CRS("+proj=longlat +datum=WGS84 +no_defs")) -> HUC2_sp

  # read in runoff raster
  import_raster(file.path(data_dir, file_paths["runoff"])) -> runoff_raster

  # read in population raster
  import_raster(file.path(data_dir, file_paths["population"])) -> population_raster

  # read in waste flow points
  get_wasteflow_points() -> wasteflow_points

  # read in EPA facilities
  get_epa_facilities() -> epa_facilities

  # read in watershed time series
  get_watershed_ts(watersheds) -> runoff_totals

  # read in NHD flow shapefile
  import_shapefile(file.path(data_dir, file_paths["nhd_flow"])) %>% st_as_sf() -> watershed_nhd_flows

  # temporary fix for Jackson MS
  if(watersheds == 3743){
    runoff_totals <- get_watershed_ts(3683) %>% mutate(watershed = 3743)
  }

  # map through all cities, computing teleconnections
  watersheds %>%
    map(function(watershed){

      message(watershed)

      # subset watersheds shapefile for target watershed
      sup(subset(catchment_shapes, DVSN_ID %in% watershed)) -> watersheds_select

      # deal with groundwater cases by returning blank result
      if(nrow(watersheds_select) == 0){
        done(paste0(watershed, " (groundwater)"))
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

        # subset hydro plants for target watershed
        hydro_points[watersheds_select, ] %>% as_tibble() %>%
          select(PLANT_NAME, NID_ID, generation) ->
          hydro_plants

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
        # TELECONNECTION - WATERSHED RUNOFF AND FLOW VALUES
        runoff_totals %>%
          filter(watershed == !!watershed) ->
          watershed_runoff

        watershed_runoff[["flow_BCM"]] %>% mean() * BCMmonth_to_m3sec ->
          historical_runoff_mean_m3sec

        watershed_runoff %>%
          group_by(month) %>% summarise(flow = mean(flow_BCM)) %>%
          .[["flow"]] %>% min() * BCMmonth_to_m3sec ->
          historical_runoff_dry_m3sec

        # Calculate flow from NHD flowline

        if(watershed %in% watershed_nhd_flows$DVSN_ID){

          watershed_nhd_flows %>%
            filter(DVSN_ID %in% watershed) %>%
            .[[nhdplus_flow_metric]] -> nhd_flow_m3sec

          if_else(is.na(nhd_flow_m3sec),
                  historical_runoff_mean_m3sec,
                  nhd_flow_m3sec) ->
            historical_flow_mean_m3sec

        }else{

          historical_runoff_mean_m3sec -> historical_flow_mean_m3sec

        }

        # not yet accounting for within year flow variation
        historical_flow_dry_m3sec <- NA_real_

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

        if(run_all == TRUE){

        #--------------------------------------------------------
        # TELECONNECTION - FIND RUNOFF VALUES FOR DEVELOPED AND CULTIVATED AREAS
        # Developed Runoff Calculation
        cropcover_USA %>%
          mask_raster_to_polygon(watersheds_select) ->
          cropcover_agg_nonproj

        runoff_raster %>%
          mask_raster_to_polygon(watersheds_select) ->
          runoff_agg

        cropcover_agg_crop <- crop(cropcover_agg_nonproj, runoff_agg)

        raster::extent(cropcover_agg_crop) <- raster::extent(runoff_agg)

        if(identical(dim(cropcover_agg_crop), dim(runoff_agg)) == FALSE){
          cropcover_agg <- resample(cropcover_agg_crop, runoff_agg)
        }else{
          cropcover_agg_crop -> cropcover_agg
        }

        sup(get_runoff_values(cropcover_agg,
                              runoff_agg,
                              developed_values,
                              polygon_area,
                              land_table)) -> dev_runoff_m3persec

        # Crop Runoff Calculation
        crop_reclass_table %>% filter(!is.na(GCAM_Class)) %>%
          .[["CDL_ID"]] -> cultivated_values

        sup(get_runoff_values(cropcover_agg,
                              runoff_agg,
                              cultivated_values,
                              polygon_area,
                              land_table)) -> cultivated_runoff_m3persec

        # Rest of land runoff calculation

        append(developed_values, cultivated_values) -> crop_dev_vals

        crop_reclass_table %>%
          filter(!CDL_ID %in% crop_dev_vals) %>%
          .[["CDL_ID"]] -> other_vals

        sup(get_runoff_values(cropcover_agg,
                                           runoff_agg,
                                           other_vals,
                                           polygon_area,
                                           land_table)) -> other_runoff_m3persec

        # Total runoff calculation
        append(crop_dev_vals, other_vals) -> all_vals

        # Calculate fractions
        dev_runoff_m3persec + cultivated_runoff_m3persec + other_runoff_m3persec -> expec_total_runoff

        }else{
          expec_total_runoff <- 0
          dev_runoff_m3persec <- 0
          cultivated_runoff_m3persec <- 0
        }

        #--------------------------------------------------------
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
          dplyr::right_join(irrigation_km2, by = "GCAM_commodity") %>%
          mutate(consumption_BCM = Irr_BCM_KM2 * irr) %>%
          .[["consumption_BCM"]] %>% sum(na.rm = T) ->
          total_irr_bcm

        #--------------------------------------------------------
        # TELECONNECTION - Count # of economic sectors within watershed
        # get raster values of land use types within the watershed.
        get_zonal_data(economic_USA, watersheds_select) %>%
          # merge ids with id table to attach class names
          get_nlud_names() -> nlud_table


        # TODO:  Potentally add in interbasin transfer from commit 1989d3877304370d034091e7e0a6598c947c6c4a
          n_transfers_into <- NA_integer_
          n_transfers_out <- NA_integer_
          n_transfers_within <- NA_integer_

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

        runoff_totals %>%
          .[["flow_BCM"]] -> runoff_total_ts

        sup(yield(Q = runoff_total_ts,
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
        wasteflow_points[watersheds_select, ] %>% .@data %>% as_tibble() ->
          wwtp_points

        ## manual catch for Norfolk | VA (WWTP falls outside of watershed 4227),
        ## and Aurora | CO (WWTP falls outside of 2349)
        ## and Colorado Springs | CO (WWTP just below dam crest for watershed 3689)

        if(watershed %in% c(4227, 2349)) wwtp_points <- tibble()
        if(watershed == 3689){
          wwtp_points <-
            wwtp_points %>% filter(lat > 39)
        }


        wwtp_points[["flow_cumecs"]] %>% sum() -> wastewater_outflow_m3persec
        nrow(wwtp_points %>% unique()) -> n_treatment_plants

        #---------------------------------------------------------
        # count the facilities present in the watershed
        epa_facilities[watersheds_select, ] %>% .@data -> facilities_present
        n_ag_crops <- facilities_present %>% filter(sector == "ag_crops") %>% nrow()
        n_ag_livestock <- facilities_present %>% filter(sector == "ag_livestock") %>% nrow()
        n_construction_and_manufacturing <- facilities_present %>% filter(sector == "construction_and_manufacturing") %>% nrow()
        n_mining <- facilities_present %>% filter(sector == "mining") %>% nrow()
        n_other <- facilities_present %>% filter(sector == "other") %>% nrow()
        n_transport_and_utilities <- facilities_present %>% filter(sector == "transport_and_utilities") %>% nrow()
        n_oil_and_gas_extraction <- facilities_present %>% filter(sector == "oil_and_gas_extraction") %>% nrow()
        n_unspecified <- facilities_present %>% filter(sector == "unspecified") %>% nrow()
        n_facilities_total <- facilities_present %>% nrow()

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
              "outlow treatment plants",       n_treatment_plants,
              "ag_crops_facilities",           n_ag_crops,
              "ag_livestock_facilities",       n_ag_livestock,
              "construction_and_manufacturing",n_construction_and_manufacturing,
              "mining",                        n_mining,
              "oil_and_gas",                   n_oil_and_gas_extraction,
              "facilities_total",              n_facilities_total
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
              "average historical runoff",       "m3/sec",   historical_runoff_mean_m3sec,
              "average historical flow",         "m3/sec",   historical_flow_mean_m3sec,
              "driest month average runoff",     "m3/sec",   historical_runoff_dry_m3sec,
              "driest month average flow",       "m3/sec",   historical_flow_dry_m3sec,
              "wastewater discharge",            "m3/sec",   wastewater_outflow_m3persec
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
