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
                                            cities = NULL,
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
                                              runoff = "water/Historical_Mean_Runoff/Historical_Mean_Runoff.tif"
                                            )){

  count_watershed_data(data_dir = data_dir,
                       cities = cities,
                       file_paths = file_paths) ->
    watershed_data

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
  import_shapefile(paste0(data_dir, file_paths["watersheds"]),
                   method = "rgdal") %>%
    # subset for desired watersheds
    subset(DVSN_ID %in% city_watershed_mapping$DVSN_ID) ->
    # remove City_Name; it's inconsistent with other shapefiles
    #dplyr::select(-City_Name) ->
    watersheds

  # # read NID point file and select only Flood Control Dams (C = Flood Control)
  # import_shapefile(paste0(data_dir, dams_file_path)) %>%
  #   subset(grepl("C", Purposes)) %>%
  #   as_Spatial() -> flood_control_dams

  get_ucs_power_plants(paste0(data_dir, file_paths["powerplants"])) ->
    power_plants_usa

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, file_paths["citypoint"]),
                   method = "rgdal") %>% st_as_sf() %>%
    rename("city_uid" = "City_ID") %>%
    subset(city_uid %in% city_watershed_mapping$city_uid) -> city_points

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, file_paths["withdrawal"]),
                   method = "rgdal") %>% st_as_sf() -> withdrawal_points

  # Get watershed usage table
  get_watershed_usage(cities) -> watershed_usage_table

  # Get Population
  get_population() -> cities_population

  # map through all cities, computing teleconnections
  cities %>%
    map_dfr(function(city){
      filter(city_watershed_mapping, city_state == !!city) -> city_select

      city_select %>% .[["DVSN_ID"]] -> city_intake_ids

      watersheds %>%
        subset(DVSN_ID %in% city_intake_ids) -> watersheds_city

      # catch cases with only groundwater points (i.e., no watershed polygons)
      if(nrow(watersheds_city) == 0){

        # city population
        cities_population %>%
          filter(city_state == !! city) %>%
          .[["Population"]] -> city_population

        done(city)
        return(
          tibble(city,
                 city_population)
        )
      }else{

        watershed_data[which(names(watershed_data) %in% city_intake_ids)] ->
          city_watershed_data

        power_plants_city <- power_plants_usa[watersheds_city,]

        # city population
        cities_population %>%
          filter(city_state == !! city) %>%
          .[["Population"]] -> city_population

        # number of watersheds
        length(city_watershed_data) -> n_watersheds

        # maximum distance between city's watershed(s)
        city_points %>%
          subset(city_uid %in% city_select$city_uid) %>%
          as_Spatial() -> city_centroid

        withdrawal_points %>%
          subset(DVSN_ID %in% city_intake_ids) %>%
          as_Spatial() -> withdrawal_city

        distGeo(withdrawal_city, city_centroid) %>%
          as_tibble() %>%
          .[["value"]] -> distance_values

        distance_values %>% max()/m_to_km -> max_withdr_dist_km

        distance_values %>% mean()/m_to_km -> avg_withdr_dis_km

        # Get number of cities
        watershed_usage_table %>%
          subset(Watershed_DVSN_ID %in% city_intake_ids) -> city_usage
        if(nrow(city_usage) == 0){
          dependent_city_pop <- 0
          n_other_cities <- 0
        }else{
          city_usage %>% select(city_state) %>% unique() -> select_cities

          nrow(select_cities) -> n_other_cities

          left_join(select_cities, cities_population, by = "city_state") %>%
            .[["Population"]] %>%
            sum() -> dependent_city_pop
        }
        # number of climate zones
        map(city_watershed_data, function(x){
          x$climate_zones
        }) %>% unlist() %>% unique() %>% length() -> n_climate_zones

        # number of transfers in
        map(city_watershed_data, function(x){
          x$counts %>% filter(item == "transfers in") %>%
            .[["count"]]
        }) %>% unlist() %>% sum() -> n_transfers_in

        # number of transfers within
        map(city_watershed_data, function(x){
          x$counts %>% filter(item == "transfers within") %>%
            .[["count"]]
        }) %>% unlist() %>% sum() -> n_transfers_within

        # number of transfers out
        map(city_watershed_data, function(x){
          x$counts %>% filter(item == "transfers out") %>%
            .[["count"]]
        }) %>% unlist() %>% sum() -> n_transfers_out

        # total_watershed_area
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "watershed area") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> watershed_area_sqkm

        # total storage
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "total reservoir storage") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> storage_BCM

        # total yield
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "total watershed yield") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> yield_BCM

        # total irrigation consumption
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "total irrigation consumption") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> irr_cons_BCM

        # population water consumption
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "population") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> pop

        # population water consumption
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "population water consumption") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> pop_cons_ltr_sqkm

        # hydro plants
        map(city_watershed_data, function(x){
          x$counts %>% filter(item == "hydro plants") %>%
            .[["count"]]
        }) %>% unlist() %>% sum() -> n_hydro_plants

        # thermal plants
        map(city_watershed_data, function(x){
          x$counts %>% filter(item == "thermal plants") %>%
            .[["count"]]
        }) %>% unlist() %>% sum() -> n_thermal_plants

        # hydro gen
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "total hydro generation") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> hydro_gen_MWh

        # thermal gen
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "total thermal generation") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> thermal_gen_MWh

        # thermal consumption
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "total thermal water consumption") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> thermal_cons_BCM

        # thermal withdrawal
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "total thermal water withdrawals") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> thermal_with_BCM

        # cultivated runoff percent
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "runoff % from cropland") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> cropland_runoff_fraction

        # developed runoff percent
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "runoff % from developed land") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> developed_runoff_fraction

        # number of utilities
        power_plants_city %>%
          subset(cooling == "Yes" | `Power Plant Type` == "Hydro") %>%
          .[["UTILITY_ID"]] %>% unique() %>%
          length() -> n_utilities

        # number of balancing authorities
        power_plants_city %>%
          subset(cooling == "Yes" | `Power Plant Type` == "Hydro") %>%
          .@data %>% dplyr::select(PLANT_CODE, CNTRL_AREA) %>% unique() %>%
          .[["CNTRL_AREA"]] -> tc_ba_na
        sum(is.na(tc_ba_na)) -> n_missing_ba
        if(n_missing_ba > 0) message(paste0("For ", city,", ", n_missing_ba,
                                            " plant(s) not assigned a Balancing Authority"))
        tc_ba_na %>% .[!is.na(.)] %>% unique() %>% length() -> n_ba


        # TELECONNECTION - NUMBER OF FLOOD CONTROL DAMS WITHIN WATERSHED.
        # flood_control_dams[watersheds_city, ] %>%
        #   length() -> tc_fcdam

        # number of GCAM crop classes present
        map(city_watershed_data, function(x){
          x$land %>% filter(is_crop == TRUE, cell_freq > 0) %>%
            .[["GCAM_ID"]]
        }) %>% unlist() %>% unique() %>% length() -> n_crop_classes

        # area of crop land as fraction of total land area
        map(city_watershed_data, function(x){
          x$land %>% .[["cell_freq"]] %>% sum(na.rm = T)
        }) %>% unlist() %>% sum() -> total_cell_area

        map(city_watershed_data, function(x){
          x$land %>% filter(is_crop == TRUE) %>%
            .[["cell_freq"]] %>% sum(na.rm = T)
        }) %>% unlist() %>% sum() -> crop_cell_area

        # developed / area
        map(city_watershed_data, function(x){
          x$land %>% filter(grepl("Developed", CDL_Class)) %>%
            .[["cell_freq"]] %>% sum(na.rm = T)
        }) %>% unlist() %>% sum() -> dev_cell_area

        crop_cell_area / total_cell_area -> cropland_fraction
        dev_cell_area / total_cell_area -> developed_fraction

        # economic sectors
        map(city_watershed_data, function(x){
          x$economic_sectors %>%
            filter(Reclass != "Water") %>%
            .[["Reclass"]] %>% unique()
        }) %>% unlist() %>% unique() %>% length() -> n_economic_sectors


        done(city)

        # output
        return(
          tibble(city,
                 city_population,
                 n_watersheds,
                 n_other_cities,
                 dependent_city_pop,
                 watershed_area_sqkm,
                 storage_BCM,
                 yield_BCM,
                 irr_cons_BCM,
                 pop,
                 n_climate_zones,
                 n_transfers_in,
                 n_transfers_out,
                 n_transfers_within,
                 n_hydro_plants,
                 n_thermal_plants,
                 hydro_gen_MWh,
                 thermal_gen_MWh,
                 thermal_cons_BCM,
                 thermal_with_BCM,
                 n_utilities,
                 n_ba,
                 n_crop_classes,
                 cropland_fraction,
                 cropland_runoff_fraction,
                 developed_fraction,
                 developed_runoff_fraction,
                 n_economic_sectors,
                 max_withdr_dist_km,
                 avg_withdr_dis_km)
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
