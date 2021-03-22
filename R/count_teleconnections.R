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
#' @import rgdal
#' @export
count_watershed_teleconnections <- function(data_dir,
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
                                              population = "land/pden2010_block/pden2010_60m.tif",
                                              runoff = "water/UWSCatch/USA_Mean_Runoff.tif",
                                              nhd_flow = "water/UWSCatCH/UWSCatCH_Intake_Flows.shp",
                                              contributions = "water/UWSCatCH/Watershed_Contributions.csv"
                                            )){

  suppressWarnings(count_watershed_data(data_dir = data_dir,
                       cities = cities,
                       run_all = run_all,
                       file_paths = file_paths)) ->
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
    watersheds


  # read power plant data
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

  # Get source contributions
  get_source_contribution(data_dir,file_paths) -> contributions

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

        # get facility counts per unit area
        map(city_watershed_data, function(x){
          x$counts %>% filter(item %in% c("ag_crops_facilities",
                                          "ag_livestock_facilities",
                                          "construction_and_manufacturing",
                                          "mining",
                                          "oil_and_gas",
                                          "facilities_total"))
        }) %>% bind_rows() %>% group_by(item) %>% summarise(n = sum(count)) ->
          total_facility_counts
        filter(total_facility_counts, item == "ag_crops_facilities") %>% .[["n"]] -> n_fac_agcrop
        filter(total_facility_counts, item == "ag_livestock_facilities") %>% .[["n"]] -> n_fac_aglivestock
        filter(total_facility_counts, item == "construction_and_manufacturing") %>% .[["n"]] -> n_fac_cnsmnf
        filter(total_facility_counts, item == "mining") %>% .[["n"]] -> n_fac_mining
        filter(total_facility_counts, item == "oil_and_gas") %>% .[["n"]] -> n_fac_oilgas
        filter(total_facility_counts, item == "facilities_total") %>% .[["n"]] -> n_fac_total

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

        # watershed_population
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "watershed population") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> watershed_pop

        # population water consumption
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "population water consumption") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> pop_cons_ltr_day

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

        # get contribution of different source watersheds to city supply
        contributions %>%
          filter(city_state == city) %>%
          mutate(DVSN_ID = as.character(DVSN_ID)) -> source_contributions

        # runoff contaminant concentrations
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "total runoff from watershed") %>%
            .[["value"]] -> total_runoff
          x$metrics %>% filter(metric == "development runoff") %>%
            .[["value"]] -> dev_runoff
          x$metrics %>% filter(metric == "cultivated runoff") %>%
            .[["value"]] -> cropland_runoff

          tibble(total_runoff, dev_runoff, cropland_runoff) %>%
            mutate(nonpristine_conc = (dev_runoff + cropland_runoff) / total_runoff,
                   ag_conc = cropland_runoff / total_runoff,
                   dev_conc = dev_runoff / total_runoff) %>%
            select(nonpristine_conc, ag_conc, dev_conc)

        }) -> runoff_contaminant_concentrations

        names(runoff_contaminant_concentrations) %>%
          map(function(ws){
            runoff_contaminant_concentrations[[ws]] %>%
              mutate(DVSN_ID = ws)
          }) %>% bind_rows() %>%
          right_join(source_contributions, by = "DVSN_ID") ->
          runoff_contaminant_and_contributions

        get_average_runoff_conc <- function(metric){
          runoff_contaminant_and_contributions %>%
            select(x = one_of(metric), contribution_to_supply) %>%
            tidyr::replace_na(list(x = 0)) %>%
            mutate(conc = x * contribution_to_supply) %>% .[["conc"]] %>% sum()
        }


        get_average_runoff_conc_exgw <- function(metric){
          runoff_contaminant_and_contributions %>%
            select(x = one_of(metric), contribution_to_supply) %>%
            filter(!is.na(x)) %>%
            mutate(conc = x * contribution_to_supply) -> concs
          sum(concs[["conc"]]) / sum(concs[["contribution_to_supply"]])
        }

        get_average_runoff_conc_exgw_unweighted <- function(metric){
          runoff_contaminant_and_contributions %>%
            select(x = one_of(metric)) %>%
            filter(!is.na(x)) %>% .[["x"]] %>% mean()
        }

        get_max_runoff_conc <- function(metric){
          runoff_contaminant_and_contributions %>%
            select(x = one_of(metric)) %>%
            .[["x"]] %>% max(na.rm = T)
        }

        ag_runoff_max <- get_max_runoff_conc("ag_conc")
        ag_runoff_av_exgw <- get_average_runoff_conc_exgw("ag_conc")
        ag_runoff_av <- get_average_runoff_conc("ag_conc")

        dev_runoff_max <- get_max_runoff_conc("dev_conc")
        dev_runoff_av_exgw <- get_average_runoff_conc_exgw("dev_conc")
        dev_runoff_av <- get_average_runoff_conc("dev_conc")

        np_runoff_max <- get_max_runoff_conc("nonpristine_conc")
        np_runoff_av_exgw <- get_average_runoff_conc_exgw("nonpristine_conc")
        np_runoff_av_exgw_unweighted <- get_average_runoff_conc_exgw_unweighted("nonpristine_conc")
        np_runoff_av <- get_average_runoff_conc("nonpristine_conc")

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


        source_contributions %>%
          filter(type == "surface water") %>%
          .[["contribution_to_supply"]] %>% sum() * 100 -> surface_contribution_pct


        # runoff, flow, and wastewater plant discharge
        map(city_watershed_data, function(x){
          x$metrics %>%
            filter(metric %in% c("average historical runoff",
                                 "average historical flow",
                                 "driest month average runoff",
                                 "driest month average flow")) %>%
            select(-unit)
        }) -> runoff_and_flows

        # wastewater plant discharge
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "wastewater discharge") %>%
            .[["value"]]
        }) %>% as_tibble() %>% gather(DVSN_ID, wastewater_discharge_m3sec) ->
          wastewater_discharge_m3sec

        names(runoff_and_flows) %>%
          map_dfr(function(DVSN){
            runoff_and_flows[[DVSN]] %>%
              mutate(DVSN_ID = !! DVSN) %>%
              left_join(wastewater_discharge_m3sec, by = "DVSN_ID")
          }) %>%
          mutate(conc_pct = 100 * (wastewater_discharge_m3sec / (value + wastewater_discharge_m3sec))) %>%
          select(DVSN_ID, metric, conc_pct) %>%
          split(.$metric) %>%
          map(function(x) right_join(x, source_contributions, by = "DVSN_ID")) ->
          flow_and_runoff_metrics_and_contributions

        # mean concentration of surface water supply
        flow_and_runoff_metrics_and_contributions %>%
          map(function(x){
            x %>% filter(type == "surface water") %>%
              summarise(conc = sum(conc_pct * contribution_to_supply) / sum(contribution_to_supply)) %>%
              .[["conc"]]
          }) -> surface_water_concentrations

        # unweighted mean concentration of surface water supply
        flow_and_runoff_metrics_and_contributions %>%
          map(function(x){
            x %>% filter(type == "surface water") %>%
              .[["conc_pct"]] %>% mean()
          }) -> surface_water_concentrations_unweighted

        # mean concentration of all water supply
        flow_and_runoff_metrics_and_contributions %>%
          map(function(x){
            x %>%
              mutate(conc_pct = if_else(is.na(conc_pct) & type != "surface water", 0, conc_pct)) %>%
              summarise(conc = sum(conc_pct * contribution_to_supply) / sum(contribution_to_supply)) %>%
              .[["conc"]]
          }) -> all_water_concentrations

      # concentration from worst-case surface water supply
        flow_and_runoff_metrics_and_contributions %>%
          map(function(x){
            x %>% filter(type == "surface water") %>%
              summarise(conc = max(conc_pct)) %>% .[["conc"]]
          }) -> surface_water_concentration_worst

        # importance of worst-case surface water supply DVSN
        flow_and_runoff_metrics_and_contributions$`average historical flow` %>%
          filter(type == "surface water") %>%
          filter(conc_pct == max(conc_pct)) %>%
          .[["contribution_to_supply"]] %>% dplyr::first() * 100 ->
          importance_of_worst_watershed_pct

        # treatment plants
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "outflow treatment plants") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> n_treatment_plants

        # population-based discharge
        pop_cons_ltr_day * ltrday_to_m3sec -> pop_cons_m3sec

        # nhd flow
        map(city_watershed_data, function(x){
          x$metrics %>% filter(metric == "nhd flow") %>%
            .[["value"]]
        }) %>% unlist() %>% sum() -> nhd_average_flow


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
                 n_climate_zones,
                 n_transfers_in,
                 n_transfers_out,
                 n_transfers_within,
                 n_hydro_plants,
                 n_thermal_plants,
                 n_fac_agcrop,
                 n_fac_aglivestock,
                 n_fac_cnsmnf,
                 n_fac_mining,
                 n_fac_oilgas,
                 n_fac_total,
                 hydro_gen_MWh,
                 thermal_gen_MWh,
                 thermal_cons_BCM,
                 thermal_with_BCM,
                 n_utilities,
                 n_ba,
                 n_crop_classes,
                 cropland_fraction,
                 #cropland_runoff_percent,
                 developed_fraction,
                 #developed_runoff_percent,
                 ag_runoff_max,
                 ag_runoff_av_exgw,
                 ag_runoff_av,
                 dev_runoff_max,
                 dev_runoff_av_exgw,
                 dev_runoff_av,
                 np_runoff_max,
                 np_runoff_av_exgw,
                 np_runoff_av_exgw_unweighted,
                 np_runoff_av,
                 n_economic_sectors,
                 max_withdr_dist_km,
                 avg_withdr_dis_km,
                 n_treatment_plants,
                 watershed_pop,
                 pop_cons_m3sec,
                 av_fl_sur_conc_pct = surface_water_concentrations[["average historical flow"]],
                 av_fl_sur_conc_pct_unweighted = surface_water_concentrations_unweighted[["average historical flow"]],
                 #dr_fl_sur_conc_pct = surface_water_concentrations[["driest month average flow"]],
                 av_ro_sur_conc_pct = surface_water_concentrations[["average historical runoff"]],
                 #dr_ro_sur_conc_pct = surface_water_concentrations[["driest month average runoff"]],
                 av_fl_all_conc_pct = all_water_concentrations[["average historical flow"]],
                 #dr_fl_all_conc_pct = all_water_concentrations[["driest month average flow"]],
                 av_ro_all_conc_pct = all_water_concentrations[["average historical runoff"]],
                 #dr_ro_all_conc_pct = all_water_concentrations[["driest month average runoff"]],
                 av_fl_max_conc_pct = surface_water_concentration_worst[["average historical flow"]],
                 #dr_fl_max_conc_pct = surface_water_concentration_worst[["driest month average flow"]],
                 av_ro_max_conc_pct = surface_water_concentration_worst[["average historical runoff"]],
                 #dr_ro_max_conc_pct = surface_water_concentration_worst[["driest month average runoff"]],
                 surface_contribution_pct,
                 importance_of_worst_watershed_pct)
        )
      }

    })
}
