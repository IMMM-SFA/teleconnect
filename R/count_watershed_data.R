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
#' @details counts teleconnections assoicated with water supply catchments associated with each city
#' @importFrom purrr map_dfr map
#' @importFrom dplyr filter group_indices left_join
#' @importFrom tibble tibble
#' @importFrom sf as_Spatial st_as_sf st_cast
#' @importFrom foreign read.dbf
#' @importFrom exactextractr exact_extract
#' @importFrom geosphere areaPolygon
#' @importFrom tmaptools set_projection
#' @importFrom lwgeom st_startpoint st_endpoint
#' @import dams
#' @export
count_watershed_data <- function(data_dir,
                                      watersheds_file_path = "water/CWM_v2_2/World_Watershed8.shp",
                                      powerplants_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
                                      crop_file_path = "land/2016_30m_cdls/cdl_lowres_usa.img",
                                      crop_attribute_path = "land/2016_30m_cdls/cdl_lowres_usa.img.vat.dbf",
                                      irrigation_file_path = "land/usa_demeter.csv",
                                      nlud_file_path = "land/usa_nlud_LR.tif",
                                      hydro_file_path = "energy/EHA_Public_PlantFY2019_GIS_6/ORNL_EHAHydroPlant_PublicFY2019final.xlsx",
                                      transfers_file_path = "water/USIBTsHUC6_Dickson.shp",
                                      watersheds = NULL){
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

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, watersheds_file_path),
                   method = "rgdal") %>%
    # subset for desired watersheds
    subset(DVSN_ID %in% watershed_mapping$DVSN_ID) ->
    # remove City_Name; it's inconsistent with other shapefiles
    #dplyr::select(-City_Name) ->
    catchment_shapes

  # read ucs plant data
  get_ucs_power_plants(paste0(data_dir, powerplants_file_path)) -> power_plants_USA

  # read croptype raster for US
  import_raster(paste0(data_dir, crop_file_path)) -> cropcover_USA

  read.dbf(paste0(data_dir,crop_attribute_path)) -> crop_cover_levels

  # read reclassified crop table
  reclassify_raster(crop_cover_levels) -> crop_reclass_table

  # read NLUD economic raster
  import_raster( paste0(data_dir, nlud_file_path)) -> economic_USA

  # read NID point file
  dams::nid_cleaned -> nid_dataset

  # read hydrosource hydropower dataset
  get_hydro_dataset(data_dir, hydro_file_path) %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = "+proj=longlat +datum=WGS84 +no_defs") %>% as_Spatial() -> hydro_points

  #read demeter irrigation file
  get_demeter_file(paste0(data_dir, irrigation_file_path)) -> usa_irrigation

  import_shapefile(paste0(data_dir, transfers_file_path), method = "rgdal",) %>%
    st_as_sf() %>% st_transform("+proj=longlat +datum=WGS84 +no_defs") -> interbasin_transfers

  # map through all cities, computing teleconnections
  watersheds %>%
    map_dfr(function(watershed){
      # subset watersheds shapefile for target city
      catchment_shapes %>%
        subset(DVSN_ID %in% watershed) -> watersheds_select
        #----------------------------------------------------
        # TELECONNECTION - NUMBER OF CITIES USING WATERSHED
        watershed_mapping %>%
          filter(DVSN_ID == watershed) %>%
          .[["DVTR_ID"]] %>% length() -> tc_n_cities

        #----------------------------------------------------
        # TELECONNECTION - HYDRO ENERGY
        # subset power plants for target watersheds
        power_plants_USA[watersheds_select, ] %>% as.data.frame() -> watershed_power_plants
        subset(watershed_power_plants, Power.Plant.Type == "Hydropower") -> hydro_plants
        if(nrow(hydro_plants) == 0){
          hydro_generation= NA_character_
          hydro_storage = NA_character_
        }else{
        # get area of watershed
        areaPolygon(watersheds_select) -> polygon_area

        # Generation per unit watershed area MWh/km2
        sum(hydro_plants$MWh) -> MWh
        MWh / polygon_area -> hydro_generation

        # subset hydro dams for target watershed
        hydro_points[watersheds_select, ] %>% as.data.frame() -> watershed_dams
        left_join(watershed_dams,nid_dataset, by = "NID_ID") %>%
          select(hydro_variable) %>% na.omit() -> hydro_dams
        # convert acre-feet to km3 and add up all storage
        sum(hydro_dams[1])*acrefeet_conv -> dam_storage_km3
        # Storage per unit watershed area km3/km2
        dam_storage_km3 / polygon_area -> hydro_storage
        }
        #---------------------------------------------------
        # TELECONNECTION - THERMAL ENERGY
        # subset power plants for target watersheds
        subset(watershed_power_plants, cooling == "Yes") %>%
          mutate(consumption = as.integer(consumption)) %>%
          mutate(withdrawal = as.integer(withdrawal))-> thermal_plants
        if(nrow(thermal_plants) == 0){
          thermal_generation = NA_character_
          thermal_consumption = NA_character_
          thermal_withdrawal = NA_character_
        }else{
        # get area of watershed
        areaPolygon(watersheds_select) -> polygon_area

        # Generation per unit watershed area MWh/km2
        sum(thermal_plants$MWh) -> MWh
        MWh / polygon_area -> thermal_generation
        # Consumption per unit area Mg/yr/km2
        sum(thermal_plants$consumption) -> total_consum
        total_consum / polygon_area -> thermal_consumption
        # Withdrawal per unit area Mg/yr/km2
        sum(thermal_plants$withdrawal) -> total_withdrawal
        total_withdrawal / polygon_area -> thermal_withdrawal
        }
        #---------------------------------------------------

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
        usa_irrigation[watersheds_select, ] %>%
          as_tibble() -> irrigation_city
        # count the number of crop types that are irrigated
        get_irrigation_count(irrigation_city) %>%
          filter(GCAM_Class %in% crop_and_landcover_types$GCAM_Class) -> irr_crops
        length(irr_crops$irr_count) -> tc_n_irrigated_crops

        # TELECONNECTION - Count # of economic sectors within watershed
        # get raster values of land use types within the watershed.
        get_zonal_data(economic_USA, watersheds_select) -> economic_ids

        # merge ids with id table to attach class names
        get_nlud_names(economic_ids) -> nlud_table
        # Filter out water and count the unique class variables within the watershed
        nlud_table %>%
          filter(Reclass != "Water") %>%
          .[["Reclass"]] %>% unique() %>%
          length() -> tc_n_economicsectors

        # -------------------------------------------------------
        # TELECONNECTION - Interbasin Transfers in Watershed
        st_as_sf(watersheds_select) -> sf_watershed
        suppressWarnings(suppressMessages(sf::st_intersection(interbasin_transfers, sf_watershed))) %>%
         as.data.frame %>% .[c("OBJECTID")] -> watershed_transfers

        subset(interbasin_transfers, OBJECTID %in% watershed_transfers$OBJECTID) -> transfers

        st_startpoint(transfers) %>%
          as.data.frame() %>% st_as_sf() -> transfer_start
        transfer_start$observation <- paste0(1:nrow(transfer_start),"start")

        st_endpoint(transfers) %>%
          as.data.frame() %>% st_as_sf() -> transfer_end
        transfer_end$observation <- paste0(1:nrow(transfer_end),"end")

        transfer_points <- rbind(transfer_start,transfer_end) %>% as.data.frame()

        suppressWarnings(suppressMessages(st_intersection(sf_watershed, transfer_start))) -> transfer_out
        suppressWarnings(suppressMessages(st_intersection(sf_watershed, transfer_end))) -> transfer_in
        suppressWarnings(suppressMessages(st_intersection(sf_watershed, transfer_points))) -> transfer_both

        library(ggplot2)
        ggplot() +
          geom_sf(data = sf_watershed) +
          geom_sf(data = transfers, col = "purple") +
          geom_sf(data = transfer_start, col = "green") +
          geom_sf(data = transfer_end, col = "red")
        #---------------------------------------------------------

        done(watershed)

                tibble(DVSN_ID = !! watershed,
                       hydro_generation,
                       hydro_storage,
                       thermal_consumption,
                       thermal_generation,
                       thermal_withdrawal,
                       wtrshd_impact = watershed_condition,
                       n_cropcover = tc_n_cropcover,
                       n_irrigatedcrops = tc_n_irrigated_crops,
                       n_economicsectors = tc_n_economicsectors) -> watershed_table
                watershed_mapping[c(7,6,4,8,5)] %>%
                  mutate(., DVSN_ID = as.character(DVSN_ID))-> water_mapping_select
                watershed_table %>% mutate(., DVSN_ID = as.character(DVSN_ID)) -> watershed_table_df
                left_join(water_mapping_select, watershed_table_df, by = "DVSN_ID") -> city_watershed_data
        return(city_watershed_data)
      }
    )

  }

