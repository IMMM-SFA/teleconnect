#' count_watershed_teleconnections
#'
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param watersheds_file_path path of watersheds shapefile within data_dir
#' @param landcover_file_path path of NLCD landcover data file
#' @param powerplants_file_path path of power plants data file
#' @param cities a vector of cities to be included in the count. If omitted, all cities will be included.
#' @details counts teleconnections assoicated with water supply catchments associated with each city
#' @importFrom purrr map_dfr
#' @importFrom dplyr filter
#' @importFrom tibble tibble
#' @export
#'
count_watershed_teleconnections <- function(data_dir,
                                            watersheds_file_path = "water/CWM_v2_2/World_Watershed8.shp",
                                            landcover_file_path = "land/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img",
                                            powerplants_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
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

  # read landcover raster for US
  import_raster(paste0(data_dir, landcover_file_path)) ->
    landcover_USA

  # map through all cities, computing teleconnections
  cities %>%
    map_dfr(function(city){
      filter(city_watershed_mapping, city_state == !!city) %>%
        .[["DVSN_ID"]] -> city_intake_ids

      # subset watersheds shapefile for target city
      watersheds %>%
        subset(DVSN_ID %in% city_intake_ids) ->
        watersheds_city

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

      # TELECONNECTION - NUMBER OF LAND USE CLASSES
      get_raster_val_classes(landcover_USA, watersheds_city) %>%
        length() ->
        tc_n_landclasses

      done(city)

      # output
      return(
        tibble(city = !! city,
               n_watersheds = tc_n_watersheds,
               n_hydro = tc_n_hydroplants,
               n_thermal = tc_n_thermalplants,
               n_landclasses = tc_n_landclasses)
      )
    })
}
