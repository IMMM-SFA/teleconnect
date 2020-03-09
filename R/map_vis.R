#' plot_watershed
#'
#' @details plots the watersheds and relevant features for a chosen city
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param watersheds_file_path path of watersheds shapefile within data_dir
#' @param landcover_file_path path of NLCD landcover data file
#' @param powerplants_file_path path of power plants data file
#' @param crop_file_path path of CDL img file
#' @param conus_file_path path of CONUS state boundaries shape file
#' @param city character. City to be analysed in the format "City | STATE" (e.g., "Phoenix | AZ")
#' @import tmap
#' @importFrom tmaptools aggregate_map
#' @importFrom dplyr filter mutate if_else right_join
#' @importFrom tibble rowid_to_column
#' @importFrom raster levels extent
#' @export
#'
plot_watershed <- function(data_dir,
                           watersheds_file_path = "water/CWM_v2_2/World_Watershed8.shp",
                           crop_file_path = "land/2016_30m_cdls/2016_30m_cdls.img",
                           powerplants_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
                           conus_file_path = "misc/states_21basic/states.shp",
                           dams_file_path = "water/nabd_fish_barriers_2012/nabd_fish_barriers_2012.shp",
                           nlud_file_path = "land/usa_nlud_raster.tif",
                           citypoint_file_path = "water/CWM_v2_2/City_Centroid.shp",
                           city){

  all_cities <- get_cities()[["city_state"]]

  # throw error if any chosen city lies outside available cities
  if(!city %in% all_cities) {
    stop(paste0(paste(city), " is not part of 'teleconnect'!"))
  }

  # get city watershed_ids
  get_cities() %>%
    filter(city_state == !!city) %>%
    subset(key_watershed == TRUE) ->
    city_watershed_mapping

  # Load city centroid file and merge with city mapping file
  import_shapefile(paste0(data_dir, citypoint_file_path)) %>% st_as_sf() %>%
    rename("city_uid" = "City_ID") %>%
    left_join(city_watershed_mapping, by = "city_uid") %>%
    select(city_state, geometry) %>% unique() %>%
    subset(city_state %in% city) -> city_point

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, watersheds_file_path),
                   method = "rgdal") %>%
    # subset for desired watersheds
    subset(DVSN_ID %in% city_watershed_mapping$DVSN_ID) ->
    watersheds_city

  # get state shape
  import_shapefile(paste0(data_dir, conus_file_path),
                   method = "rgdal") %>%
    subset(STATE_ABBR == substr(city, nchar(city) - 1, nchar(city))) ->
    state_shape

  # read ucs plant data
  get_ucs_power_plants(paste0(data_dir, powerplants_file_path)) ->
    power_plants_USA

  power_plants_USA@data <- power_plants_USA@data %>%
    mutate(cooling_tech = if_else(cooling_tech == "#not designated",
                                  "None", cooling_tech)) %>%
    mutate(`Power Plant Type` = if_else(cooling == "Yes", "Thermal",
                                        `Power Plant Type`)) %>%
    mutate(`Power Plant Type` = as.factor(`Power Plant Type`),
           `Cooling Technology` = as.factor(cooling_tech))

  power_plants_USA[watersheds_city, ] -> power_plants_city

  # read NID point file and select only Flood Control Dams (C = Flood Control)
  import_shapefile(paste0(data_dir, dams_file_path)) %>%
    #subset(grepl("C", Purposes)) %>%
    as_Spatial() -> flood_control_dams

  # read crop raster for US
  import_raster(paste0(data_dir, crop_file_path)) ->
    cropcover_USA

  # get crop classes and color scheme
  levels(cropcover_USA)[[1]] %>%
    filter(Class_Names != "") %>%
    mutate(color = rgb(Red, Green, Blue, maxColorValue = 255)) %>%
    dplyr::select(ID, cropclass = Class_Names, color) %>%
    mutate(color = if_else(ID == 0, "#ffffff", color)) %>%
    mutate(cropclass = as.character(cropclass)) ->
    cropclasses

  # make color palette
  setNames(cropclasses$color, cropclasses$cropclass) -> crop_palette

  # mask to required city
  cropcover_USA %>%
    mask_raster_to_polygon(watersheds_city) ->
    cropcover_city_watersheds
  rm(cropcover_USA)

  rat_cc <- levels(cropcover_city_watersheds)[[1]]
  rat_cc <- right_join(rat_cc, cropclasses, by = "ID") %>%
    dplyr::select(ID, cropclass) %>%
    mutate(cropclass = if_else(ID == 131, "Barren_", cropclass),
           cropclass = if_else(ID == 152, "Shrubland_", cropclass)) %>%
    mutate(cropclass = as.factor(cropclass))
  levels(cropcover_city_watersheds) <- rat_cc

  # aggregate raster to increase plot speed (and also convert to format usable with tmap)
  cropcover_city_watersheds %>%
    aggregate_map(fact = 4, agg.fun = "modal") ->
    cropcover_city_watersheds_agg

  # correct raster attribute table in aggregated raster
  levels(cropcover_city_watersheds_agg) %>% .[[1]] %>%
    mutate(ID = rat_cc$ID) -> rat_cc_agg
  levels(cropcover_city_watersheds_agg) <- rat_cc_agg

  # get the coordinate system of the crop raster and transform watershed
  r_crs <- st_crs(projection(cropcover_city_watersheds_agg))
  watersheds_city_trans <- st_as_sf(watersheds_city) %>% st_transform(crs = r_crs)
  state_shape_trans <- st_as_sf(state_shape) %>% st_transform(crs = r_crs)

  # create extents to define bounding boxes of maps
  watersheds_city_trans %>% extent() -> ws_extent
  state_shape_trans %>% extent() -> state_extent
  extent_y_diff <- (ws_extent@ymax - ws_extent@ymin) * 0.05
  extent_x_diff <- (ws_extent@xmax - ws_extent@xmin) * 0.05

  extent(
    ws_extent@xmin - extent_x_diff,
    ws_extent@xmax + extent_x_diff,
    ws_extent@ymin - extent_y_diff,
    ws_extent@ymax + extent_y_diff
  ) -> map_extent

  extent(
    min(ws_extent@xmin, state_extent@xmin),
    max(ws_extent@xmax, state_extent@xmax),
    min(ws_extent@ymin, state_extent@ymin),
    max(ws_extent@ymax, state_extent@ymax)
  ) -> state_ws_extent


  # create map visuals
  tmap_options(max.categories = 132)
  tmap_arrange(
    # 3. Power plant map
    # tm_shape(watersheds_city_trans, bbox = map_extent) +
    #   tm_borders() + tm_fill("lightgrey") +
    #   tm_shape(power_plants_city) +
    #   tm_bubbles(shape = "Cooling Technology",
    #              col = "Power Plant Type", size = 0.5) +
    #   tm_layout(frame = TRUE, legend.position = c("left", "top"),
    #             legend.bg.color = "white",
    #             legend.bg.alpha = 0.5,
    #             legend.width = -0.32),
    # 1. map of state with watersheds defined
    tm_shape(state_shape_trans, bbox = state_ws_extent) +
      tm_borders() + tm_fill(col = "white") +
      tm_shape(watersheds_city_trans) +
      tm_fill(col = "dodgerblue") + #tm_borders("darkgrey") +
      tm_shape(city_point) +
      tm_bubbles(col = "hotpink", shape = 22, border.col = "black", size = 2, alpha = 0.7) +
      tm_layout(title = "",
                frame = FALSE),
    # 2. Flood control dams
    # tm_shape(watersheds_city_trans, bbox = map_extent) +
    #   tm_borders() + tm_fill("lightgrey") +
    #   tm_shape(flood_control_dams) +
    #   tm_bubbles(size = 0.5) +
    #   tm_layout(frame = TRUE, legend.position = c("left", "top"),
    #             legend.bg.color = "white",
    #             legend.bg.alpha = 0.5,
    #             legend.width = -0.32),

    # 4. cropcover raster map
    tm_shape(cropcover_city_watersheds_agg, bbox = map_extent) +
      tm_raster(palette = crop_palette,
                legend.show = F) +
      tm_shape(watersheds_city_trans, bbox = map_extent) +
      tm_borders(col = "black") +
      tm_shape(subset(power_plants_city, `Power Plant Type` %in% c("Hydropower", "Thermal"))) +
      tm_bubbles(col = "Power Plant Type", size = 0.5,
                 palette = c("yellow", "red")) +
      tm_shape(flood_control_dams[watersheds_city, ]) +
      tm_bubbles(size = 0.2, col = "lightblue", shape = 25) +
      tm_layout(frame = FALSE, legend.position = c("left", "top"),
                legend.bg.color = "white",
                legend.bg.alpha = 0.5,
                legend.width = -0.32,
                legend.show = FALSE),
    ncol = 2
  )

}
