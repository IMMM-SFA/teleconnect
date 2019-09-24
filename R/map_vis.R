#' plot_watershed
#'
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param watersheds_file_path path of watersheds shapefile within data_dir
#' @param landcover_file_path path of NLCD landcover data file
#' @param powerplants_file_path path of power plants data file
#' @param city character. City to be analysed in the format "City | STATE" (e.g., "Phoenix | AZ")
#' @details plots the watersheds and relevant features for a chosen city
#' @import tmap
#' @importFrom tmaptools aggregate_map
#' @importFrom dplyr filter mutate if_else right_join
#' @importFrom tibble rowid_to_column
#' @importFrom raster levels extent
#' @export
#'
plot_watershed <- function(data_dir,
                           watersheds_file_path = "water/CWM_v2_2/World_Watershed8.shp",
                           landcover_file_path = "land/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img",
                           crop_file_path = "land/2016_30m_cdls/2016_30m_cdls.img",
                           powerplants_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
                           conus_file_path = "misc/states_21basic/states.shp",
                           city){

  all_cities <- get_cities()[["city_state"]]

  # throw error if any chosen city lies outside available cities
  if(!city %in% all_cities) {
    stop(paste0(paste(city), " is not part of 'teleconnect'!"))
  }

  # get city watershed_ids
  get_cities() %>%
    filter(city_state == !!city) ->
    city_watershed_mapping

  # read shapefiles for watersheds
  import_shapefile(paste0(data_dir, watersheds_file_path),
                   method = "rgdal") %>%
    # subset for desired watersheds
    subset(DVSN_ID %in% city_watershed_mapping$DVSN_ID) ->
    watersheds_city

  import_shapefile("../spatial data/water/CWM_v2_2/City_Centroid.shp",
                   method = "rgdal") %>%
    subset(grepl("USA", ISOURBID)) %>%
    subset(NAME == city)

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
    mutate(`Power Plant Type` = as.factor(`Power Plant Type`),
           `Cooling Technology` = as.factor(cooling_tech))

  # read landcover raster for US
  import_raster(paste0(data_dir, landcover_file_path)) ->
    landcover_USA

  # read crop raster for US
  import_raster(paste0(data_dir, crop_file_path)) ->
    cropcover_USA

  # get land classes and color scheme
  levels(landcover_USA)[[1]] %>%
    filter(NLCD.Land.Cover.Class != "") %>%
    mutate(color = rgb(Red, Green, Blue, maxColorValue = 255)) %>%
    dplyr::select(ID, landclass = NLCD.Land.Cover.Class, color) %>%
    mutate(color = if_else(ID == 0, "#ffffff", color)) ->
    landclasses

  # get crop classes and color scheme
  levels(cropcover_USA)[[1]] %>%
    filter(Class_Names != "") %>%
    mutate(color = rgb(Red, Green, Blue, maxColorValue = 255)) %>%
    dplyr::select(ID, cropclass = Class_Names, color) %>%
    mutate(color = if_else(ID == 0, "#ffffff", color)) %>%
    mutate(cropclass = as.character(cropclass)) ->
    cropclasses

  # make color palette
  setNames(landclasses$color, landclasses$landclass) -> land_palette
  setNames(cropclasses$color, cropclasses$cropclass) -> crop_palette

  # mask to required city
  landcover_USA %>%
    mask_raster_to_polygon(watersheds_city) ->
    landcover_city_watersheds
  rm(landcover_USA)

  cropcover_USA %>%
    mask_raster_to_polygon(watersheds_city) ->
    cropcover_city_watersheds
  rm(cropcover_USA)

  # get raster attribute table
  rat_lc <- levels(landcover_city_watersheds)[[1]]
  rat_lc <- right_join(rat_lc, landclasses, by = "ID") %>%
    dplyr::select(ID, landclass)
  levels(landcover_city_watersheds) <- rat_lc

  rat_cc <- levels(cropcover_city_watersheds)[[1]]
  rat_cc <- right_join(rat_cc, cropclasses, by = "ID") %>%
    dplyr::select(ID, cropclass) %>%
    mutate(cropclass = if_else(ID == 131, "Barren_", cropclass),
           cropclass = if_else(ID == 152, "Shrubland_", cropclass)) %>%
    mutate(cropclass = as.factor(cropclass))
  levels(cropcover_city_watersheds) <- rat_cc

  # aggregate raster to increase plot speed (and also convert to format usable with tmap)
  landcover_city_watersheds %>%
    aggregate_map(fact = 4, agg.fun = "modal") ->
    landcover_city_watersheds_agg

  # aggregate raster to increase plot speed (and also convert to format usable with tmap)
  cropcover_city_watersheds %>%
    aggregate_map(fact = 4, agg.fun = "modal") ->
    cropcover_city_watersheds_agg

  # correct raster attribute table in aggregated raster
  levels(landcover_city_watersheds_agg) %>% .[[1]] %>%
    mutate(ID = rat_lc$ID) -> rat_lc_agg
  levels(landcover_city_watersheds_agg) <- rat_lc_agg

  levels(cropcover_city_watersheds_agg) %>% .[[1]] %>%
    mutate(ID = rat_cc$ID) -> rat_cc_agg
  levels(cropcover_city_watersheds_agg) <- rat_cc_agg

  # get the coordinate system of the landclass raster and transform watershed
  r_crs <- st_crs(projection(landcover_city_watersheds_agg))
  watersheds_city_trans <- st_as_sf(watersheds_city) %>% st_transform(crs = r_crs)
  state_shape_trans <- st_as_sf(state_shape) %>% st_transform(crs = r_crs)

  #power_plants_city_trans <- st_as_sf(power_plants_city) %>% st_transform(crs = r_crs)

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

  # plot NCLD raster

  tmap_options(max.categories = 132)

  tmap_arrange(
    tm_shape(state_shape_trans, bbox = state_ws_extent) +
      tm_borders() + tm_fill(col = "white") +
      tm_shape(watersheds_city_trans) +
      tm_fill(col = "lightgrey") + tm_borders("black") +
      tm_layout(title = paste0("Water supply catchments \n", city),
                frame = FALSE),
    tm_shape(landcover_city_watersheds_agg, bbox = map_extent) +
      tm_raster(style = "cat",
                palette = land_palette,
                title = "NLCD Land Class") +
        tm_layout(frame = TRUE,
                  legend.frame = FALSE,
                  legend.bg.color = "white",
                  legend.bg.alpha = 0.5,
                  legend.width = -0.32),
    tm_shape(watersheds_city_trans, bbox = map_extent) +
      tm_borders() + tm_fill("lightgrey") +
      tm_shape(power_plants_USA) +
      tm_bubbles(shape = "Cooling Technology",
                 col = "Power Plant Type", size = 0.5) +
      tm_layout(frame = TRUE, legend.position = c("left", "top"),
                legend.bg.color = "white",
                legend.bg.alpha = 0.5,
                legend.width = -0.32),
    tm_shape(cropcover_city_watersheds_agg, bbox = map_extent) +
      tm_raster(palette = crop_palette,
                legend.show = F) +
      tm_layout(frame = TRUE,
                title = "Crop Class"),
    ncol = 2
  )

}
