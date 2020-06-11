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
#' @importFrom foreign read.dbf
#' @export
#'
plot_watershed <- function(data_dir,
                           watersheds_file_path = "water/CWM_v2_2/World_Watershed8.shp",
                           crop_file_path = "land/2016_90m_cdls/cdl_lowres_usa.img",
                           crop_attributes = "land/2016_90m_cdls/cdl_lowres_usa.img.vat.dbf",
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

  raster::crs(cropcover_USA) <- raster::crs(watersheds_city)

  # read crop table
  foreign::read.dbf(paste0(data_dir, crop_attributes)) ->crop_cover_levels

  crop_cover_levels %>%
    filter(Class_Name != "") %>%
    mutate(color = rgb(Red, Green, Blue, maxColorValue = 255)) %>%
    dplyr::select(Value, cropclass = Class_Name, color) %>%
    mutate(color = if_else(Value == 0, "#ffffff", color)) %>%
    mutate(cropclass = as.character(cropclass)) -> cropclasses

  cropclasses %>%
    dplyr::rename(ID = Value) %>%
    select(ID, cropclass) %>%
    mutate(cropclass = if_else(ID == 131, "Barren_", cropclass),
           cropclass = if_else(ID == 152, "Shrubland_", cropclass)) -> crop_tbl
  crop_tbl$ID <- as.integer(crop_tbl$ID)
  crop_tbl$cropclass <- as.factor(crop_tbl$cropclass)
  levels(cropcover_USA) <- crop_tbl

  # make color palette
  setNames(cropclasses$color, cropclasses$cropclass) -> crop_palette

  # mask to required city
  cropcover_USA %>%
    mask_raster_to_polygon(watersheds_city) ->
    cropcover_city_watersheds_agg
  cropcover_city_watersheds_agg@legend@colortable <- logical(0)

  rm(cropcover_USA)

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

  subset(power_plants_city, `Power Plant Type` %in% c("Hydropower", "Thermal")) -> power_plant_list
  if(nrow(power_plant_list) == 0){
    if(nrow(flood_control_dams) == 0){
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
          tm_layout(title = paste0(city),
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
          tm_layout(frame = FALSE, legend.position = c("left", "top"),
                    legend.bg.color = "white",
                    legend.bg.alpha = 0.5,
                    legend.width = -0.32,
                    legend.show = FALSE),
        ncol = 2
      )
    }else{

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
          tm_layout(title = paste0(city),
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
  }else{

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
        tm_layout(title = paste0(city),
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
}


#' plot_watershed_ggplot
#'
#' @details plots the watersheds and relevant features for a chosen city
#' @param data_dir root directory for the spatial data ("/pic/projects/im3/teleconnections/data/")
#' @param watersheds_file_path path of watersheds shapefile within data_dir
#' @param landcover_file_path path of NLCD landcover data file
#' @param powerplants_file_path path of power plants data file
#' @param crop_file_path path of CDL img file
#' @param conus_file_path path of CONUS state boundaries shape file
#' @param city character. City to be analysed in the format "City | STATE" (e.g., "Phoenix | AZ")
#' @importFrom tmaptools aggregate_map
#' @importFrom dplyr filter mutate if_else right_join
#' @importFrom tibble rowid_to_column
#' @importFrom raster levels extent rasterToPoints
#' @importFrom foreign read.dbf
#' @importFrom vroom vroom
#' @import ggplot2
#'
#' @export
#'
plot_watershed_ggplot <- function(data_dir,
                           watersheds_file_path = "water/CWM_v2_2/World_Watershed8.shp",
                           crop_file_path = "land/2016_90m_cdls/cdl_lowres_usa.img",
                           crop_attributes = "land/2016_90m_cdls/cdl_lowres_usa.img.vat.dbf",
                           powerplants_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
                           conus_file_path = "misc/states_21basic/states.shp",
                           dams_file_path = "water/nabd_fish_barriers_2012/nabd_fish_barriers_2012.shp",
                           nlud_file_path = "land/usa_nlud_raster.tif",
                           citypoint_file_path = "water/CWM_v2_2/City_Centroid.shp",
                           city){

  all_cities <- get_cities()[["city_state"]]

  get_teleconnect_table() -> all_data
  all_data[complete.cases(all_data), ] -> all_data
  all_data$city -> plotted_cities

  # throw error if any chosen city lies outside available cities
  if(!city %in% plotted_cities) {
    stop(paste0(paste(city), " has no watershed and cannot be plotted!"))
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

  raster::crs(cropcover_USA) <- raster::crs(watersheds_city)

  # read crop table
  foreign::read.dbf(paste0(data_dir, crop_attributes)) ->crop_cover_levels

  crop_cover_levels %>%
    filter(Class_Name != "") %>%
    mutate(color = rgb(Red, Green, Blue, maxColorValue = 255)) %>%
    dplyr::select(Value, cropclass = Class_Name, color) %>%
    mutate(color = if_else(Value == 0, "#ffffff", color)) %>%
    mutate(cropclass = as.character(cropclass)) -> cropclasses

  cropclasses %>%
    dplyr::rename(ID = Value) %>%
    select(ID, cropclass) %>%
    mutate(cropclass = if_else(ID == 131, "Barren_", cropclass),
           cropclass = if_else(ID == 152, "Shrubland_", cropclass)) -> crop_tbl
  crop_tbl$ID <- as.integer(crop_tbl$ID)
  crop_tbl$cropclass <- as.factor(crop_tbl$cropclass)
  levels(cropcover_USA) <- crop_tbl

  # make color palette
  setNames(cropclasses$color, cropclasses$cropclass) -> crop_palette

  # mask to required city
  cropcover_USA %>%
    mask_raster_to_polygon(watersheds_city) ->
    cropcover_city_watersheds_agg
  cropcover_city_watersheds_agg@legend@colortable <- logical(0)

  rm(cropcover_USA)

  # get the coordinate system of the crop raster and transform watershed
  r_crs <- st_crs(projection(cropcover_city_watersheds_agg))
  watersheds_city_trans <- st_as_sf(watersheds_city) %>% st_transform(crs = r_crs)
  state_shape_trans <- st_as_sf(state_shape) %>% st_transform(crs = r_crs)

#--------------CONVERT RASTER--------------------------------------------------------------------

  rasterToPoints(cropcover_city_watersheds_agg) %>%
    data.frame() -> df
  colnames(df) <- c("Longitude", "Latitude", "Value")
  left_join(df, cropclasses, by = "Value") -> df2

  flood_control_dams %>% st_as_sf() -> flood_control_sf
  watersheds_city %>% st_as_sf() -> watershed_sf
  flood_control_sf[watershed_sf, ] -> select_dams

  subset(power_plants_city, `Power Plant Type` %in% c("Hydropower", "Thermal")) %>%
    as.data.frame() %>%
    mutate(color2 = if_else(Power.Plant.Type == "Hydropower", "yellow1", "red")) -> plants

#----------------GET TELE DATA--------------------------------------------------------------------

  all_data %>% filter(city == !!city) -> city_data

  city_data %>% pull(yield_BCM) -> yield_point
  city_data %>% pull(irr_cons_BCM) -> cons_point
  cons_point / yield_point -> cons_per_yield


  city_data %>% pull(n_other_cities) -> number_cities

  city_data %>% pull(n_hydro_plants) -> hydro
  city_data %>% pull(n_thermal_plants) -> thermal
  hydro + thermal -> total_plants

  city_data %>% pull(cropland_fraction) -> crop_point

#------------------GRAPH PLOTS----------------------------------------------------------

  all_data %>%
    ggplot(aes(n_other_cities)) +
    stat_ecdf(pad = F) -> cdf_plot
  cdf_plot +
    geom_point(data = layer_data(cdf_plot) %>% filter(x == number_cities),
               aes(x, y), col = "red", size = 5) +
    xlab("Other Cities") +
    ylab("") +
    theme_classic() -> plot3

  all_data %>%
    ggplot(aes(irr_cons_BCM / yield_BCM)) +
    stat_ecdf(pad = F) -> cdf_plot
  cdf_plot +
    geom_point(data = layer_data(cdf_plot) %>% filter(x == cons_per_yield),
               aes(x, y), col = "red", size = 5) +
    xlab("IrrCons/Yield") +
    ylab("") +
    theme_classic() -> plot4

  all_data %>%
    ggplot(aes(n_hydro_plants + n_thermal_plants)) +
    stat_ecdf(pad = F) -> cdf_plot
  cdf_plot +
    geom_point(data = layer_data(cdf_plot) %>% filter(x == total_plants),
               aes(x, y), col = "red", size = 5) +
    xlab("Total Plants") +
    ylab("") +
    theme_classic() -> plot5

  all_data %>%
    ggplot(aes(cropland_fraction)) +
    stat_ecdf(pad = F) -> cdf_plot
  cdf_plot +
    geom_point(data = layer_data(cdf_plot) %>% filter(x == crop_point),
               aes(x, y), col = "red", size = 5) +
    xlab("Crop Fraction") +
    ylab("") +
    theme_classic() -> plot6

  cowplot::plot_grid(plot3,plot5) -> graph_plot1
  cowplot::plot_grid(plot4,plot6) -> graph_plot2

  ggplot(state_shape_trans)+
    geom_sf(fill = "white")+
    geom_sf(data = watersheds_city_trans, fill = "dodgerblue3", colour = "darkgrey") +
    geom_sf(data = city_point,fill = "hotpink", colour = "black", size = 7, alpha = 0.7, shape = 22) +
    theme_void() +
    ggtitle(paste0(city)) -> plot1

#---------------PLOT ALL TOGETHER----------------------------------------------------------------------

  subset(power_plants_city, `Power Plant Type` %in% c("Hydropower", "Thermal")) -> power_plant_list
  if(nrow(power_plant_list) == 0){
    if(nrow(flood_control_dams) == 0){

      ggplot(df2) +
        geom_raster(aes(x= Longitude, y=Latitude, fill=color)) +
        scale_fill_identity() +
        geom_sf(data = watersheds_city_trans, fill = "white",colour = "black", alpha = 0) +
        theme_void() -> plot2

    }else{

      ggplot(df2) +
        geom_raster(aes(x= Longitude, y=Latitude, fill=color)) +
        scale_fill_identity() +
        geom_sf(data = watersheds_city_trans, fill = "white",colour = "black", alpha = 0) +
        geom_sf(data = select_dams, size = 2, fill = "lightblue", colour = "black", shape = 25) +
        theme_void() -> plot2
}
  }else{

    ggplot(df2) +
      geom_raster(aes(x= Longitude, y=Latitude, fill=color)) +
      scale_fill_identity() +
      geom_sf(data = watersheds_city_trans, fill = "white",colour = "black", alpha = 0) +
      geom_point(data = plants, aes(x = lon.1, y = lat.1, color = Power.Plant.Type),size=3) +
      scale_color_manual(values = c("Hydropower" = "yellow", "Thermal" = "red")) +
      geom_sf(data = select_dams, size = 2, fill = "lightblue", colour = "black", shape = 25) +
      theme_void() -> plot2
  }
  cowplot::plot_grid(plot1,plot2,graph_plot1,graph_plot2)
}
