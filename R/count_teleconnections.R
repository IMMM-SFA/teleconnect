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
                                            crop_file_path = "land/2016_30m_cdls/2016_30m_cdls.img",
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

  # read croptype raster for US
  import_raster(paste0(data_dir, crop_file_path)) ->
    cropcover_USA

      # START GCAM CROP CLASS BINDING

      # Load in the GCAM classes CSV.
      gcam_csv <- get_crop_mapping()

      # Delete the rows that have crop types not present in CDL classes.
      gcam_csv <- gcam_csv[-c(70, 100, 131, 45, 104, 182, 181),]

      # Extract only the GCAM Class column
      as.data.frame(gcam_csv[,c(3)]) %>%
        mutate(GCAM_ID = group_indices_(gcam_csv, .dots=c("GCAM_commodity"))) ->
        gcam_classes

      # Shrink the Column so that there is only one of each Class
      gcam_class <- distinct(gcam_classes)

      # Obtain the attribute table for the raster.
      levels(cropcover_USA)[[1]] %>%
        filter(Class_Names != "") %>%
        dplyr::select(ID, CDL_Class = Class_Names) %>%
        mutate(CDL_Class = as.character(CDL_Class))->
            cdl_classes

      # This HUGE chunk assigns the GCAM IDs to matching CDL Class. Most likely a more efficient way of doing this.
      # Since there were some differences in classes, I did the best I could to assign each class properly. All Double Crops are under MiscCrop.
      cond1 <- cdl_classes$ID %in% c(1,12,13,225)
        cdl_classes$GCAM_ID[cond1] <- 1
      cond2 <- cdl_classes$ID %in% c(2,32,36)
        cdl_classes$GCAM_ID[cond2] <- 2
      cond3 <- cdl_classes$ID %in% c(59,60)
        cdl_classes$GCAM_ID[cond3] <- 3
      cond4 <- cdl_classes$ID %in% c(37,58,224)
        cdl_classes$GCAM_ID[cond4] <- 4
      cond5 <- cdl_classes$ID %in% c(10,11,14,42,44,47:57,61,66:72,74:77,204,206:210,212:223,227,229:254,226)
        cdl_classes$GCAM_ID[cond5] <- 5
      cond6 <- cdl_classes$ID %in% c(5,6,26,31,33:35,38,211)
        cdl_classes$GCAM_ID[cond6] <- 6
      cond7 <- cdl_classes$ID %in% c(4,21,25,27:29,39,205)
        cdl_classes$GCAM_ID[cond7] <- 7
      cond8 <- cdl_classes$ID %in% c(3)
        cdl_classes$GCAM_ID[cond8] <- 8
      cond9 <- cdl_classes$ID %in% c(43,46)
        cdl_classes$GCAM_ID[cond9] <- 9
      cond10 <- cdl_classes$ID %in% c(41,45)
        cdl_classes$GCAM_ID[cond10] <- 10
      cond11 <- cdl_classes$ID %in% c(22:24,30)
        cdl_classes$GCAM_ID[cond11] <- 11
      cond12 <- cdl_classes$ID %in% c(82,121:124)
        cdl_classes$GCAM_ID[cond12] <- 12
      cond13 <- cdl_classes$ID %in% c(87,190,195)
        cdl_classes$GCAM_ID[cond13] <- 13
      cond14 <- cdl_classes$ID %in% c(92)
        cdl_classes$GCAM_ID[cond14] <- 14
      cond15 <- cdl_classes$ID %in% c(141:143,63)
        cdl_classes$GCAM_ID[cond15] <- 15
      cond16 <- cdl_classes$ID %in% c(64,152:176)
        cdl_classes$GCAM_ID[cond16] <- 16
      cond17 <- cdl_classes$ID %in% c(112)
        cdl_classes$GCAM_ID[cond17] <- 17
      cond18 <- cdl_classes$ID %in% c(131,65)
        cdl_classes$GCAM_ID[cond18] <- 18
      cond19 <- cdl_classes$ID %in% c(81,88,83,111,0)
        cdl_classes$GCAM_ID[cond19] <- 0


      # Once all of the classes are assigned, this next part merges the GCAM Class name with the ID.
      crop_reclass_table<- merge(cdl_classes, gcam_class, by = "GCAM_ID", all.x = TRUE)


      # Since there are some Land Cover classes in this data set, not all matched up with the GCAM Classes and were NA. If they are NA, this code fills it in with the CDL Class
      crop_reclass_table$GCAM_commodity <- if_else(is.na(crop_reclass_table$GCAM_commodity),
                                                    paste(crop_reclass_table$CDL_Class),
                                                    paste(crop_reclass_table$GCAM_commodity))

      # Rename the Columns for simplicity.
      crop_reclass_table <- rename(crop_reclass_table, CDL_ID = ID)

      crop_reclass_table <- rename(crop_reclass_table, GCAM_Class = GCAM_commodity)

      # This can get removed, but this clumps all of the Land covers into more generalized classes, and makes background, water, and undefined all NA.
      # Again, there could be a more efficient way of doing this.
      cond22 <- crop_reclass_table$GCAM_ID %in% 15
        crop_reclass_table$GCAM_Class[cond22] <- "Forest"
      cond23 <- crop_reclass_table$GCAM_ID %in% 12
        crop_reclass_table$GCAM_Class[cond23] <- "Developed"
      cond24 <- crop_reclass_table$GCAM_ID %in% 13
        crop_reclass_table$GCAM_Class[cond24] <- "Wetlands"
      cond25 <- crop_reclass_table$GCAM_ID %in% 16
        crop_reclass_table$GCAM_Class[cond25] <- "Rangeland"
      cond25 <- crop_reclass_table$GCAM_ID %in% 0
        crop_reclass_table$GCAM_Class[cond25] <- "N/A"

      # END GCAM CROP CLASS BINDING


  # map through all cities, computing teleconnections
  cities %>%
    map_dfr(function(city){
      filter(city_watershed_mapping, city_state == !!city) %>%
        .[["DVSN_ID"]] -> city_intake_ids

      # subset watersheds shapefile for target city
      watersheds %>%
        subset(DVSN_ID %in% city_intake_ids) ->
        watersheds_city

      # catch cases with only groundwater points (i.e., no watershed polygons)
      if(nrow(watersheds_city) == 0){
        done(city)
        return(
          tibble(city = !! city,
                 n_watersheds = 0,
                 n_hydro = 0,
                 n_thermal = 0,
                 n_landclasses = 0,
                 n_cropcover = 0)
        )
      }else{

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


        # TELECONNECTION - NUMBER OF CROP TYPES BASED ON GCAM CLASSES

        get_raster_val_classes(cropcover_USA, watersheds_city) -> cropcover_count

            # Create into data frame, merge GCAM IDs, keep only the column with the GCAM Ids, filter out only unique values, count values.
            as.data.frame(cropcover_count) -> cropcover_count_df

            cropcover_count_df <- rename(cropcover_count_df, CDL_ID = cropcover_count)

            crop_merge<- merge(cropcover_count_df, crop_reclass_table, by = "CDL_ID", all.x = TRUE)

            remove_water <- 0

            crop_merge2 <- subset(crop_merge, crop_merge[,2] > remove_water)

            crop_merge3 <- crop_merge2[,-c(1,3,4)]
            crop_list <- unique.numeric_version(crop_merge3)

            tc_n_cropcover <- length(crop_list)

        done(city)

        # output
        return(
          tibble(city = !! city,
                 n_watersheds = tc_n_watersheds,
                 n_hydro = tc_n_hydroplants,
                 n_thermal = tc_n_thermalplants,
                 n_landclasses = tc_n_landclasses,
                 n_cropcover = tc_n_cropcover)
        )
      }

    })
}
