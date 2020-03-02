
# projection
proj4_string <- "+proj=longlat +datum=WGS84 +no_defs"

# non land CDL classes
non_land_cdl_classes <- c(0, 81, 83, 88, 111)

non_devcrop_class <- c(0,63:65,81,83:92,111,112,131:195)

hydro_variable <- "Normal_Storage"

acrefeet_conv <- 1.2335e-6

c(watersheds_file_path = "water/CWM_v2_2/World_Watershed8.shp",
  withdrawal_file_path = "water/CWM_v2_2/Snapped_Withdrawal_Points.shp",
  citypoint_file_path = "water/CWM_v2_2/City_Centroid.shp",
  powerplants_file_path = "water/UCS-EW3-Energy-Water-Database.xlsx",
  crop_file_path = "land/2016_90m_cdls/cdl_lowres_usa.img",
  crop_attribute_path = "land/2016_90m_cdls/cdl_lowres_usa.img.vat.dbf",
  irrigation_file_path = "land/usa_demeter.csv",
  nlud_file_path = "land/usa_nlud_LR.tif",
  hydro_file_path = "energy/EHA_Public_PlantFY2019_GIS_6/ORNL_EHAHydroPlant_PublicFY2019final.xlsx",
  transfers_file_path = "water/USIBTsHUC6_Dickson.shp",
  climate_file_path = "land/kop_climate_classes.tif") -> teleconnect_files
