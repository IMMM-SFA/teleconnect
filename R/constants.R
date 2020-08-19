# projection
proj4_string <- "+proj=longlat +datum=WGS84 +no_defs"

# non land CDL classes
non_land_cdl_classes <- c(0, 81, 83, 88, 111)

developed_values  <- c(82,121,122,123,124)

non_devcrop_class <- c(0,63:65,81,83:92,111,112,131:195)

hydro_variable <- "Normal_Storage"

AF_to_BCM <- 1233.48e-9

Mgallons_to_BCM <- 0.00378541 * 1e-3  # Million gallons to billion cubic meters

m2_to_km2 <- 1000000

m_to_km <- 1000

mm_to_m <- 0.001

day_to_sec <- 86400

missing_watersheds <- c(3589, 3850, 3678, 3582, 3615, 4226, 4227, 2336, 3448, 2162, 3474) # for runoff

avg_wateruse_ltr_per_day <- 200

MGD_to_m3sec <- 0.0438

BCMyr_to_m3sec <- 31.71

BCMmonth_to_m3sec <- 373.36

ltrday_to_m3sec <- 1.15741e-8
