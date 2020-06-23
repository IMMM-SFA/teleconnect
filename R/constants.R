# projection
proj4_string <- "+proj=longlat +datum=WGS84 +no_defs"

# non land CDL classes
non_land_cdl_classes <- c(0, 81, 83, 88, 111)

non_devcrop_class <- c(0,63:65,81,83:92,111,112,131:195)

hydro_variable <- "Normal_Storage"

AF_to_BCM <- 1233.48e-9

Mgallons_to_BCM <- 0.00378541 * 1e-3  # Million gallons to billion cubic meters

m_to_km <- 1000

missing_watersheds <- c(3589, 3850, 3678, 3582, 3615, 4226, 4227, 2336, 3448, 2162, 3474) # for runoff

avg_wateruse_ltr <- 333
