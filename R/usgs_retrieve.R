#' update_USGS_data
#' @details download latest USGS gage readings to update flows
#' @param UWB_USGS_Gages_path full path to WB_USGS_Gages.shp
#' @importFrom dplyr filter one_of mutate as_tibble select group_by summarise tibble
#' @importFrom purrr map pmap_dfr
#' @importFrom tidyr gather spread complete
#' @importFrom dataRetrieval readNWISdv
#' @importFrom lubridate month year
#' @author Sean Turner (sean.turner@pnnl.gov)
update_USGS_data <- function(UWB_USGS_Gages_path){

import_shapefile(UWB_USGS_Gages_path) ->
  watersheds

watersheds %>% as_tibble() %>%
  select(-geometry, -DVSN_Name) %>%
  pmap_dfr(function(DVSN_ID, USGS_1, USGS_2, USGS_3,
                USGS_4, USGS_5, USGS_6,
                USGS_7, USGS_8){

    c(USGS_1, USGS_2, USGS_3,
      USGS_4, USGS_5, USGS_6,
      USGS_7, USGS_8) %>%  .[!is.na(.)] ->
      USGS_gage_ids

    USGS_gage_ids %>%
      map(function(id){
        ifelse(nchar(id) == 7, paste0("0", id), as.character(id))
      }) %>% unlist() -> USGS_gage_ids_fixed

    message(USGS_gage_ids_fixed)

    readNWISdv(siteNumbers = USGS_gage_ids_fixed,
               startDate = "1980-01-01",
               endDate = "2019-12-31",
               parameterCd = "00060"
               ) %>%
      as_tibble() -> gage_data

    if(nrow(gage_data) == 0){
      return(
        tibble()
        )
    }else{
      return(
        gage_data %>%
          select(site_no, date = Date, flow_cfs = X_00060_00003) %>%
          group_by(date) %>% summarise(flow_cfs = sum(flow_cfs)) %>%
          mutate(month = month(date, label = T),
                 year = year(date)) %>%
          group_by(year, month) %>% summarise(flow_cfs = mean(flow_cfs)) %>%
          ungroup() %>%
          complete(expand.grid(year = 1980:2019, month = month.abb)) %>%
          mutate(DVSN = DVSN_ID)
      )
    }

  }) -> total_flows

total_flows %>%
  mutate(flow_cfs = round(flow_cfs, 4)) %>% spread(DVSN, flow_cfs) %>%
  write_csv(
    paste0(system.file("extdata", package = "teleconnect"),
           "/teleconnect_flows_USGS_cfs.csv")
  )

}


