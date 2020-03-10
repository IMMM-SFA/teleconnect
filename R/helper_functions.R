#' suppress_intersect
#' @param method tool used to intersect 2 objects and remove unnecessary warnings
#' @return intersection of 2 objects
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @importFrom sf st_intersection
#' @export
suppress_intersect <- function(object1, object2) {

  suppressWarnings(suppressMessages(st_intersection(object1, object2))) -> object3

  return(object3)
}


#' get_startpoints
#' @param method retrieve startpoints of all lines in a file
#' @return start points of all lines
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @importFrom lwgeom st_startpoint
#' @importFrom sf st_as_sf
#' @export
get_startpoints <- function(transfers) {

  sup(st_startpoint(transfers)) %>%
    as.data.frame() %>% st_as_sf() -> transfer_start

  return(transfer_start)
}


#' get_endpoints
#' @param method retrieve endpoints of all lines in a file
#' @return endpoints of all lines
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @importFrom lwgeom st_endpoint
#' @importFrom sf st_as_sf
#' @export
get_endpoints <- function(transfers) {

  sup(st_endpoint(transfers)) %>%
    as.data.frame() %>% st_as_sf() -> transfer_end

  return(transfer_end)
}



#' sup
#' @param method shorten suppressMessages
#' @return action
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)

#' @export
sup <- function(sup) {

  suppressWarnings(suppressMessages(sup)) -> short_sup

  return(short_sup)
}
