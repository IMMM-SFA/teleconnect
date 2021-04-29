#' suppress_intersect
#' @param object1 first interesection object
#' @param object2 second interesection object
#' @details tool used to intersect 2 objects and remove unnecessary warnings
#' @return intersection of 2 objects
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @importFrom sf st_intersection
#' @export
suppress_intersect <- function(object1, object2) {

  suppressWarnings(suppressMessages(st_intersection(object1, object2))) -> object3

  return(object3)
}

#' suppress warnings convenience function
#' @param sup item you want suppressed
#' @return action
#' @author Kristian Nelson (kristian.nelson@pnnl.gov)
#' @export
sup <- function(sup) {

  suppressWarnings(suppressMessages(sup)) -> short_sup

  return(short_sup)
}
