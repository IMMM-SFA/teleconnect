# Run the function below to download the `gamut` datasets

#' get_gamut_datasets
#'
#' @param download_path path to the downloaded file
#' @param unzip_path path to the unzipped file
#' @details function to download the gamut input datasets
#' @export
get_gamut_datasets <- function(download_path, unzip_path){
  #Download `gamut` zip file
  download.file("https://zenodo.org/record/4662993/files/Geospatial_Input_Datasets.zip?download=1",download_path)
  #Unzip folder to final data directory
  unzip(download_path,exdir = unzip_path)

}
