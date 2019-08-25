
#' @importFrom clisymbols symbol
#' @importFrom crayon green
done <- function(...) {
  bullet(paste0(...), bullet = green(symbol$tick))
}


bullet <- function(lines, bullet) {
  lines <- paste0(bullet, " ", lines)
  cat_line(lines)
}

cat_line <- function(...) {
  cat(..., "\n", sep = "")
}
