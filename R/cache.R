#' Cache and retrieve intermediate steps if path exists
#'
#' @param file_path File path at which to cache
#' @param fn Function whose reslt to cache
#'
#' @return Cached object
#' @export
#'
#' @examples
#' NULL
cache = function(file_path, fn) {
  if (grepl(".csv", file_path)) {
    write = readr::write_csv
    read = readr::read_csv
  } else {
    write = function(object, path) qs::qsave(object,  path, nthreads = 32)
    read = function(path) qs::qread(path, nthreads = 32)
  }

  if (!file.exists(file_path)) {
    result = fn()
    write(result, file_path)
  } else {
    result = read(file_path)
  }
  return(result)
}