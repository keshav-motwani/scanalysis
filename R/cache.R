#' Cache and retrieve intermediate steps if path exists
#'
#' @param file_path file path at which to cache
#' @param fn function whose result to cache
#'
#' @importFrom readr write_csv read_csv
#' @importFrom qs qsave qread
#' @importFrom parallel detectCores
#'
#' @return cached object
#' @export
#'
#' @examples
#' NULL
cache = function(file_path, fn) {

  if (grepl(".csv", file_path)) {
    write = readr::write_csv
    read = readr::read_csv
  } else {
    write = function(object, path)
      qs::qsave(object,  path, nthreads = detectCores())
    read = function(path)
      qs::qread(path, nthreads = detectCores())
  }

  if (!file.exists(file_path)) {
    result = fn()
    write(result, file_path)
  } else {
    result = read(file_path)
  }

  return(result)

}