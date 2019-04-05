#' Compute in-air downward PAR
#'
#' @param yday Numeric. Day of the year.
#' @param hour Numeric. Hour of the day.
#' @param lat Numeric. Latitude.
#' @param lon Numeric. Longitude.
#' @param tcl Numeric. Cloud optical thickness.
#' @param o3 Numeric. Ozone.
#' @param cf Numeric. Cloud fraction.
#' @param albedo Numeric. Surface albedo.
#'
#' @return Numeric. PAR calculated bwetween 400 and 700 nm expressed in umol photons m-2 s-1.
#' @export
#'
#' @examples
#' edpar(100, 12, 67.47973, -63.78953, 3, 330, 1, 0.05)
edpar <- function(yday, hour, lat, lon, tcl, o3, cf, albedo = 0.05) {

  check_args(yday, hour, lat, lon, tcl, o3, cf, albedo)

  lut_file <- system.file("extdata", "Ed0moins_LUT_5nm_v2.dat", package = "ed0")

  edpar_(yday, hour, lat, lon, tcl, o3, cf, albedo, lut_file)

}
