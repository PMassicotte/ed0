#' Compute in-air downward irradiance (Ed0)
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
#' @return Numeric vector containing downward irradiance spectra from 290 nm to 700 nm at 5 nm resolution.
#' @export
#'
#' @examples
#' ed0(100, 12, 67.47973, -63.78953, 3, 330, 1, 0.05)
ed0 <- function(yday, hour, lat, lon, tcl, o3, cf, albedo = 0.05) {

  check_args(yday, hour, lat, lon, tcl, o3, cf, albedo)

  ed0_(yday, hour, lat, lon, tcl, o3, cf, albedo)

}
