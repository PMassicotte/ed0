#' Compute in-air downward PAR
#'
#' @param yday Numeric. Day of the year.
#' @param hour Numeric. Hour of the day in UTC (0 to 24).
#' @param lat Numeric. Latitude.
#' @param lon Numeric. Longitude.
#' @param tcl Numeric. Cloud optical thickness (0 to 64.0).
#' @param o3 Numeric. Ozone concentration (100 to 550 DU).
#' @param cf Numeric. Cloud fraction (0 to 1.0).
#' @param albedo Numeric. Surface albedo (0 to 1.0).
#' @param lut_type Character. Either "ed0+" for in-air irradiance or "ed0-" for underwater irradiance.
#'
#' @return Numeric. PAR calculated bwetween 400 and 700 nm expressed in umol photons m-2 s-1.
#' @export
#'
#' @examples
#' edpar(100, 12, 67.47973, -63.78953, 3, 330, 1, 0.05, "ed0+")
edpar <- function(yday, hour, lat, lon, tcl, o3, cf, albedo = 0.05, lut_type) {
  check_args(yday, hour, lat, lon, tcl, o3, cf, albedo, lut_type)

  lut_file <-
    ifelse(
      lut_type == "ed0-",
      system.file("extdata", "Ed0moins_LUT_5nm_v2.dat", package = "ed0"),
      system.file("extdata", "Ed0plus_LUT_5nm_v2.dat", package = "ed0")
    )

  edpar_(yday, hour, lat, lon, tcl, o3, cf, albedo, lut_file)
}
