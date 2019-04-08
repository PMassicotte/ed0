#' Compute in-air downward irradiance (Ed0)
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
#' @return Numeric vector containing downward irradiance spectra from 290 nm to 700 nm at 5 nm resolution.
#' @export
#'
#' @examples
#' ed0(100, 12, 67.47973, -63.78953, 3, 330, 1, 0.05, "ed0+")
ed0 <- function(yday, hour, lat, lon, tcl, o3, cf, albedo = 0.05, lut_type = "ed0-") {
  check_args(yday, hour, lat, lon, tcl, o3, cf, albedo, lut_type)

  lut_file <-
    ifelse(
      lut_type == "ed0-",
      system.file("extdata", "Ed0moins_LUT_5nm_v2.dat", package = "ed0"),
      system.file("extdata", "Ed0plus_LUT_5nm_v2.dat", package = "ed0")
    )

  ed0_(yday, hour, lat, lon, tcl, o3, cf, albedo, lut_file)
}
