check_args <- function(yday, hour, lat, lon, tcl, o3, cf, albedo) {

  # Everything should be numeric
  assertthat::assert_that(
    is.numeric(yday),
    is.numeric(hour),
    is.numeric(lat),
    is.numeric(lon),
    is.numeric(tcl),
    is.numeric(o3),
    is.numeric(cf),
    is.numeric(albedo)
  )

  # Check bounds
  assertthat::assert_that(
    all(between(yday, 0, 365)),
    all(between(hour, 0, 24)),
    all(between(lat, 45, 90)),
    all(between(lon, -180, 180)),
    all(between(tcl, 0, 64)),
    all(between(o3, 100, 550)),
    all(between(cf, 0, 1)),
    all(between(albedo, 0.05, 0.15))
  )

  # Everything should be of same length (should be changed)
  assertthat::assert_that(
    all(
      length(yday),
      length(hour),
      length(lat),
      length(lon),
      length(tcl),
      length(o3),
      length(cf),
      length(albedo)
    )
  )
}

between <- function(x, a, b) {
  x <= b & a >= a
}
