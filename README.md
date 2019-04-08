
# ed0 [![CRAN status](https://www.r-pkg.org/badges/version/ed0)](https://cran.r-project.org/package=ed0) [![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](http://www.gnu.org/licenses/gpl-2.0.html) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/PMassicotte/ed0?branch=master&svg=true)](https://ci.appveyor.com/project/PMassicotte/ed0) [![Travis build status](https://travis-ci.org/PMassicotte/ed0.svg?branch=master)](https://travis-ci.org/PMassicotte/ed0)

The goal of `ed0` is to compute in-air downward irradiance for the
pan-Arctic region. `ed0` uses a lookup table (LUT) generated using the
radiative transfer model SBDART (Santa Barbara DISORT Atmospheric
Radiative Transfer; [Ricchiazzi, Yang, Gautier, &
Sowle, 1998](https://journals.ametsoc.org/doi/abs/10.1175/1520-0477%281998%29079%3C2101%3ASARATS%3E2.0.CO%3B2)).
The LUT and the original Fortran code to interpolate the values inside
the LUT have been produced by Dr. Simon Bélanger.

## Installation

The development version from GitHub can be installed with:

``` r
# install.packages("devtools")
devtools::install_github("pmassicotte/ed0")
```

## Example

``` r
library(tidyverse)
library(ed0)

df <- tibble(
  wavelength = seq(290, 700, by = 5),
  ed0 = ed0(
    yday = 100,
    hour = 12,
    lat = 67.47973,
    lon = -63.78953,
    tcl = 3,
    o3 = 330,
    cf = 1,
    albedo = 0.05
  )
)

df
#> # A tibble: 83 x 2
#>    wavelength         ed0
#>         <dbl>       <dbl>
#>  1        290 0.000000200
#>  2        295 0.000000533
#>  3        300 0.0000422  
#>  4        305 0.00112    
#>  5        310 0.0107     
#>  6        315 0.0429     
#>  7        320 0.0728     
#>  8        325 0.132      
#>  9        330 0.221      
#> 10        335 0.235      
#> # … with 73 more rows

df %>% 
  ggplot(aes(x = wavelength, y = ed0)) +
  geom_line() +
  geom_point() +
  theme_bw() + 
  xlab("Wavelength (nm)") +
  ylab(bquote(Downward~irradiance~(mu*mol~photons~s^{-1}~m^{-2})))
```

<img src="man/figures/README-example-1.svg" width="100%" />

## Code of conduct

Please note that the \[‘ed0’ project is released with a [Contributor
Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.
