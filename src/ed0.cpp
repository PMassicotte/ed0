// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <stdio.h>
#include <strings.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "math.h"

#define NVIS 61
#define NBWL 83
#define NO3 10
#define NTAUCLD 8
#define NTHETAS 19
#define NALB 7

double sun_zenithal_angle(int yday, double hour, double minute, double second, double latitude, double longitude) {

    double solar_zenith_angle = 0;

    double d2r = M_PI / 180;
    double r2d = 1 / d2r;

    double d = 23.45 * d2r * sin(d2r * 360 * (284 + yday) / 365);

    double E_qt = 0;

    if (yday <= 106) {
        E_qt = -14.2 * sin(M_PI * (yday + 7) / 111);
    } else {

        if (yday <= 166) {
            E_qt = 4 * sin(M_PI * (yday - 106) / 59);
        } else {

            if (yday <= 246) {
                E_qt = -6.5 * sin(M_PI * (yday - 166) / 80);
            } else {
                E_qt = 16.4 * sin(M_PI * (yday - 247) / 113);
            }
        }
    }

    double T = hour + minute / 60 + second / 3600;

    double T_solar = T + longitude / 15 + E_qt / 60;
    double w = M_PI * (12 - T_solar) / 12;
    double l = latitude * d2r;

    solar_zenith_angle =
        90 - asin(sin(l) * sin(d) + cos(l) * cos(d) * cos(w)) * r2d;

    return (solar_zenith_angle);
}

void get_indice(std::vector<float> vec, float target, int &ii, float &rr, float overflow_max_value) {

    // std::cout << target << "****" << std::endl;

    // std::cout << target << ":" <<  vec[0] << std::endl;

    if (target < vec[0]) {
        ii = 0;
        rr = 0;
    } else if (target >= vec[vec.size() - 1]) {
        ii = vec.size() - 2;
        rr = (overflow_max_value - vec[ii]) / (vec[ii + 1] - vec[ii]);
    }

    else {

        for (unsigned int i = 0; i < vec.size() - 1; i += 1) {

            if (target >= vec[i] && target < vec[i + 1]) {
                ii = i;
            }
        }

        rr = (target - vec[ii]) / (vec[ii + 1] - vec[ii]);
    }
}

void read_ed0moins_lut_(const char *filename, float downward_irradiance_table_as_output[NBWL][NTHETAS][NO3][NTAUCLD][NALB]) {


    std::ifstream infile;
    infile.open(filename);
    float tmp;
    int iteration = 0;
    for (int theta = 0; theta < NTHETAS; theta++) {
        for (int ozone = 0; ozone < NO3; ozone++) {
            for (int taucl = 0; taucl < NTAUCLD; taucl++) {
                for (int albedo = 0; albedo < NALB; albedo++) {
                    for (int wavelength = 0; wavelength < NBWL; wavelength++) {
                        infile >> tmp;
                        downward_irradiance_table_as_output[wavelength][theta][ozone][taucl][albedo] = tmp;
                    }
                }
            }
        }
    }

    // Close file
    infile.close();
}

std::vector<float> interpol_ed0moins(float ed_lut[83][19][10][8][7], double thetas, double ozone, double taucl, double alb) {

    int nwl = 83;

    float ed_tmp4[nwl][2][2][2];
    float ed_tmp3[nwl][2][2];
    float ed_tmp2[nwl][2];
    std::vector<float> ed(nwl);

    // Thetas
    std::vector<float> xthetas;
    for (int i = 0; i <= 90; i += 5) {
        xthetas.push_back(i);
    }

    // Ozone
    std::vector<float> xozone;
    for (int i = 100; i <= 550; i += 50) {
        xozone.push_back(i);
    }

    // Taucl
    std::vector<float> xtaucl;
    xtaucl.push_back(0);
    xtaucl.push_back(1);
    xtaucl.push_back(2);
    xtaucl.push_back(4);
    xtaucl.push_back(8);
    xtaucl.push_back(16);
    xtaucl.push_back(32);
    xtaucl.push_back(64);

    // Albedo
    std::vector<float> xalb;
    for (float i = 0.05; i <= 1; i += 0.15) {
        xalb.push_back(i);
    }

    // Find the index where the requested values are.
    int ithetas, iozone, itaucl, ialb;
    float rthetas, rozone, rtaucl, ralb;

    get_indice(xthetas, thetas, ithetas, rthetas, 89.99);
    get_indice(xozone, ozone, iozone, rozone, 549.99);
    get_indice(xtaucl, taucl, itaucl, rtaucl, 63.99);
    get_indice(xalb, alb, ialb, ralb, 0.9499);

    // std::cout << thetas << " " << ithetas << " " << rthetas << std::endl;
    // std::cout << ozone << " " << iozone << " " << rozone << std::endl;
    // std::cout << taucl << " " << itaucl << " " << rtaucl << std::endl;
    // std::cout << alb << " " << ialb << " " << ralb << std::endl;

    // Start interpolation
    int zthetas, zozone, ztaucl;

    // Remove the dimension on Surface Albedo
    for (int i = 0; i <= 1; i++) {

        zthetas = ithetas + i; // Need to fix

        // std::cout << zthetas << std::endl;

        for (int j = 0; j <= 1; j++) {

            zozone = iozone + j; // Need to fix
            // std::cout << zozone << std::endl;

            for (int k = 0; k <= 1; k++) {

                ztaucl = itaucl + k; // Need to fix
                // std::cout << ztaucl << std::endl;

                for (int l = 0; l < nwl; l++) {
                    // Line 128
                    ed_tmp4[l][i][j][k] =
                        ((1 - ralb) * ed_lut[l][zthetas][zozone][ztaucl][ialb]) +
                        (ralb * ed_lut[l][zthetas][zozone][ztaucl][ialb + 1]);
                }
            }
        }
    }

    // Remove the dimension on taucl
    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= 1; j++) {
            for (int l = 0; l < nwl; l++) {
                ed_tmp3[l][i][j] =
                    (1 - rtaucl) * ed_tmp4[l][i][j][0] + rtaucl * ed_tmp4[l][i][j][1];
                // std::cout << ed_tmp3[l][i][j] << std::endl;
            }
        }
    }

    // Remove the dimension on ozone
    for (int i = 0; i <= 1; i++) {
        for (int l = 0; l < nwl; l++) {
            ed_tmp2[l][i] =
                (1 - rozone) * ed_tmp3[l][i][0] + rozone * ed_tmp3[l][i][1];
            // std::cout << ed_tmp2[l][i] << std::endl;
        }
    }

    // Remove the dimention on sunzenith angle
    for (int l = 0; l < nwl; l++) {
        ed[l] = (1 - rthetas) * ed_tmp2[l][0] + rthetas * ed_tmp2[l][1];
        // std::cout << ed[l] << std::endl;
    }

    return (ed);
}


Rcpp::NumericVector compute_ed0(double o3, double tcl, double cf, double thetas, double albedo, float ed_lut[NBWL][NTHETAS][NO3][NTAUCLD][NALB]) {

    std::vector<float> ed_cloud = interpol_ed0moins(ed_lut, thetas, o3, tcl, albedo);
    std::vector<float> ed_clear = interpol_ed0moins(ed_lut, thetas, o3, 0, albedo);

    Rcpp::NumericVector ed(NBWL);

    const int nwl = ed_cloud.size();

    for (int i = 0; i < nwl; i++) {

        if (thetas < 90) {
            ed[i] = (ed_cloud[i] * cf) + (ed_clear[i] * (1 - cf));
        } else {
            ed[i] = 0;
        }
    }

    return ed;
}

// [[Rcpp::export]]
Rcpp::NumericVector ed0_(int yday, double hour, double lat, double lon, double tcl, double o3, double cf, double albedo) {

    const char* filename = "/home/pmassicotte/Desktop/ed0/inst/extdata/Ed0moins_LUT_5nm_v2.dat";
    float downward_irriadiance_table_as_output[NBWL][NTHETAS][NO3][NTAUCLD][NALB];
    read_ed0moins_lut_(filename, downward_irriadiance_table_as_output);

    double minute = 0;
    double second = 0;

    double theta = sun_zenithal_angle(yday, hour, minute, second, lat, lon);

    Rcpp::NumericVector ed = compute_ed0(o3, tcl, cf, theta, albedo, downward_irriadiance_table_as_output);

    return ed;

}

// [[Rcpp::export]]
std::vector<double> edpar_(Rcpp::IntegerVector yday, Rcpp::NumericVector hour, Rcpp::NumericVector lat, Rcpp::NumericVector lon, Rcpp::NumericVector tcl, Rcpp::NumericVector o3, Rcpp::NumericVector cf, Rcpp::NumericVector albedo) {

  const char* filename = "/home/pmassicotte/Desktop/ed0/inst/extdata/Ed0moins_LUT_5nm_v2.dat";
  float downward_irriadiance_table_as_output[NBWL][NTHETAS][NO3][NTAUCLD][NALB];
  read_ed0moins_lut_(filename, downward_irriadiance_table_as_output);

  double minute = 0;
  double second = 0;

  int n = yday.size();

  // arma::vec par = arma::zeros<arma::vec>(n);

  std::vector<double> par(n);


  for(int i = 0; i < n; i++) {

    double theta = sun_zenithal_angle(yday[i], hour[i], minute, second, lat[i], lon[i]);

    // **************
    // Problème ici, si un vecteur est de longeur 1, quand i = 2 ca sort un chiffre bidon.
    //std::cout << cf[i] << " " << albedo[i] << std::endl;
    // **************

    Rcpp::NumericVector ed = compute_ed0(o3[i], tcl[i], cf[i], theta, albedo[i], downward_irriadiance_table_as_output);

    arma::vec v(ed.begin(), ed.size(), false, true);
    arma::vec wl = arma::regspace(290, 5, 700);
    arma::uvec ids = find(wl >= 400);

    arma::vec wl2 = wl.elem(ids);
    arma::vec v2 = v.elem(ids);

    // std::cout << i << std::endl;

    // arma::vec vv = arma::trapz(wl2, v2);
    par[i] = (as_scalar(arma::trapz(wl2, v2)));

    // Rcpp::Rcout << par[i] << std::endl ;

    // std::cout << vv << std::endl;
    // std::cout << as_scalar(arma::trapz(wl2, v2));
  }

  return par;

}

/*** R
library(tidyverse)

dd <- ed0_(200, 12, 67.47973, -63.78953, 3, 330, 1, 0.05)

edpar_(100, 12, 67.47973, -63.78953, 3, 330, 1, 0.05)
#
n <- 2
res <- edpar_(rep(200, n), rep(12, n), rep(67.47973, n), rep(-63.78953, n), rep(3, n), rep(330, n), rep(1, n), rep(0.05, n))

edpar_(rep(200, n), rep(12, n), rep(67.47973, n), rep(-63.78953, n), rep(3, n), rep(330, n), 1, 0.05)

#
#
# df <- tibble(wl = seq(290, 700, by = 5), dd) %>%
  # filter(between(wl, 400, 700))

# pracma::trapz(df$wl, df$dd)


# for(i in 1:10) {
#   ed0(200, 12, 67.47973, -63.78953, 3, 330, 1, 0.05)
# }


*/
