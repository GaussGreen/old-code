/* ======================================================
   FILENAME:  utgreeks.h

   PURPOSE:   Greeks
   ====================================================== */

#ifndef UTGREEKS_H
#define UTGREEKS_H

typedef int Message;

struct greek_struct
{
    double fwd_delta_underlying;
    double spot_delta_underlying;
    double fwd_delta_strike;
    double spot_delta_strike;
    double fwd_gamma;
    double spot_gamma;
    double fwd_gamma_strike;
    double spot_gamma_strike;
    double spot_gamma_cross;
    double vega;
    double theta;
    double rho;
    double zeta;
    double theta_stochastic_vol;
};
typedef struct greek_struct Greekstruct;

typedef enum SrtGreekType
{
    FIRSTGREEKTYPE = -1,
    PREMIUM,
    DELTA,
    DELTA_FWD,
    DELTA_FWD1,
    DELTA_FWD2,
    DELTAX,
    DELTAY,
    DELTAZ,
    GAMMA,
    GAMMAX,
    GAMMAY,
    GAMMAZ,
    GAMMAXY,
    GAMMAYZ,
    GAMMAXZ,
    GAMMA_FWD,
    GAMMA_FWD1,
    GAMMA_FWD2,
    THETA,
    THETA_1D,
    THETAGAMMA,
    THETAVOLGA,
    THETAVANNA,
    THETACARRY,
    VEGA,
    VEGAX,
    VEGAY,
    VEGAZ,
    VEGAXY,
    VEGA1,
    VEGA2,
    VANNA,
    VOLGA,
    RHO,
    ZETA,
    RHORHO,
    RHODELTA,
    LASTGREEKTYPE,
    ALPHA,
    BETA,
    DU1,      /* sensitivity to the positive amplitude, in the Merton model */
    DLAMBDA1, /* sensitivity to the positive intensity*/
    DU2,      /* sensitivity to the negative amplitude */
    DLAMBDA2, /*  sensitivity to the negative intensity*/
    DENSITY,
    CUMDENSITY
} SrtGreekType;

#endif
