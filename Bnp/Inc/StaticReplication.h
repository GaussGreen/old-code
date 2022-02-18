#ifndef __STATICREPLICATION_H
#define __STATICREPLICATION_H

#include "srt_h_all.h"

Err Static_payoff_replication(
    double* tangent_points,
    int     n_tangent_points,
    double  cvx_strike,
    int     is_left_flat,
    int     is_right_flat,
    Err (*payoff)(double S, double K, double* out),
    Err (*deriv)(double S, double K, double* out),
    double* constant,
    double* strikes,
    double* notionals,
    int*    is_call);

Err convex_payoff(double S, double K, double* out);
Err convex_deriv(double S, double K, double* out);

Err Convex_replication(  // Market
    long    Today,
    long    fx_spot_date,
    double  Spot,
    double* fx_mkt_vol,
    long*   fx_mkt_vol_date,
    double* fx_mkt_smile_alpha,
    double* fx_mkt_smile_beta,
    double* fx_mkt_smile_rho,
    double* fx_mkt_smile_pi,
    int     num_fx_mkt_vol,
    int     smile_spec_type,  //	0: lognormal vol + SABR params, 1: sigma-beta + SABR params, 2:
                              //lognormal vol + BMM, 3: lognormal vol + BMM
    char* for_yc,
    char* dom_yc,
    int   is_reversed_fx,
    // Convex Data
    double Notional,
    long   Fix_date,
    long   Settle_date,
    double Kmin,
    double Kmax,
    double Convex_strike,
    // algo data
    int nb_strikes,
    // Output
    double** Replication,
    double*  Price);

Err Convex_replication_simple(  // Market
    long   Today,
    double Forward,
    double fx_mkt_vol,
    double fx_mkt_smile_alpha,
    double fx_mkt_smile_beta,
    double fx_mkt_smile_rho,
    double fx_mkt_smile_pi,
    int    smile_spec_type,  //	0: lognormal vol + SABR params, 1: sigma-beta + SABR params, 2:
                             //lognormal vol + BMM, 3: lognormal vol + BMM
    int is_reversed_fx,
    // Convex Data
    double Notional,
    long   Fix_date,
    long   Settle_date,
    double Kmin,
    double Kmax,
    double Convex_strike,
    // algo data
    int nb_strikes,
    // Output
    double** Replication,
    double*  Price);
#endif
