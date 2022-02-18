#ifndef __BREAK_EVEN_H
#define __BREAK_EVEN_H

Err BreakEvenFromMkt(
    char*    yc_name,
    char*    vc_name,
    long     spot_date,
    int      spot_lag,
    char*    swap_freq,
    char*    swap_basis,
    char*    ref_rate_name,
    int      NMaturity,
    char**   maturity_tenor,
    int      NUnderlying,
    char**   underlying_tenor,
    int      NbOfBDayPerYear,
    double** BEMatrix);

Err BreakEvenFromVol(
    long     spot_date,
    int      spot_lag,
    int      NMaturity,
    char**   maturity_tenor,
    int      NUnderlying,
    char**   underlying_tenor,
    int      NbOfBDayPerYear,
    double** VolMatrix,
    double** BEMatrix);

#endif