#include "srt_h_all.h"

Err srt_f_get_likely_most_expensive_european(
    Date*          midat_start_dates,
    Date           midat_theo_end_date,
    String         ref_rate_code,
    SrtCompounding fixed_freq,
    BasisCode      fixed_basis,
    long           nEx,
    Err (*get_vol)(Date, Date, double, SRT_Boolean, double*),
    String           yc_name,
    SrtDiffusionType swp_vol_type,
    char*            pay_rec_str,
    double*          swp_strikes,

    double** swp_fwd_vol,

    long** likely_most_expensive);