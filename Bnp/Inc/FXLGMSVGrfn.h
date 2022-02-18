#ifndef FXLGMSVGRFN
#define FXLGMSVGRFN

#include "LGMSVUtil.h"
#include "grf_h_mdlcomm.h"
#include "math.h"
#include "srt_h_all.h"

typedef struct
{
    GRFNCOMMSTRUCT global;
    FIRSTMktAtT*   local;

    int fx_idx;
    int dom_idx; /* -1 if none */
    int for_idx; /* -1 if none */

    long    num_fx_df;
    double* fx_df_tms;
    long*   fx_df_dts;
    double* fx_dff;
    double* fx_gam1;
    double* fx_gam2;
    double* fx_gam1_2;
    double* fx_gam2_2;
    double* fx_gam12;
    int     do_fx;

    long    num_dom_df;
    double* dom_df_tms;
    long*   dom_df_dts;
    double* dom_dff;
    double* dom_gam1;
    double* dom_gam2;
    double* dom_gam1_2;
    double* dom_gam2_2;
    double* dom_gam12;
    int     do_dom;

    long    num_for_df;
    double* for_df_tms;
    long*   for_df_dts;
    double* for_dff;
    double* for_gam1;
    double* for_gam2;
    double* for_gam1_2;
    double* for_gam2_2;
    double* for_gam12;
    int     do_for;

} grfn_parm_mc_fxlgmsv, *GRFNPARMMCFXLGMSV;

Err payoff_fxlgmsv_mc(
    long   path_index,
    double evt_date,
    double evt_time,
    void*  func_parm,

    // Domestic
    double domft1,
    double domft2,
    double dompsi1,
    double dompsi2,
    double dompsi12,
    double domv,

    // Foreign
    double forft1,
    double forft2,
    double forpsi1,
    double forpsi2,
    double forpsi12,
    double forv,

    // FX
    double fx_spot,

    /* Vector of results to be updated */
    int nprod,
    /* Result	*/
    double* prod_val,
    int*    stop_path);

Err payoff_qtolgmsv1f_mc(
    long   path_index,
    double evt_date,
    double evt_time,
    void*  func_parm,

    // Domestic
    double domft,
    double dompsi,

    // Foreign
    double forft,
    double forpsi,
    double forv,

    /* Vector of results to be updated */
    int nprod,
    /* Result	*/
    double* prod_val,
    int*    stop_path);

Err payoff_qtolgmsv2f_mc(
    long   path_index,
    double evt_date,
    double evt_time,
    void*  func_parm,

    // Domestic
    double domft,
    double dompsi,

    // Foreign
    double forft1,
    double forft2,
    double forpsi1,
    double forpsi2,
    double forpsi12,
    double forv,

    /* Vector of results to be updated */
    int nprod,
    /* Result	*/
    double* prod_val,
    int*    stop_path);

#endif
