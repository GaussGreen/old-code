#ifndef FxSabrGrfnH
#define FxSabrGrfnH

#include "grf_h_mdlcomm.h"

typedef struct {
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  int fx_idx;
  int dom_idx; /* -1 if none */
  int for_idx; /* -1 if none */

  long num_fx_df;
  double *fx_df_tms;
  long *fx_df_dts;

  long num_dom_df;
  double *dom_df_tms;
  long *dom_df_dts;

  long num_for_df;
  double *for_df_tms;
  long *for_df_dts;

} grfn_parm_fx_sabr, *GRFNPARMFXSABR;

Err payoff_fx_sabr_adi(double evt_date, double evt_time, void *func_parm,

                       /* Market data	*/
                       long today, double spot_fx, /*	The cash one */
                       void *dom_yc, void *for_yc,

                       /* Grid data	*/
                       int l1, int u1, int l2, int u2, double *x,

                       /* Vector of results to be updated */
                       int nprod, double ***prod_val);

Err payoff_fx_sabr_adi_ATM_opt(double evt_date, double evt_time,
                               void *func_parm,

                               /* Market data	*/
                               long today, double spot_fx, /*	The cash one */
                               void *dom_yc, void *for_yc,

                               /* Grid data	*/
                               int l1, int u1, int l2, int u2, double *x,

                               /* Vector of results to be updated */
                               int nprod, double ***prod_val);

Err payoff_fx_sabr_adi_opt(double evt_date, double evt_time, void *func_parm,

                           /* Market data	*/
                           long today, double spot_fx, /*	The cash one */
                           void *dom_yc, void *for_yc,

                           /* Grid data	*/
                           int l1, int u1, int l2, int u2, double *x,

                           /* Vector of results to be updated */
                           int nprod, double ***prod_val);

#endif