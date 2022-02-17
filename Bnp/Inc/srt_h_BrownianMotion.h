
#include "pde_h_types.h"

typedef struct {
  double equity_vol;
  double forward;
  double yr_to_exp;
  SrtReceiverType pay_rec;
  double equity_strike;

} srt_black_scholes_parms;

Err srt_f_pde_eur_option(Date eval_date, Date exp_date, double forward,
                         double equity_strike, double equity_vol,
                         char *rec_pay_str, char *pde_solving_scheme_str,
                         long min_node, long num_mesh,

                         /*OUTPUTS */
                         double *option_pv);

Err srt_f_pde_amer_option(Date eval_date, Date *tEx, long nEx, double forward,
                          double equity_strike, double equity_vol,
                          char *rec_pay_str, char *pde_solving_scheme_str,
                          char *PDEBoundaryCondStr, long min_node,
                          long num_mesh,

                          /*OUTPUTS */
                          double *option_pv);