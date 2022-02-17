
#include "pde_h_types.h"

Err srt_f_SABW_Opt(Date EvalDate, Date ExpDate, double OptStrike,

                   char *rec_pay_str, char *greek_str,

                   char *yield_crv_name, char *grow_crv_name,

                   /* LOCAL VOLATILITY PARAMETERS */
                   double S, double A, double B, double W,

                   /* PDE PARAMETERS */
                   long MinNode, long MinNumMesh, double X0, double Xmin,
                   double Xmax, double dXmax,

                   double *PV);