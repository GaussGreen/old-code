/* ====================================================================================
   FILENAME:  srt_h_etabeta.h

   PURPOSE:   provide a few useful functions when using the Beta-Power model
   ====================================================================================
 */

#ifndef SRT_H_ETABETAMODEL_H
#define SRT_H_ETABETAMODEL_H

/* -------------------------------------------------------------------------------
   Computes M(t  ,x  ,T) defined by log E[exp-h(T)[X_T-x_t] | x_t]
   -------------------------------------------------------------------------------
 */

Err M_t_x_T_function(double h_T, double power);

/* ------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------
   The h(t) function is defined as :
                  h(x  ,t) = h0(t) * h1(x)
   where:
           h0(t) =                             --> h_time_coeff_fct
   --------------------------------------------------------------------- */

double h_t_etabeta_function(double time, TermStruct *ts);

#endif
