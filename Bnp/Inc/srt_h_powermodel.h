/* ====================================================================================
   FILENAME:  srt_h_powermodel.h

   PURPOSE:   provide a few useful function swhen using Pat's model
   ====================================================================================
 */

#ifndef SRT_H_POWERMODEL_H
#define SRT_H_POWERMODEL_H

/* -------------------------------------------------------------------------------
   Computes log_S(t  ,x  ,T) defined by log E[exp-h(T)*X_T|x_t]
   -------------------------------------------------------------------------------
 */

Err logS_function(double h_T, double cond_mean_X_T, double cond_var_X_T,
                  double power, double *log_S);

double x_power_fct(double x, double power);

double h_t_power_function(double t_in_y, TermStruct *ts);

#endif
