/* ----------------------------------------------------------------------------

   FILENAME:          srt_h_ts_FX_int.h

   PURPOSE:           functions to initialise an FX underlying TermStruct

   ----------------------------------------------------------------------------
 */

#ifndef SRT_H_TS_EQ_H
#define SRT_H_TS_EQ_H

#include "srt_h_ts.h"

/* -------------------------------------------------------------------------------
   The MAIN function to call when a TS has to be initialised for an EQD
   underlying
   -------------------------------------------------------------------------------
 */

Err srt_f_init_EQ_TermStruct(Date today, double **sig, /* sig[col][row] */
                             int sig_cols, int nsig,

                             char *undname, char *iror_yc_name,
                             SrtMdlType mdl_type,

                             /* SRVGS Model Parameter */
                             double omega, double beta, double gamma,

                             double voldrift, double vovol, double rho,

                             /* Output */
                             TermStruct **ts);

/* -------------------------------------------------------------------------------
   The MAIN function to free a Term Structure of Volatility initialised
   for an EQD underlying
   -------------------------------------------------------------------------------
 */

Err srt_f_free_EQ_TermStruct(TermStruct **l);

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------------
                              NON SMILE SPECIFIC FUNCTIONS
   -------------------------------------------------------------------------------
 */

/* Function required to free a TermStructObject once attached into a linked list
 */
Err srt_f_EquityTermStructValFree(void *tsvalptr);

/* -------------------------------------------------------------------------- */
/* Initialises the Volatility Term Structure for Non smile Eq Models */
Err srt_f_init_EQ_NORMLOG_TermStruct(char *undname, char *ir_or_yc_name,
                                     SrtMdlType mdl_type, TermStruct **ts,
                                     Date today, double **sig_data,
                                     int sig_cols, int num_sig,

                                     double omega, double beta, double gamma,

                                     double voldrift, double vovol, double rho);

/* -------------------------------------------------------------------------- */
/* A function to output an existing TermStructure (Memory Allocation done
 * inside) */

Err srt_f_display_EQ_TermStruct(TermStruct *ts, double **sigma_date,
                                double **sigma, long *sigma_n);

/* -------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------
                 IN SRT_F_TS_EQ_FCT.c
   -------------------------------------------------------------------------- */

/* Finds the Value of the loval volatility at one given point in time */
double find_eq_sig(double time, TermStruct *l);

/* Returns the Cumulative Volatility of the Equity Forward */
double eq_cum_vol_func(double time, SRT_Boolean sigsq, TermStruct *l);

Err srt_f_eq_implied_vol(double expiry, char *und_name, double *eq_bs_vol);

Err srt_f_eq_calib(long n_ex, double *dates, double *vols,

                   char *ir_und_name, char **und_name);

Err srt_f_eq_forward(Date date, SrtUndPtr und, char *yc_name,

                     double *forward);

double find_eq_omega(double time, TermStruct *l);
double find_eq_beta(double time, TermStruct *l);
double find_eq_gamma(double time, TermStruct *l);
double find_eq_basevol(double time, TermStruct *l);
double find_eq_voldrift(double time, TermStruct *l);
double find_eq_vovol(double time, TermStruct *l);
double find_eq_rho(double time, TermStruct *l);

#endif
