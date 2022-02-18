/* ----------------------------------------------------------------------------

   FILENAME:          srt_h_ts_FX.h

   PURPOSE:           functions to initialise and work with an FX underlying
                      TermStruct

   ---------------------------------------------------------------------------- */

#ifndef SRT_H_TS_FX_H
#define SRT_H_TS_FX_H

#include "srt_h_ts.h"

/* To be find with "find in file":  PVMPI*/
/* The following lines must be included in changes */
#ifdef __cplusplus
extern "C"
{
#endif

    /* -------------------------------------------------------------------------------
       The MAIN function to call when a TS has to be initialised for an FX underlying
       ------------------------------------------------------------------------------- */

    Err srt_f_init_FX_TermStruct(
        Date     today,
        double** sig, /* sig[col][row] */
        int      sig_cols,
        int      nsig,

        SrtMdlType mdl_type,

        /* SMILE */
        double beta,

        /* STOCH VOL */
        double vovol,

        /* For FX_STOCH_RATE */
        char* undName,
        char* domirundName,
        char* forirundName,

        /* OUTPUT */
        TermStruct** ts);

    Err srt_f_display_FX_TermStruct(
        char* szFxUndName,

        long*    sigma_n,
        double** sigma_date,
        double** sigma,

        long*    corr_date_n,
        double** corr_date,
        double** corr);

    Err srt_f_display_FXBS_TermStruct(
        char* szFxUndName,

        long*    sigma_n,
        double** sigma_date,
        double** sigma);

    /* -------------------------------------------------------------------------- */

    /* -------------------------------------------------------------------------------
       The MAIN function to free a Term Structure of Volatility initialised
       for a FOREIGN EXCHANGE underlying
       ------------------------------------------------------------------------------- */

    Err srt_f_free_FX_TermStruct(TermStruct** l);

    /* -------------------------------------------------------------------------------
                                  NON SMILE SPECIFIC FUNCTIONS
       ------------------------------------------------------------------------------- */

    /* Function required to free a TermStructObject once attached into a linked list */
    Err srt_f_fxtsvalfree(void* tsvalptr);

    /* Initialises the Volatility Term Structure for Non smile FX Models */

    Err srt_f_init_FX_NORMLOG_TermStruct(
        TermStruct** ts, Date today, double** sig_data, int sig_cols, int num_sig);

    Err srt_f_init_FX_STOCH_RATES_TermStruct(
        char*        undName,
        SrtUndPtr    dom_und,
        SrtUndPtr    for_und,
        TermStruct** ts,
        Date         today,
        double**     sig_data,
        int          sig_cols,
        int          num_sig);

    /* -------------------------------------------------------------------------- */

    /* --------------------------------------------------------------------------
                     IN SRT_F_TS_FX_FCT.c
       -------------------------------------------------------------------------- */

    /* ------------------------------------------------------------------------ */

    /* Finds the Value of the loval volatility at one given point in time */

    double find_fx_sig(double time, TermStruct* l);

    /* ---------------------------------------------------------------------------- */

    double fx_cum_vol_func(double time, SRT_Boolean sigsq, TermStruct* l);

    /* ---------------------------------------------------------------------------------------
                                                            For FX Stoch Rates with Jumping
       Numeraire
       -----------------------------------------------------------------------------------------*/

    /* Compute the different functions needed and hang them to each point of the term structure  */

    double I_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double K_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double L_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double O_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double Q_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double R_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double S_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double V_dx_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double W_dx_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double O_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double Phi_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double Phi_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double H_fd_func(double time, TermStruct* domts, TermStruct* forts, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double M_fx_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double N_fx_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double P_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double R_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double M_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double Q_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double S_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double T_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double U_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double V_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double W_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double X_fd_func(double time, TermStruct* l);
    /* --------------------------------------------------------------------------------------- */

    double Y_fd_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double V_fx_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    double W_fx_func(double time, TermStruct* l);

    /* --------------------------------------------------------------------------------------- */

    Err srt_f_get_fx_implied_vol(double yr_to_exp, String fx_und_name, double* fx_bs_implied_vol);

    Err srt_f_get_fx_stoch_rates_ind_correl(
        double  start_date,
        double  end_date,
        String  first_und_name,
        String  second_und_name,
        String  fx_und_name,
        double* corr);
/* To be find with "find in file":  PVMPI*/
/* The following lines must be included in changes */
#ifdef __cplusplus
}
#endif

#endif
