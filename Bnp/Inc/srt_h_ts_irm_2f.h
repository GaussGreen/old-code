/* ============================================================================

   FILENAME:   SRT_H_TS_IRM_2F_INIT.H

   PURPOSE:    Initialise a volatility term structure for a two factor model  ,
               straight from raw data (dates + values)

   ============================================================================
 */

#ifndef SRT_H_TS_IRM_2F_H
#define SRT_H_TS_IRM_2F_H

typedef double SrtTFTSMat[2][2];
typedef double SrtTFTSVec[2];

/* ----------------------------------------------------------------------------------
                                                                PART I
                                MAIN FUNCTION THAT DEFINES A TERM STRUCT
   ----------------------------------------------------------------------------------
 */
Err srt_f_init_IRM_TwoFac_TermStruct(
    TermStruct **list, /* Pointer to the Term Struct (Filled on OUTPUT) */

    Date today, /* Today's date */

    double **sig_data, /* Sigma Curve: [0]Date [1]sig1 [2]sig2 [3]rho */
    int sig_cols,      /* Number of Columns (4) */
    int num_sigs,      /* Number of dates for discretisation */
    double **tau_data, /* Tau Curve: [0]Date [1]tau1 [2]tau2 */
    int tau_cols,      /* Number of Coulmns (3) */
    int num_taus,      /* Number of dates for discretisation */

    SrtMdlType mdl_type,

    /* SMILE */
    double beta, /* For Cheyette Beta */

    /* CORRELATION */
    double alpha, double gamma, double rho,

    /* MIXED SMILES */
    double omega);
/* ----------------------------------------------------------------------------
 */

/* Function required to free the TermStruct once attached into a linked list */

Err srt_f_irm2ftsvalfree(void *tsvalptr);

/* ----------------------------------------------------------------------------
 */

/*  In srt_f_ts_irm_2f_fct. c */

Err find_2f_tau(double time, TermStruct *l, double *tau1, double *tau2);

Err find_2f_sig(double time, TermStruct *l, double *sig1, double *sig2);

Err find_2f_rho(double time, TermStruct *l, double *rho);

Err find_tf_beta(double time, TermStruct *l, double *beta1, double *beta2);

double find_tf_omega(double time, TermStruct *l);

Err srt_f_display_IRM_TwoFac_TermStruct(TermStruct *ts, double **sigma_date,
                                        double **sigma1, double **beta1,
                                        double **sigma2, double **beta2,
                                        double **rho, long *plNumSigmas,
                                        double **tau_date, double **tau1,
                                        double **tau2, long *plNumTaus);

Err get_2f_F_funcs(double time, TermStruct *l, SrtTFTSVec *F);

Err get_2f_Psi_funcs(double time, TermStruct *l, SrtTFTSVec *Psi);

Err get_2f_sig2_rho_interp(double time1, double time2, TermStruct *l,
                           double *sig1_s, double *sig2_s, double *rho);

Err get_2f_G_funcs(double time, TermStruct *l, SrtTFTSMat *val_G);

Err get_2f_H_funcs(double time, TermStruct *l, SrtTFTSMat *val_H);

Err get_2f_J_funcs(double time, TermStruct *l, SrtTFTSVec *J);

Err make_2f_lgm_cum_vols(void *tts, double fixing_time, double start_time,
                         double end_time, double *Cum_Vol3, double *Cum_Vol4,
                         double *Var);

Err find_tf_beta(double time, TermStruct *l, double *beta1, double *beta2);

Err make_2f_F_Psi_vectors(TermStruct *l);

Err make_2f_G_H_vectors(TermStruct *l);

#endif
