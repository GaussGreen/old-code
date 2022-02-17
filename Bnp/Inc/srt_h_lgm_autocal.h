#ifndef SRT_H_LGM_AUTOCAL_H
#define SRT_H_LGM_AUTOCAL_H

#define DEFAULT_TAU 20

typedef enum srtstrikemethod {
  IRRSTRIKE,
  STDEVSTRIKE,
  MIXEDSTRIKE
} SrtStrikeMethod;

static Err interp_strike_method(String strike_method_str,
                                SrtStrikeMethod *strikemethod);

Err srt_f_lgm_autocal(Err (*GetVol)(), String vol_type_str, SrtUndPtr und,
                      Date *date, long ndates, long lastexdateindex,
                      double *main_strike, double *aux_strike,
                      SrtCompounding compounding, SrtBasisCode main_basis,
                      SrtBasisCode aux_basis,
                      /* OUTPUT */
                      Date **g_date, double **g, long *num_g, Date **zeta_date,
                      double **zeta, long *num_zeta, long update_ts_flag,
                      TermStruct **ts);

static Err eval_lgm_main_swaption_sqzeta_g(double sqzeta, double g,
                                           double *price, double *dp_dz,
                                           double *dp_dg);

static Err eval_lgm_aux_swaption_sqzeta_g(double sqzeta, double g,
                                          double *price, double *dp_dz,
                                          double *dp_dg);

static Err eval_lgm_swaption_sqzeta(double sqzeta, double *price,
                                    double *dp_dz);

static Err eval_lgm_swaption_g(double g, double *price, double *dp_dg);

static Err eval_lgm_swaption_fct(double sqzeta, double g, double *price,
                                 double *dp_dz, double *dp_dg);

static Err eval_lgm_swap_fct(double phi, double *price, double *deriv);

static Err new_g_func(double x, double *g_x, double *dg_dx);

static void quickly_fill_in_g_array(double *g, double *mat, long index_start,
                                    long index_end);

Err srt_f_update_ts_from_zeta_g(TermStruct **ts, Date tnow, Date *zeta_date,
                                double *zeta, long num_zeta, Date *g_date,
                                double *g, long num_g);

#endif
