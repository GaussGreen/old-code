#ifndef SRT_MIDAT_AUTOCAL_H
#define SRT_MIDAT_AUTOCAL_H

String srt_f_lgm_midat_autocal(
    long In_ex, Date It_ex[], Date It_start[], double Iprem[], long npay,
    Date tfix_start[], Date tfix_end[], Date tpay[], double fixpymnt[],
    double Iredpay[], SrtReceiverType p_r_flag, SrtUndPtr und,
    long end_of_day_flag, String (*GetVol)(), String vol_type_str,
    int update_ts_flag, int find_exer_bdry, String outfile, double *LGMvalue,
    long *n_exer_bdryptr, double *r_exer_bdryptr[], Date *t_exer_bdryptr[]);

String LGMamerican(long nfix, Date tfix_start[], Date tfix_end[],
                   Date tfix_pay[], double full_pay_fix[], double Inprem[],
                   BasisCode fix_basis, int early_flag_fix, long nflt,
                   Date t_flt[], BasisCode flt_basis, int early_flag_flt,
                   Date t_first_fix, int lag_exer_start, int cal_or_bus,
                   BusDayConv conv_start, SrtReceiverType rec_pay,
                   SrtUndPtr und, long end_of_day_flag, String (*GetVol)(),
                   String vol_type_str, int update_ts_flag, int find_exer_bdry,
                   String outfile, double *LGMamervalue, long *n_exer_bdryptr,
                   double *r_exer_bdryptr[], Date *t_exer_bdryptr[]);

#endif
