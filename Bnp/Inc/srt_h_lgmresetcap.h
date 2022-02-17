#ifndef _LGMRESETCAP_H
#define _LGMRESETCAP_H

#define ERROR_IN_RESETCAP 9999.9999

Err srt_f_lgm_phifunc(SrtUndPtr und, Ddate fixing, double phi[2][2]);

Err srt_f_lgm_resetcaplet(SrtUndPtr und, Ddate str_fix, Ddate str_start_end[2],
                          double str_cvg, Ddate spot_fix,
                          Ddate spot_start_end[2], double spot_cvg,
                          Ddate pay_date, double pay_level,
                          SrtReceiverType rec_pay, double *premium,
                          double fixed_strike);

double gen_ir_resetcap(SrtUndPtr und, long *fixing_date,
                       long *start_end_pay_dts, double *df, double *cvg,
                       long nfp, SrtReceiverType rec_pay, double first_fixing);

#endif
