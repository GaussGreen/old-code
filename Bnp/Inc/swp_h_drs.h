#ifndef SWP_H_DRS_H
#define SWP_H_DRS_H

Err srt_f_DRS(long start, long end_or_nfp, String basisStr, double strike,
              String recPayStr, /*call or put DRS*/
              Err (*GetVol)(char *vc_id, double exercise, double end,
                            double strike, double *bs_vol),
              char *vc_id, String vol_type, SRT_Boolean spread_adjust_vol,
              String refRateCode, String yc_id,
              int num_caplets,     /*nb of caplets used*/
              double delta_strike, /*gap between two strikes*/
              double *ans);
#endif
