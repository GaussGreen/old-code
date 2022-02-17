/*--------------------------------------------------------------
        FILE: LGM1FCalibToProduct.h
        PURPOSE: LGM1f calibration to any product
        AUTHOR: Dimitri Mayevski
        DATE: 11/07/2002
  --------------------------------------------------------------*/

#ifndef __LGM1FCALIBTOPRODUCT_H__
#define __LGM1FCALIBTOPRODUCT_H__

Err srt_f_lgm1f_calib_to_product(char **yc_names, int n_mkts, int mkt_idx,
                                 SrtProduct *product, int nex, long *ex_dates,
                                 double lam, double *sig, int nx,
                                 double qto_fudge);

#endif /* #ifndef __LGM1FCALIBTOPRODUCT_H__ */