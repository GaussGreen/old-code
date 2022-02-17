
#include "srt_h_lgmtypes.h"

double H2_Func(Date tRef, Date *tStart, Date tNow, long nEx, double *H1,
               double gamma);
double H1_Func(Date tRef, Date *tStart, Date tNow, long nEx, double *H1);
double XLH1_Func(Date tRef, Date *tStart, Date tNow, Date tLast, long nEx,
                 double *H1);
double XLH2_Func(Date tRef, Date *tStart, Date tNow, Date tLast, long nEx,
                 double *H1, double gamma);

Err lgm_2f_autocal_calibrate(
    void *dealPtr, LGM_TS **lgm_ts_ptr_ptr,
    SrtLgmRefSwptnData *lgm_ref_swp_data, LGMCalParm *cal_req,
    Err (*srt_f_get_vol)(Date, Date, double, SRT_Boolean, double *),
    Err (*srt_f_get_beta)(Date, Date, double *), String yc_name);

Err XLLGMZeta1Func(double t, char *UndName, double *Zeta);
Err XLLGMZeta2Func(double t, char *UndName, double *Zeta);
Err XLLGMZeta12Func(double t, char *UndName, double *Zeta);

Err XLLGMH1Func(double t, char *UndName, double *H);
Err XLLGMH2Func(double t, char *UndName, double *H);

double XL_New_Zeta2_Func(Date tRef, Date *tEx, long nEx, Date tNow,
                         double *Zeta1, double alpha, double gamma, double rho);

double XL_New_Zeta12_Func(Date tRef, Date *tEx, long nEx, Date tNow,
                          double *Zeta1, double alpha, double gamma,
                          double rho);

double New_Zeta2_Func(Date tRef, Date *tEx, long nEx, Date tNow, double *Zeta1,
                      double *Zeta2, double alpha, double gamma, double rho);
double New_Zeta12_Func(Date tRef, Date *tEx, long nEx, Date tNow, double *Zeta1,
                       double *Zeta12, double alpha, double gamma, double rho);

double Zeta1_Func(Date tRef, Date *tEx, long nEx, Date tNow, double *Zeta1);
