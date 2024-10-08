/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        Description	:



        Author		:

        Created		:

        History		:

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
// prevent multiple inclusions
#ifndef LGM2F_CMS_closedform
#define LGM2F_CMS_closedform

//#include "srt_h_lgmtypes.h"

#include "CPDCalib.h"
#include "srt_h_und_struct.h"

#pragma once

double compute_sum_sigma(double t, double *sigma, double *sigma_date,
                         int nb_sigma, double lambda);
double compute_sum_diagterm(double t, double *sigma_1, double *sigma_2,
                            double *sigma_date, int nb_sigma, double lambda_1,
                            double lambda_2);
double compute_A(double t, double T, double sum_sigma_1, double sum_sigma_2,
                 double sum_diagterm, double rho, double lambda_1,
                 double lambda_2, double Tpay);
double compute_c1(double t, double _T, double Tpay, double *sigma_1,
                  double *sigma_2, double *sigma_date, int nb_sigma, double rho,
                  double lambda_1, double lambda_2);
double compute_c2(double t, double _T, double Tpay, double *sigma_2,
                  double *sigma_date, int nb_sigma, double rho,
                  double lambda_2);
double compute_CMS_GaussLeg(double Tf, double Tpay, double *DF_fix,
                            double *fix_coupon_date, double *fix_cvg,
                            int nb_fix_coupon, double *DF_float, double DF_TF,
                            double *float_coupon_date, double *float_cvg,
                            int nb_float_coupon, double *spread,
                            double *sigma_1, double *sigma_2,
                            double *sigma_date, int nb_sigma, double lambda_1,
                            double lambda_2, double rho, int nb_integration,
                            double bound, SrtReceiverType ProductType,
                            double strike);
double compute_CMScov_GL(
    double Tf1, double Tpay1, double *DF_fix1, double *fix_coupon_date1,
    double *fix_cvg1, int nb_fix_coupon1, double *DF_float1, double DF_TF1,
    double *float_coupon_date1, double *float_cvg1, int nb_float_coupon1,
    double *spread1, double Tf2, double Tpay2, double *DF_fix2,
    double *fix_coupon_date2, double *fix_cvg2, int nb_fix_coupon2,
    double *DF_float2, double DF_TF2, double *float_coupon_date2,
    double *float_cvg2, int nb_float_coupon2, double *spread2, double *sigma_11,
    double *sigma_21, double *sigma_date1, int nb_sigma1, double lambda_11,
    double lambda_21, double rho1, double *sigma_12, double *sigma_22,
    double *sigma_date2, int nb_sigma2, double lambda_12, double lambda_22,
    double rho2, int nb_integration, double bound);
Err compute_CMScorrel_GL(
    double Tf1, double Tpay1, double *DF_fix1, double *fix_coupon_date1,
    double *fix_cvg1, int nb_fix_coupon1, double *DF_float1, double DF_TF1,
    double *float_coupon_date1, double *float_cvg1, int nb_float_coupon1,
    double *spread1, double Tf2, double Tpay2, double *DF_fix2,
    double *fix_coupon_date2, double *fix_cvg2, int nb_fix_coupon2,
    double *DF_float2, double DF_TF2, double *float_coupon_date2,
    double *float_cvg2, int nb_float_coupon2, double *spread2, double *sigma_11,
    double *sigma_21, double *sigma_date1, int nb_sigma1, double lambda_11,
    double lambda_21, double rho1, double *sigma_12, double *sigma_22,
    double *sigma_date2, int nb_sigma2, double lambda_12, double lambda_22,
    double rho2, int nb_integration, double bound, double **price);

Err LGM2FCMS(char *yc_name, char *vol_curve_name, SrtUndPtr und, double alpha,
             double gamma, double rho, int do_calib, long lTfixing,
             long lStartDate, long lEndDate, char *fixFreq, char *fixBasis,
             char *refRateName, long lTpay, long lToday, int nb_integration,
             double bound, char *ProductTypeStr, double strike,
             char *(GetCpdAutocalCashVol)(char *lpszVolCurveName, double dStart,
                                          double dEnd, double dStrike,
                                          int bIsACap, char *lpszRefRateCode,
                                          double *pdVol, double *pdPower),
             double *price);
Err LGM2FDRS(char *yc_name, char *vol_curve_name, SrtUndPtr und, double alpha,
             double gamma, double rho, int do_calib, long lTfixing,
             long lStartDate, long lEndDate, char *fixFreq, char *fixBasis,
             char *refRateName, long lTpay, long lToday, int nb_integration,
             double bound,
             char *(GetCpdAutocalCashVol)(char *lpszVolCurveName, double dStart,
                                          double dEnd, double dStrike,
                                          int bIsACap, char *lpszRefRateCode,
                                          double *pdVol, double *pdPower),
             double **output);
Err LGM2FCMSCorrel(char *yc_name, char *vol_curve_name, SrtUndPtr und,
                   double alpha, double gamma, double rho, int do_calib,
                   long lTfixing, long lStartDate, long lEndDate1,
                   char *fixFreq1, char *fixBasis1, char *refRateName,
                   long lTpay, long lToday, long lEndDate2, char *fixFreq2,
                   char *fixBasis2, int nb_integration, double bound,
                   char *(GetCpdAutocalCashVol)(char *lpszVolCurveName,
                                                double dStart, double dEnd,
                                                double dStrike, int bIsACap,
                                                char *lpszRefRateCode,
                                                double *pdVol, double *pdPower),
                   double **price);
Err LGM2FCMSCorrelBuild(char *yc_name, char *und, long lTfixing,
                        long lStartDate, long lEndDate1, long lEndDate2,
                        long lTpay, long lToday, char *fixFreq1,
                        char *fixBasis1, char *refRateName, char *fixFreq2,
                        char *fixBasis2, double *sigma_1, double *sigma_2,
                        double *sigma_date, int sigma_n, double lambda1,
                        double lambda2, double alpha, double gamma, double rho,
                        int nb_integration, double bound, double **output);

#endif