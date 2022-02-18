
#ifndef __CPD_CALIB_H__
#define __CPD_CALIB_H__

#define CALPRES 1.0e-08

void static_lgmsetupG_tauts(
    int    nlam,
    double lam_time[],
    double lam[],
    int    ncpn,       /*	Number of cash-flows */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_G[],    /*	Output: G at cash-flow dates
                                                       G(T) = (1.0 - exp (- lambda * T )) */
    int    nex,        /*	Number of exercise dates */
    double ex_time[],  /*	Exercise times */
    double ex_G[]);

void static_interpolate_zeta(
    int    nb_old_zeta,
    double old_zeta_time[],
    double old_zeta[],
    int    nlam,
    double lam_time[],
    double lam[],
    int    nb_new_zeta,
    double new_zeta_time[],
    double new_zeta[]);

char* lgmcalibzeta1F_tauts2(
    int    ncpn,       /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G[],    /*	G at cash-flow dates */
    int    nex,        /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int    ex_cpn[],   /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G[],      /*	G at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double mkt_vega[],
    double ex_zeta[], /*	Output: zetas */
    int    nlam,
    double lam_time[],
    double lam[],
    int skip_last, /*	If 1, the last option is disregarded and the forward volatility is flat from
                      option n-1 */
    double prec,
    int    vega_prec,
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps);  /*	Use Jumps */

char* lgmcalibzeta2F_tauts2(
    int    ncpn,       /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G1[],   /*	G1 at cash-flow dates */
    double cpn_G2[],   /*	G2 at cash-flow dates */
    int    nex,        /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int    ex_cpn[],   /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[],     /*	G2 at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double mkt_vega[],
    double ex_zeta[], /*	Output: zetas (1) */
    /*	Lambda, Alpha, gamma, rho */
    int    nlam,
    double lam_time[],
    double lam[],
    double alpha,
    double gamma,
    double rho,
    int skip_last, /*	If 1, the last option is disregarded and the forward volatility is flat from
                      option n-1 */
    double prec,
    int    vega_prec,
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps);  /*	Use Jumps */

double static_lgmcalcexpfact_tauts(double T1, double T2, int nlam, double lam_time[], double lam[]);

void static_lgmcalczeta2zeta12_tauts(
    int    n,       /*	Number of dates */
    double t[],     /*	Times */
    double zeta1[], /*	Zeta1 */
    /*	Lambda, Alpha, Beta, Rho */
    int    nlam,
    double lam_time[],
    double lam[],
    double alpha,
    double gamma,
    double rho,
    /*	Output */
    double zeta2[],
    double zeta12[]);

/*	Setup G2 function (2F) */
/*	Tau-TS enabled version */
void static_lgmsetupG2_tauts(
    int    nlam,
    double lam_time[],
    double lam[],
    double gamma,
    int    ncpn,       /*	Number of cash-flows */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_G2[],   /*	Output: G2 at cash-flow dates
                                                       G2(T) = (1.0 - exp (- lambda2 * T )) */
    int    nex,        /*	Number of exercise dates */
    double ex_time[],  /*	Exercise times */
    double ex_G2[]);

/*  Useful functions */
double static_lgmystar_stochvol(
    int    ncpn,  /*	Num coupons */
    double cpn[], /*	Discounted coupons */
    double beta[],
    double psi,
    int*   dir); /*	Output, whether iv is increasing in y
                                   1: decreasing, -1: increasing */

/*	Just a safer Gaussian */
double static_lgmsafenorm(double z);

/*	Value of European option on a general set of cash flows within LGM 1F */
double lgmopval1F(
    int    ncpn,    /*	Number of cash-flows */
    double cpn[],   /*	Discounted Cash-Flows */
    double cpn_G[], /*	G at cash-flow dates */
    double ex_zeta, /*	Zeta at exercise date */
    double ex_G);   /*	G at exercise date */

/*	Value of European option on a general set of cash flows within LGM Stoch Vol */
double lgmopval_stochvol(
    int    ncpn,   /*	Number of cash-flows */
    double cpn[],  /*	Discounted Cash-Flows */
    double beta[], /*	beta(T*, T) at cash-flow dates */
    double psi,    /*	Psi at exercise date */
    double driftf, /*  Std of f(t,T*) */
    double rho2,   /*	1 - rho * rho */
    double sqrho2);

/*	Value of European option on a general set of cash flows within LGM 2F */
double lgmopval2F(
    int    ncpn,      /*	Number of cash-flows */
    double cpn[],     /*	Discounted Cash-Flows */
    double cpn_G1[],  /*	G1 at cash-flow dates */
    double cpn_G2[],  /*	G2 at cash-flow dates */
    double ex_zeta1,  /*	Z1 at exercise date */
    double ex_zeta2,  /*	Z2 at exercise date */
    double ex_zeta12, /*	Z12 at exercise date */
    double ex_G1,     /*	G1 at exercise date */
    double ex_G2);    /*	G2 at exercise date */

/*	Value of European Swap option within LGM 1F */
double lgmsopval1F(
    int ncpn,       /*	Number of cash-flow dates, including
                                            start and end date */
    double df[],    /*	Df to cash flow dates */
    double cvg[],   /*	cvg from i-1 to i */
    double cpn_G[], /*	G at cash-flow dates */
    double ex_zeta, /*	Z at exercise date */
    double ex_G,    /*	G at exercise date */
    double strike); /*	Strike */

/*	Value of European Swap option within LGM 2F */
double lgmsopval2F(
    int ncpn,         /*	Number of cash-flow dates, including
                                              start and end date */
    double df[],      /*	Df to cash flow dates */
    double cvg[],     /*	cvg from i-1 to i */
    double cpn_G1[],  /*	G1 at cash-flow dates */
    double cpn_G2[],  /*	G2 at cash-flow dates */
    double ex_zeta1,  /*	Z1 at exercise date */
    double ex_zeta2,  /*	Z2 at exercise date */
    double ex_zeta12, /*	Z12 at exercise date */
    double ex_G1,     /*	G1 at exercise date */
    double ex_G2,     /*	G2 at exercise date */
    double strike);   /*	Strike */

/*	Value of European Cap within LGM 1F */
double lgmcapval1F(
    int ncpn,        /*	Number of cash-flow dates, including
                                             start and end date */
    double df[],     /*	Df to cash flow dates */
    double cvg[],    /*	cvg from i-1 to i */
    double cpn_G[],  /*	G at cash-flow dates */
    int    nex,      /*	Number of exercise dates */
    int    ex_cpn[], /*	For each exercise date, first coupon
                                             to be exercised */
    int ex_ncpn[],   /*	For each exercise date, number of coupons
                                             to be exercised */
    double ex_weight[],
    double ex_zeta[], /*	Z at exercise date */
    double ex_G[],    /*	G at exercise date */
    double strike[]); /*	Strikes */

/*	Value of European Cap within LGM 2F */
double lgmcapval2F(
    int ncpn,           /*	Number of cash-flow dates, including
                                                start and end date */
    double df[],        /*	Df to cash flow dates */
    double cvg[],       /*	cvg from i-1 to i */
    double cpn_G1[],    /*	G1 at cash-flow dates */
    double cpn_G2[],    /*	G2 at cash-flow dates */
    int    nex,         /*	Number of exercise dates */
    int    ex_cpn[],    /*	For each exercise date, first coupon
                                                to be exercised */
    int ex_ncpn[],      /*	For each exercise date, number of coupons
                                                to be exercised */
    double ex_weight[], /*	Weights on each caplet */
    double ex_zeta1[],  /*	Z1 at exercise date */
    double ex_zeta2[],  /*	Z2 at exercise date */
    double ex_zeta12[], /*	Z12 at exercise date */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[],     /*	G2 at exercise date */
    double strike[]);   /*	Strikes */

/*	Calibrate zeta to diagonal given G: 1F case */
char* lgmcalibzeta1F(
    int    ncpn,        /*	Total number of cash-flow dates */
    double cpn_time[],  /*	Cash-Flow times */
    double cpn_df[],    /*	Df to cash-flow dates */
    double cpn_cvg[],   /*	cvg from i-1 to i */
    double cpn_G[],     /*	G at cash-flow dates */
    int    nex,         /*	Total number of exercise dates */
    double ex_time[],   /*	Exercise times */
    int    ex_cpn[],    /*	Index of the first cash-flow to be exercised */
    double ex_G[],      /*	G at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas */
    double lambda,
    int skip_last,   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    double prec,     /*	Precision on primary instruments */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps);  /*	Allow vol term structure to jump */

char* lgmcalibzeta1F_tauts(
    int    ncpn,        /*	Total number of cash-flow dates */
    double cpn_time[],  /*	Cash-Flow times */
    double cpn_df[],    /*	Df to cash-flow dates */
    double cpn_cvg[],   /*	cvg from i-1 to i */
    double cpn_G[],     /*	G at cash-flow dates */
    int    nex,         /*	Total number of exercise dates */
    double ex_time[],   /*	Exercise times */
    int    ex_cpn[],    /*	Index of the first cash-flow to be exercised */
    double ex_G[],      /*	G at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas */
    int    nlam,
    double lam_time[],
    double lam[],
    int skip_last,   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps);  /*	Allow vol term structure to jump */

/*	Calibrate zeta to swaptions given G: 1F case */
/*	New version: not necessarily the diagonal */
char* lgmcalibzeta1F_2(
    int    ncpn,       /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G[],    /*	G at cash-flow dates */
    int    nex,        /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int    ex_cpn[],   /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G[],      /*	G at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas */
    double lambda,
    int    skip_last); /*	If 1, the last option is disregarded
                                      and the forward volatility is flat from option
                                      n-1 */

/*	Calibrate zeta to diagonal given G: 2F case */
char* lgmcalibzeta2F(
    int    ncpn,        /*	Total number of cash-flow dates */
    double cpn_time[],  /*	Cash-Flow times */
    double cpn_df[],    /*	Df to cash-flow dates */
    double cpn_cvg[],   /*	cvg from i-1 to i */
    double cpn_G1[],    /*	G1 at cash-flow dates */
    double cpn_G2[],    /*	G2 at cash-flow dates */
    int    nex,         /*	Total number of exercise dates */
    double ex_time[],   /*	Exercise times */
    int    ex_cpn[],    /*	Index of the first cash-flow to be exercised */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[],     /*	G2 at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas (1) */
    /*	Lambda, Alpha, gamma, rho */
    double lambda,
    double alpha,
    double gamma,
    double rho,
    int skip_last,   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    double prec,     /*	Precision on primary instruments */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps);  /*	Allow vol term structure to jump */

char* lgmcalibzeta2F_tauts(
    int    ncpn,        /*	Total number of cash-flow dates */
    double cpn_time[],  /*	Cash-Flow times */
    double cpn_df[],    /*	Df to cash-flow dates */
    double cpn_cvg[],   /*	cvg from i-1 to i */
    double cpn_G1[],    /*	G1 at cash-flow dates */
    double cpn_G2[],    /*	G2 at cash-flow dates */
    int    nex,         /*	Total number of exercise dates */
    double ex_time[],   /*	Exercise times */
    int    ex_cpn[],    /*	Index of the first cash-flow to be exercised */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[],     /*	G2 at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas (1) */
    /*	Lambda, Alpha, gamma, rho */
    int    nlam,
    double lam_time[],
    double lam[],
    double alpha,
    double gamma,
    double rho,
    int skip_last,   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps);  /*	Allow vol term structure to jump */

/*	Calibrate zeta to diagonal given G: 2F case */
/*	New version: not necessarily the diagonal */
char* lgmcalibzeta2F_2(
    int    ncpn,       /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G1[],   /*	G1 at cash-flow dates */
    double cpn_G2[],   /*	G2 at cash-flow dates */
    int    nex,        /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int    ex_cpn[],   /*	Index of the first cash-flow to be exercised */
    /*	The extra argument */
    int ex_endcpn[], /*	Index of the last cash-flow to be exercised */
    /* */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[],     /*	G2 at exercise date */
    double strike[],    /*	Strikes */
    double mkt_price[], /*	Market prices */
    double ex_zeta[],   /*	Output: zetas (1) */
    /*	Lambda, Alpha, gamma, rho */
    double lambda,
    double alpha,
    double gamma,
    double rho,
    int    skip_last); /*	If 1, the last option is disregarded
                                      and the forward volatility is flat from option
                                      n-1 */

/*	Calibrate and price cap given lambda */
char* lgmprcapgivenlambda(
    int    ncpn,         /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int    nex,          /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int    ex_cpn[],     /*	Index of the first cash-flow to be exercised */
    int    ex_sncpn[],   /*	Number of coupons in each caplet */
    double ex_lstrike[], /*	Strikes for diagonal */
    double ex_lprice[],  /*	Market prices for diagonal */
    double ex_sstrike[], /*	Strikes for cap */
    double ex_sweight[], /*	Weights on each caplet */
    double ex_zeta[],    /*	Output: zetas */
    double lambda,       /*	Lambda */
    int    one2F,        /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,
    int skip_last, /*	If 1, the last option is disregarded and the forward volatility is flat from
                      option n-1 */
    double  prec,  /*	Precision on primary instruments */
    double  min_fact,   /*	Maximum down jump on variance */
    double  max_fact,   /*	Maximum up jump on variance */
    int     use_jumps,  /*	Allow vol term structure to jump */
    int     price_cap,  /*	0: just calibrate */
    double* ex_sprice); /*	Cap price as output */

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
/*	Tau-TS enabled version
 *** no lambda calibration so far for this version *** */
char* lgmprcapgivenlambda_tauts(
    int    ncpn,         /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int    nex,          /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int    ex_cpn[],     /*	Index of the first cash-flow to be exercised */
    double ex_lstrike[], /*	Strikes */
    double ex_lprice[],  /*	Market prices */
    double ex_zeta[],    /*	Output: zetas */
    int    one2F,        /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    int    nlam, /*	Lambda TS: may NOT be changed in the process */
    double lam_time[],
    double lam[],
    double alpha,
    double gamma,
    double rho,
    int skip_last,   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps);  /*	Allow vol term structure to jump */

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
/*	Tau-TS enabled version + calibration */
char* lgmprcapgivenlambda_tauts2(
    int    ncpn,       /*	Total number of cash-flow dates */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    int    nex,        /*	Total number of exercise dates */
    double ex_time[],  /*	Exercise times */
    int    ex_cpn[],   /*	Index of the first cash-flow to be exercised */
    int    ex_lendcpn[],
    double ex_lstrike[], /*	Strikes */
    double ex_lprice[],  /*	Market prices */
    double ex_lvega[],   /*	Market vegas */
    int    ex_sendcpn[],
    double ex_sstrike[], /*	Strikes */
    double ex_zeta[],    /*	Output: zetas */
    int    one2F,        /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    int    nlam, /*	Lambda TS: may NOT be changed in the process */
    double lam_time[],
    double lam[],
    double alpha,
    double gamma,
    double rho,
    int skip_last, /*	If 1, the last option is disregarded and the forward volatility is flat from
                      option n-1 */
    double prec,
    int    vega_prec,
    double min_fact,        /*	Maximum down jump on variance */
    double max_fact,        /*	Maximum up jump on variance */
    int    use_jumps,       /*	Use Jumps */
    int    price_cap,       /*	0: just calibrate */
    double ex_spriceres[]); /*	Cap price as output */

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
/*	New version: calibrates not necessarily to digonal
 *** no lambda calibration so far for this version *** */
char* lgmprcapgivenlambda_2(
    int    ncpn,         /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int    nex,          /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int    ex_cpn[],     /*	Index of the first cash-flow to be exercised */
    int    ex_endcpn[],  /*	Index of the last cash-flow to be exercised */
    double ex_lstrike[], /*	Strikes */
    double ex_lprice[],  /*	Market prices */
    double ex_zeta[],    /*	Output: zetas */
    double lambda,       /*	Lambda: may NOT be changed in the process */
    int    one2F,        /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,
    int    skip_last); /*	If 1, the last option is disregarded
                                      and the forward volatility is flat from option
                                      n-1 */

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
char* lgmcalibzetalambda(
    int    ncpn,         /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int    nex,          /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int    ex_cpn[],     /*	Index of the first cash-flow to be exercised */
    int    ex_sncpn[],   /*	Number of coupons in each caplet */
    double ex_lstrike[], /*	Strikes for diagonal */
    double ex_lprice[],  /*	Market prices for diagonal */
    double ex_sstrike[], /*	Strikes for cap */
    double ex_sweight[], /*	Weights on each caplet */
    double ex_sprice,    /*	Market price for cap */
    double ex_zeta[],    /*	Output: zetas */
    int    fix_lambda,   /*	0: calib lambda to cap, 1: fix lambda calib
                                                 to diagonal */
    double* lambda,      /*	Lambda: may be changed in the process */
    int     one2F,       /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,
    int skip_last, /*	If 1, the last option is disregarded and the forward volatility is flat from
                      option n-1 */
    double long_prec,  /*	Precision on primary instruments */
    double short_prec, /*	Precision on secondary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int    use_jumps);    /*	Allow vol term structure to jump */

/* Levenberg parameters */
typedef struct
{
    int    use_moment;
    int    nb_moment;
    int    break_moment;
    int    nb_iter;
    int    vega_weight;
    int    freq_short;
    int    shift_freq;
    double precision;

    int use_new;

} diag_calib_lm_params, *DIAG_CALIB_LM_PARAMS;

/*	Set defaults */
void diag_calib_lm_params_set_default_param(DIAG_CALIB_LM_PARAMS param);

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
char* lgmcalibzetalambda_tauts(
    int    ncpn,         /*	Total number of cash-flow dates */
    double cpn_time[],   /*	Cash-Flow times */
    double cpn_df[],     /*	Df to cash-flow dates */
    double cpn_cvg[],    /*	cvg from i-1 to i */
    int    nex,          /*	Total number of exercise dates */
    double ex_time[],    /*	Exercise times */
    int    ex_cpn[],     /*	Index of the first cash-flow to be exercised */
    int    ex_lncpn[],   /*	Number of coupons in each caplet */
    double ex_lstrike[], /*	Strikes for diagonal */
    double ex_lprice[],  /*	Market prices for diagonal */
    int    ex_sncpn[],   /*	Number of coupons in each caplet */
    double ex_sstrike[], /*	Strikes for cap */
    double ex_sprice[],  /*	Market price for cap */
    double ex_svega[],   /*	Market vega for cap */
    double ex_zeta[],    /*	Output: zetas */
    int    fix_lambda,   /*	0: calib lambda to cap, 1: fix lambda calib
                                                 to diagonal */
    int    nlam,         /*	Lambda TS: may NOT be changed in the process */
    double lam_time[],
    double lam[],
    int    one2F, /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double               alpha,
    double               gamma,
    double               rho,
    int                  skip_last,
    DIAG_CALIB_LM_PARAMS lm_params); /*	If 1, the last option is disregarded
                                                                     and the forward volatility is
                                        flat from option n-1 */

/*	Structures that hold calibration instrument data */
typedef struct
{
    /* Primary instruments */
    int     num_inst;
    long*   exer_dates_long;
    long*   start_dates;
    long*   end_dates;
    double* long_strikes;
    double* market_prices_long;

    /* Secondary instruments */
    int     num_insts;
    long*   exer_dates_short;
    long*   start_datess;
    long*   end_datess;
    double* short_strikes;
    double* short_weights;
    double* market_prices_short;

    /* Smile instruments */
    int     num_inst_smile;
    long*   start_dates_smile;
    long*   end_dates_smile;
    double* alpha;
    double* rho;

} cpd_calib_inst_data, *CPD_CALIB_INST_DATA;

char* get_end_date(long ex_date, long struct_end_date, char* tenor, int theo_act, long* end_date);

void cpd_init_hermite_for_calib(int one2F);

/*	Calibrate lgm: main function */
char* cpd_calib_diagonal(
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */
    char* ref_rate_name,  /*	Name of the reference rate */
    char* (*get_cash_vol)(/*	Function to get cash vol from the market */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),
    double vol_shift,
    int    shift_type,    /*	0:	Additive
                                          1:	Multiplicative */
                          /*	If ex_date is NULL,
                          exercise dates will be generated 2bd before start */
    int num_ex_dates,     /*	Exercise dates,
                                                          all supposed to be on or after today */
    long* ex_date,        /*	Supposed to be sorted
                                                          NULL = 2bd before each coupon */
    long    end_date,     /*	End date for diagonal */
    double* long_strike,  /*	Diagonal swaption strikes
                                                                  NULL = ATM */
    double* short_strike, /*	Short swaption strikes
                                                                  NULL = ATM */
    int strike_type,      /*	0: ATM
                                                  1: CASH
                                                  2: SWAP
                                                  3: STD */
    double max_std_long,
    double max_std_short,
    char*  swaption_freq, /*	Frequency and basis of underlying swaptions */
    char*  swaption_basis,
    int    fix_lambda,    /*	0: calib lambda to cap, 1: fix lambda calib
                                                  to diagonal */
    int one_f_equi,       /*	1F equivalent flag:
                                                  if set to 1, then 2F lambda will calibrate
                                                  to the cap priced within calibrated 1F
                                                  with the given lambda */
    int skip_last,        /*	If 1, the last option is disregarded
                                                  and the forward volatility is flat from option
                                                  n-1 */
    double  long_prec,    /*	Precision on primary instruments */
    double  short_prec,   /*	Precision on secondary instruments */
    double  min_fact,     /*	Maximum down jump on variance */
    double  max_fact,     /*	Maximum up jump on variance */
    int     use_jumps,    /*	Allow vol term structure to jump */
    int     proba_weight, /*	Proba weighting for caplet */
    double* proba,

    double* lambda, /*	Lambda: may be changed in the process */
    int     one2F,  /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double   alpha,
    double   gamma,
    double   rho,
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data); /*	NULL = don't save calibration instrument data */

/*	Calibrate lgm: main function */
/*	Tau-TS enabled version
 *** no lambda calibration so far for this version *** */
char* cpd_calib_diagonal_tauts(
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */
    char* ref_rate_name,  /*	Name of the reference rate */
    char* (*get_cash_vol)(/*	Function to get cash vol from the market */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),
    double vol_shift,
    int    shift_type,   /*	0:	Additive
                                         1:	Multiplicative */
                         /*	If ex_date is NULL,
                         exercise dates will be generated 2bd before start */
    int num_ex_dates,    /*	Exercise dates,
                                                         all supposed to be on or after today */
    long* ex_date,       /*	Supposed to be sorted
                                                         NULL = 2bd before each coupon */
    long    end_date,    /*	End date for diagonal */
    double* long_strike, /*	Diagonal swaption strikes
                                                                 NULL = ATM */
    int strike_type,     /*	0: ATM
                                                 1: CASH
                                                 2: SWAP
                                                 3: STD */
    double max_std_long,
    char*  swaption_freq, /*	Frequency and basis of underlying swaptions */
    char*  swaption_basis,
    int skip_last,   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int    use_jumps, /*	Allow vol term structure to jump */

    int    nlam, /*	Lambda TS: may NOT be changed in the process */
    double lam_time[],
    double lam[],
    int    one2F, /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double   alpha,
    double   gamma,
    double   rho,
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data); /*	NULL = don't save calibration instrument data */

typedef enum
{
    NORMAL_VOL,
    LOGNORMAL_VOL,
    LGM_VOL

} DIAGCALIB_VOLTYPE;

char* diagcalib_interp_voltype(const char* constStr, DIAGCALIB_VOLTYPE* vol_type);

typedef enum
{
    ADDITIVE,
    MULTIPLICATIVE

} DIAGCALIB_SHIFTTYPE;

char* diagcalib_interp_shifttype(const char* constStr, DIAGCALIB_SHIFTTYPE* shift_type);

/*	Parameter structure */
typedef struct cpd_diag_calib_param_
{
    /*	Shifts */
    double              vol_shift;
    DIAGCALIB_VOLTYPE   vol_type;   /*	0: normal, 1: lognormal */
    DIAGCALIB_SHIFTTYPE shift_type; /*	0:	Additive
                                                            1:	Multiplicative */
    double lambda_shift;
    double lambda_min;
    double lambda_max;

    int transform_vol;

    /*	Strikes */
    int strike_type; /*	0: ATM
                                     1: CASH
                                     2: SWAP
                                     3: STD
                                     4: LONGSTD: same number of STD
                                             as the long instruments */
    /*	Exclusions */
    double max_std;
    double min_time; /*	Minimum time between 2 instruments */
    int skip_last;   /*	If 1, the last option is disregarded and the forward volatility is flat from
                        option n-1 */
    int    keep_first;     /*	Min time is not applied on "keep_first" first options */
    double min_calib_time; /*	Minimum time from today for an option to be calibrated */
    double min_fact;       /*	Maximum down jump on variance */
    double max_fact;       /*	Maximum up jump on variance */
    int    use_jumps;      /*	1: we allow jumps on vol, 0: we don't */

    /* For Newton */
    double precision;
    int    vega_prec;
    int    nb_iter_max;

    /* For Smile */
    double              smile_vol_shift;
    DIAGCALIB_VOLTYPE   smile_vol_type;
    DIAGCALIB_SHIFTTYPE smile_shift_type;
    int                 smile_strike_type;

    int    smile_use_jumps;
    double smile_precision;
    int    smile_nb_iter_max;

    /*	Fx */
    double fx_vol_shift;

} cpd_diag_calib_param, *CPD_DIAG_CALIB_PARAM;

/*	Set defaults */
void cpd_calib_set_default_param(CPD_DIAG_CALIB_PARAM param);

/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to digonal
 *** no lambda calibration so far for this version *** */
char* cpd_calib_diagonal_2(
    /*	Market */
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */
    char* ref_rate_name,  /*	Name of the reference rate */
    char* (*get_cash_vol)(/*	Function to get cash vol from the market */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),
    char* instr_freq, /*	Frequency and basis of instruments */
    char* instr_basis,
    /*	Structure */
    /*	If ex_date is NULL,
    exercise dates will be generated 2bd before start */
    int num_ex_dates,  /*	Exercise dates,
                                                       all supposed to be on or after today */
    long*  ex_date_,   /*	Supposed to be sorted */
    int*   cal_date,   /*	1: use ex_date as calibration date, 0: don't */
    char** end_tenor_, /*	Tenors of the underlying instruments
                                                               or "DIAG" */
    long    end_date,  /*	End date for diagonal */
    double* strike_,   /*	Strikes
                                               0: ATM */
    /*	Model */
    double lambda, /*	Lambda: may NOT be changed in the process */
    int    one2F,  /*	Number of factors */
    double alpha,  /*	Alpha, Gamma, Rho (2F only) */
    double gamma,
    double rho,
    /*	Output */
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data); /*	NULL = don't save calibration instrument data */

/*	Initialise calibration instrument data */
void cpd_init_calib_inst_data(CPD_CALIB_INST_DATA inst_data);

/*	Copy calibration instrument data */
void cpd_copy_calib_inst_data(CPD_CALIB_INST_DATA dest, CPD_CALIB_INST_DATA src);

/*	Free calibration instrument data */
void cpd_free_calib_inst_data(CPD_CALIB_INST_DATA inst_data);

/*	Calibrate a 3f model: main function */
char* cpd_calib_all(
    /*	Today */
    long today,
    /*	Get Cash Vol function */
    char* (*get_cash_vol)(/*	Function to get cash vol from the market */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),
    /*	Domestic market */
    char* dom_yc_name,        /*	Name of the yield curve */
    char* dom_vol_curve_name, /*	Name of the market vol curve */
    char* dom_ref_rate_name,  /*	Name of the reference rate */
    char* dom_instr_freq,     /*	Frequency and basis of instruments */
    char* dom_instr_basis,
    /*	Domestic Model */
    double dom_lambda, /*	Lambda: may NOT be changed in the process */
    /*	Foreign market */
    char* for_yc_name,        /*	Name of the yield curve */
    char* for_vol_curve_name, /*	Name of the market vol curve */
    char* for_ref_rate_name,  /*	Name of the reference rate */
    char* for_instr_freq,     /*	Frequency and basis of instruments */
    char* for_instr_basis,
    /*	Domestic Model */
    double for_lambda, /*	Lambda: may NOT be changed in the process */
    /*	Fx market */
    long*   fx_mkt_vol_date, /*	Option maturity dates */
    double* fx_mkt_vol,      /*	Option BS vol */
    int     num_fx_mkt_vol,  /*	Number of options */
    /*	Fx model */
    double* corr_times,
    double* correl_dom_for,
    double* correl_dom_fx,
    double* correl_for_fx,
    long    corr_n_times,
    /*	Structure */
    /*	If ex_date is NULL,
    exercise dates will be generated 2bd before start */
    int num_ex_dates,      /*	Exercise dates,
                                                           all supposed to be on or after today */
    long*   ex_date,       /*	Supposed to be sorted */
    int*    cal_date,      /*	1: use ex_date as calibration date, 0: don't */
    char**  dom_end_tenor, /*	Tenors of the underlying instruments or "DIAG" */
    char**  for_end_tenor,
    long    end_date,   /*	End date for diagonal */
    double* dom_strike, /*	Domestic strikes 0: ATM */
    double* for_strike, /*	Foreign strikes 0: ATM */
    /*	Output */
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** dom_sig,
    double** for_sig,
    int*     num_fx_vol,
    double** fx_vol_time,
    double** fx_vol,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA dom_inst_data,  /*	NULL = don't save calibration instrument data */
    CPD_CALIB_INST_DATA for_inst_data); /*	NULL = don't save calibration instrument data */

/*	Make underlying out of the results of the previous function */
char* cpd_calib_all_makeund(
    long    today,
    char*   dom_ccy,
    char*   dom_und_name,
    char*   dom_yc_name,
    double* dom_sig,
    double  dom_lambda,
    char*   for_ccy,
    char*   for_und_name,
    char*   for_yc_name,
    double* for_sig,
    double  for_lambda,
    double* corr_times,
    double* correl_dom_for,
    double* correl_dom_fx,
    double* correl_for_fx,
    long    corr_n_times,
    int     num_sig,
    double* sig_time,
    int     num_fx_vol,
    double  spot_fx,
    double* fx_vol_time,
    double* fx_vol,
    /*	Output only */
    char fx_und_name[]);

/*	Calibrate lgm: main function */
/*	New version: calibrates not necessarily to digonal
                with lambda calibration */
char* cpd_calib_diagonal_3(
    /*	Market */
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */
    char* ref_rate_name,  /*	Name of the reference rate */
    char* (*get_cash_vol)(/*	Function to get cash vol from the market */
                          char*   vol_curve_name,
                          double  start_date,
                          double  end_date,
                          double  cash_strike,
                          int     zero,
                          char*   ref_rate_name,
                          double* vol,
                          double* power),
    char* instr_freq, /*	Frequency and basis of instruments */
    char* instr_basis,
    /*	If ex_date is NULL,
    exercise dates will be generated 2bd before start */
    /*	Structure */
    int num_ex_dates,   /*	Exercise dates,
                                                        all supposed to be on or after today */
    long*  ex_date_,    /*	Supposed to be sorted */
    int*   cal_date,    /*	1: use ex_date as calibration date, 0: don't */
    char** end_tenorl_, /*	Tenors of the underlying instruments
                                                                or "DIAG" */
    char** end_tenors_, /*	Tenors of the underlying instruments
                                                                or "DIAG" */
    long    end_date,   /*	End date for diagonal */
    double* strikel_,   /*	Strikes
                                                0: ATM */
    double* strikes_,   /*	Strikes
                                                0: ATM */
    /*	Model */
    int    fix_lambda,
    int    nlam, /*	Lambda TS: may NOT be changed in the process */
    double lam_time[],
    double lam[],
    int    one2F, /*	Number of factors */
    double alpha, /*	Alpha, Gamma, Rho (2F only) */
    double gamma,
    double rho,
    /*	Output */
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param,
    DIAG_CALIB_LM_PARAMS lm_params,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data); /*	NULL = don't save calibration instrument data */

double export_lgmcalcexpfact_tauts(
    double T1, double T2, int nlam, double lam_time[], double lam[], double gamma);

void export_lgmcalczeta1_tauts(
    int    nsig,       /*	Number of Sigmas */
    double sig_time[], /*	Sigma Times */
    double sig[],      /*	Sigmas */
    /*	Lambda, Alpha, Beta, Rho */
    int    nlam,
    double lam_time[],
    double lam[], /*	Lambdas */
    double alpha,
    double gamma,
    double rho,
    /*	Output */
    int    nzeta,
    double zeta_t[], /*	Zeta times */
    double zeta[]);  /*	Output */

void export_lgmcalczeta2zeta12(
    int    n,       /*	Number of dates */
    double t[],     /*	Times */
    double zeta1[], /*	Zeta1 */
    /*	Lambda, Alpha, Beta, Rho */
    double lambda,
    double alpha,
    double gamma,
    double rho,
    /*	Output */
    double zeta2[],
    double zeta12[]);

void export_lgmcalczeta2zeta12_ts(
    int    n,       /*	Number of dates */
    double t[],     /*	Times */
    double zeta1[], /*	Zeta1 */
    /*	Lambda, Alpha, Beta, Rho */
    int    nlambda,
    double lambda_time[],
    double lambda[],
    double alpha,
    double gamma,
    double rho,
    /*	Output */
    double zeta2[],
    double zeta12[]);

void export_lgmsetupG(
    double lambda,
    int    ncpn,       /*	Number of cash-flows */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_G[],    /*	Output: G at cash-flow dates
                                                       G(T) = (1.0 - exp (- lambda * T )) */
    int    nex,        /*	Number of exercise dates */
    double ex_time[],  /*	Exercise times */
    double ex_G[]);    /*	Output: G at exercise dates */

/// added by albert wang with lambda ts
void export_lgmsetupG_ts(
    int    nlambda,
    double lambda_time[],
    double lambda[],

    int    ncpn,       /*	Number of cash-flows */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_G[],    /*	Output: G at cash-flow dates
                                                       G(T) = (1.0 - exp (- lambda * T )) */
    int    nex,        /*	Number of exercise dates */
    double ex_time[],  /*	Exercise times */
    double ex_G[]);    /*	Output: G at exercise dates */

void export_lgmsetupG2(
    double lambda,
    double gamma,
    int    ncpn,       /*	Number of cash-flows */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_G2[],   /*	Output: G2 at cash-flow dates
                                                       G2(T) = (1.0 - exp (- lambda2 * T )) */
    int    nex,        /*	Number of exercise dates */
    double ex_time[],  /*	Exercise times */
    double ex_G2[]);   /*	Output: G2 at exercise dates */

void export_lgmsetupG2_ts(
    int    nlambda,
    double lambda_time[],
    double lambda[],

    double gamma,
    int    ncpn,       /*	Number of cash-flows */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_G2[],   /*	Output: G2 at cash-flow dates
                                                       G2(T) = (1.0 - exp (- lambda2 * T )) */
    int    nex,        /*	Number of exercise dates */
    double ex_time[],  /*	Exercise times */
    double ex_G2[]);   /*	Output: G2 at exercise dates */

//// integrate lambda function from 0 to cpn_time[i]
double _average_lambda_(double dEndTime, int nlambda, double lambda_time[], double lambda[]);

#endif