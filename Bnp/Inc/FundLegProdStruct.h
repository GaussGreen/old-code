#ifndef __FUND_LEG_STRUCT_H
#define __FUND_LEG_STRUCT_H

#define FUND_LEG_NCPN 512

/*	Structures and functions for the funding leg */

/*	Funding cpn */
typedef struct {
  long start_date;   /*	Coupon start date */
  double start_time; /*	Coupon start time */
  long pay_date;     /*	Coupon pay date */
  double pay_time;   /*	Coupon pay time */
  double cvg;        /*	Cvg */
  double cpn;        /*	Notional * (fwd_spr + margin) * cvg */
} funding_cpn, *FUNDING_CPN;

/*	Funding leg */
typedef struct {
  double notional; /*	Notional */
  int num_cpn;     /*	Number of coupons */
  int spot_lag;    /*	Spot Lag	*/
  /*	0..num_cpn-1 */
  funding_cpn *cpn;
} funding_leg, *FUNDING_LEG;

Err fill_funding_leg(
    /*	Coupons that started before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I      , 1: E */
    int spot_lag, double fund_not, int fund_ncpn, long *fund_fix,
    long *fund_start, long *fund_pay, char **fund_basis, double *fund_spr,
    double *fund_mrg, FUNDING_LEG fund_leg);

/*	Check dates consistency */
Err check_funding_leg(FUNDING_LEG fund_leg);

/*	Free */
Err free_funding_leg(FUNDING_LEG fund_leg);

Err FundLeg_RequestDfDates(funding_leg *fundleg, long date, int *n_dates,
                           long **pdates);

Err FundLeg_Payoff(funding_leg *fundleg, long today, long date, double *dfs,
                   double *payoff);

#endif