/*          LGM autocal
Evaluate a midatlantic or american swaption after matching to market.
For american swaptions:
    Evaluates the deal as a midatlantic with exercises every month
    apart  , re-evaluates the deal as a midatlantic with exercises every
    two months apart  , and then extrapolates to get the value for exercises
    every day. The value of the american swaption is *LGMamervalue.

For midatlantics:
 1. Find the LGM parameters by matching swaptions to market (BS) price
    A. Fit model parameters by matching a series of long swaptions and a series
    of short (1 year) swaptions

    B. Swaptions fitted:
    "Long" swaptions: Vanilla swaptions which start every six months in
semiannual markets (USD) and every year in annual markets (Europe)  , beginning
with the initial start date of the mid-atlantic  , All these "long" reference
swaptions end at the approximate end date of mid-atlantic "Short" swaptions:
Vanilla swaptions of one year duration starting every six months in semiannual
markets (USD) and every year in annual markets (Europe)  , beginning with the
    initial start date of the mid-atlantic.
                Caplets fitted:	Caplets that start every three or six months
depending on the market
        ***If fix_tau_flag!=1  , then the code fits the long swaption and  ,
whenever possible  , the short swaptions. It is very rare that a deal cannot
have all short swaptions fitted exactly (4 out of 150) as well as all the long
swaptions.
        ***If fix_tau_flag==1 and usecaps!=1  , then only the long swaptions are
fitted
        ***If fix_tau_flag==1 and usecaps==1  , then only the caplets are fitted
       Strikes of the long reference swaptions:
    The strikes are chosen so that the swaption has the same ratio of PV(fixed
leg)/PV(floating leg) as the swap that would be recieved if the mid-atlantic was
exercised at the starting date of the reference swaption. Strikes for the short
reference swaptions: The strike is the "effective interest rate" of the
mid-atlantic for that year  , including adjustments for exercise premiums  ,
etc. Strikes for the caplets: The strikes are chosen like the short reference
swaption strikes


    C. Numerical method used: Various global & local Newton routines for
fitting. Whenever fit  , swaption prices are matched to relative error ofless
than 1.e-8.

 2. Uses the fitted parameters and evaluates the midatlantic by convolution
method. Numerical accuracy of the convolution is better than a relative error
of 1.E-5. The midatlantic value is *LGMvalue.

 3. If the update_ts_flag is on  , the code will determine the sigma and tau
values which generate the integrals (zeta[j] and g[j] used by the code and put
them in the market denoted by SrtCurvePtr. IMPORTANT NOTE: Any method using LGM
calibration should NOT interpolate on sigmas and tau values. This is an unstable
procedure  , and may lead to inaccurate prices (it hasn't yet  , however.)
    Should interpolate on the zeta(t) and G(t) values instead.
*/

/*
Payment rules implemented:
 For mid-atlantics (LGMautocal)  , the payment rule is assumed to be:

    If the deal is exercised at It_ex[j]  , which has a settlement ("as of")
date It_start[j]  , then the first fixed leg payment is fixpymnt[i]-Iredpay[j]
paid on tpay[i]  , and the other fixed leg payments are fixpymnt[i+1] on
tpay[i+1]  , ...  , fixpymnt[npay-1] on tpay[npay-1]. Here coupon i is the first
coupon whose pay date tpay[i] is strictly after the settlement date It_start[j].

    Note that the final payment should include the notional.

    The floating payments  , fees  , etc. are assumed to be equivalent to a
payment of Iprem[j] paid on the settlement date It_start[j]. Note that this
includes the notional  , and should be equal to the notional apart from fees and
adjustments.

    Any alterations in the payment schedule due to e.g. intraperiod starts
should be accounted for in Iredpay[j] (reduction in first fixed coupon payment
upon exercise j) and by adjusting the premium Iprem[j] before the code is
called. Any alterations in the floating side due to  , e.g. the floating pay
date being several days after the end date of the floating period need to be
accounted for in Iprem[j] before the code is called.

 For American swaptions THE NOTIONAL IS ASSUMED TO BE 1. The payment rules are:

    First possible exercise date is the latter of tfirst_ex and either today
(calculation date) if end_of_day_flg is off or the next business day if
end_day_flg is on.

    If the deal is exercised at T  , the theoretical start (settlement or "as
of") date S is lag_exer calendar or business days after T. The actual start date
is obtained from the theoretical start date using the Business Day Rule
conv_start.

    Let i be the first fixed coupon whose end date tfix_end[i] is strictly after
S. If early_flag_fix is off  , it is assumed that the first fixed coupon will be
reduced by the interest accrued between the coupon's start date tfix_start[i]
and the settlement date S (so that the first payment corresponds to the interest
for the period between S and tfix_end[i]). If early_fix_flag is on  , it is
assumed that the interest accrued between tfix_start[i] and S is added to the
floating leg at S  , effectively increasing the value of the floating leg  , and
the entire fixed coupon is paid at i.

        The last fixed rate payment full_pay_fix[npay-1] should include the
notional (which should be 1).

    On the floating side  , suppose that t_flt[k] is the first floating pay date
after the settlement date S  , and let coupon i be the first fixed coupon whose
pay date is strictly after S. The code assumes that all floating leg payments
that are made AFTER t_flt[k] plus any exercise premiums are equivalent to a
payment of Inprem[i] made at t_flt[k]. So like the mid-ats  , Inprem[i] should
include the notional plus (for receivers) or minus (for payers) any exercise
premium.

    Unless the deal is exercised exactly at the end of a floating period  , the
code will adjust the value of the floating leg for the accrued intereston the
floating leg. If the early_flag_flt is OFF  , then the code assumes that the
interest  , rflt*CVG(S  ,t_flt[k]) that accrues between settlement and the next
        floating pay date t_flt[k] is paid on t_flt[k] in addition to Inprem[i].
If the early_flag_flt is ON then the code assumes that rflt*CVG(t_flt[k-1]  ,S)
the interest that accrues between t_flt[k-1] and S is paid at settlement (goes
to reduce the floating leg's value) and that the full floating coupon
        rflt*CVG(t_flt[k-1]  ,t_flt[k]) is paid at t_flt[k] in addition to
Inprem[i].

    In the above  , rflt is the floating rate for the period from t_flt[k-1] to
t_flt[k]. If t_flt[k-1] <= tnow < t_flt[k]  , so that today is in the interval
, then the code will use the input argument rflt_current for this already-fixed
interest rate. If rflt_current is zero or negative  , or if tnow does not fall
in the interval  , then the code will calculate rflt from the current discount
        curve.
        */

/******** WARNINGS ***********
1) The payment rules must be repected when entering deals with weird frequencies
, long payment or start lags  , etc. 2) If the BS vols are bumped one at a time
for vega calculations  , best results are obtained if the vols are bumped by a
twentieth to a fiftieth of a per cent or less. If the BS vols are all bumped
   simultaneously  , bumping by a per cent is easily OK.
3) The code hasn't been tested as rigorously for "European style" markets in
which swaps have annual coupons instead of semi-annual
*******************************/

/* EXTERNAL ROUTINES USED
double swp_disc_mkt(Date t1  , Date t2  , SrtCurvePtr yldcrv  , int 1)
                           discount factor between t1 and t2 for mkt
double coverage(Date t1  , Date t2  , BasisCode basis)
                           coverage from t1 to t2 according to basis
String GetVol(Date t1  , Date t2  , double R  , double *sig)
                           get market Black-Scholes volatility (sig) for
                           swaption with exercise date t1 and end date t1+t2 */

/* Headers and prototypes */
#include "SRT_H_ALL.H> 
#include "SrtAccess.h"
#include "math.h"
#include "srt_h_lgmautocal.h"
#include "srt_h_lgmtypes.h"

static void set_lgm_volcurve(const char *cMarketId, const char *cVolCurveName);

static String mid_atl_fit_market(String (*GetVol)(),
                                 SrtDiffusionType srt_vol_type, long n,
                                 long nlast, int icmpnd,
                                 SrtCompounding compounding, SrtBasisCode basis,
                                 int spot_lag,
                                 SrtLgmRefSwptnData *lgmRefSwptnData);

static String mid_atl_fit_market_fix_tau(
    String (*GetVol)(), SrtDiffusionType srt_vol_type, long n, long nlast,
    int icmpnd, SrtCompounding compounding, SrtBasisCode basis, int usecaps,
    int spot_lag, SrtLgmRefSwptnData *lgmRefSwptnData);

static String Newt_Semi(long j, long n, double *philgptr, double *phishptr);

static String Newt_Annual(long j, long n, double *philgptr, double *phishptr);

static String pvcashflow(Date tnow, String ycName, long nex, Date t_ex[],
                         Date t_start[], long npay, Date tfix_start[],
                         Date tfix_end[], Date tpay[], double fixpymnt[],
                         long ifirst[], double pvpay[], double redpay[]);

static double comp_intrin_val(long nex, long ifirst[], double redpay[],
                              long npay, double pvpay[]);

static void choose_dates(Date tnow, Date tstart, long *n, long npay,
                         Date tpay[], long nex, Date t_ex[], int icmpnd,
                         long num_months, BusDayConv natural_conv);

static void smooth_strikes(long nlast, double smoothing);

static String Olivier_lg(long n, long nlast, int icmpnd, long nex, long npay,
                         long ifirst[], Date t_start[], double pvpay[],
                         double redpay[]);

static String Olivier_sh(long n, long nlast, long nex, long npay, long ifirst[],
                         Date t_start[], double pvpay[], double redpay[],
                         int icmpnd, int fix_tau_flag);

static String Findsigtau(long n, Date tnow, long spot_lag);

static String Findexerbdry(Date tnow, long nex, Date t_ex[], Date tend,
                           double phi_cr[], long n, String ycName, int icmpnd,
                           BasisCode natural_basis, BusDayConv natural_conv,
                           int natural_spot_lag,
                           SrtLgmExerBdryData *lgmExerBdryData);

static double Gnew(double th);

static double Gnew_der(double th);

static double Rpar_jn(double phi);

static double Rbs_par_jn(void);

static double Vbs_jn(double strike, double vol, SrtDiffusionType srt_vol_type,
                     int spot_lag);

static double pay_payer(double phi);

static double pay_payer_der(double phi);

static void Vjn_all(double *val, double *der1, double *der2);

static String R_restrictions(String (*GetVol)(), SrtDiffusionType srt_vol_type,
                             long n, long nlast, int icmpnd,
                             SrtCompounding compounding, SrtBasisCode basis,
                             int fix_tau_flag, int usecaps);

static double Vjn(double par, int whichpar);

static double Vjn_der(double par, int whichpar);

static double Vjn1(double par, int whichpar);

static double Vjn1_der(double par, int whichpar);

static String mid_atl_eval(Date tnow, long n, SrtCurvePtr crv, long nex,
                           Date t_ex[], Date t_start[], long npay, Date tpay[],
                           double pvpay[], double redpay[], long ifirst[],
                           double *LGMvalue, double phi_cr[]);

static double iterp_integrand(double exk, double rik, long ik, double Q1[],
                              long nphi, double v_cr, double exk_cr,
                              double rik_cr, long di);

static void interpolate_par(Date tnow, long n, long nex, Date t_ex[],
                            double zeta_ex[], Date t_start[], double G_start[],
                            long npay, Date tpay[], double G_pay[]);

static void genweights(long num_q, double wid, double wgt[]);

static double EXPAY(double phi, long npay, long i0, double pvpay[],
                    double redpayj, double strike, double G_pay[], double Gst,
                    double zetaj);

static void init_point(double *x, int *check, double target, long i0, long npay,
                       double pvpay[], double redpayj, double G_pay[],
                       double Gst, double zetaj);

static void ScNewt(double *x, int *check, double target, double (*func1)(),
                   double (*deriv1)());

static void ScSlowNewt(double *x, double xmin, double xmax, int *check,
                       double target, int gzfind, double (*func1)(),
                       double (*deriv1)());

static void FreeMost(void *, void *, void *, void *, void *, void *);

static void FreeMore(void *, void *, void *, void *, void *, void *, void *,
                     void *, void *, void *);

static void FreeLots(void *, void *, void *, void *, void *, void *, void *,
                     void *, void *);

static void FreeRest(void *, void *);

static Err get_ccy_defaults(SrtCurvePtr crvptr, int fix_tau_flag, int usecaps,
                            BasisCode *basis, SrtCompounding *compnd,
                            BusDayConv *conv, int *spot_lag);
/* S. R. may-07-96 */
/* MACROS */

#define ERROR_CHECK                                                            \
  if (error != NULL) {                                                         \
    FreeMost(t_ex, t_start, ifirst, pvpay, redpay, phi_cr);                    \
    return (error);                                                            \
  }

#define ERROR_CHECK1                                                           \
  if (error != NULL) {                                                         \
    FreeMore(zeta_ex, G_start, G_pay, wgt, wgt2, Q0, Q1, eff_pay, postfac,     \
             effred);                                                          \
    return (error);                                                            \
  }

#define ERROR_CHECK2                                                           \
  if (error != NULL) {                                                         \
    FreeLots(tnatst, paydate, zeta_ex, G_pay, G_st, B_st, aB_pay, aB_first,    \
             ifirstpay);                                                       \
    srt_free(exerBdryArr);                                                     \
    return (error);                                                            \
  }

/* Definitions */

#define MAX_NUM_DATES                                                          \
  500 /* allow only midatlantics with fewer                                    \
         than 241 exer. and pay dates       */

#define MIN_EX_INT                                                             \
  10 /* ignore exercise dates closer than                                      \
         10 days apart                     */

#define MAX_REF_DATES                                                          \
  500 /* allow only midatlantics shorter than                                  \
         20 years from first exercise to end                                   \
         date                               */

#define FAR_OUT                                                                \
  2.5 /* strikes of swaptions used to match                                    \
         market limited to FAR_OUT std dev.s                                   \
         from par  , regardless of mid-atlantic                                \
         strike                             */

#define WIDTH                                                                  \
  3 /* integration domain runs from - WIDTH                                    \
        to + WIDTH std deviations in                                           \
        mid-atlantic evaluator             */

#define NDIV                                                                   \
  360 /* number of points used to discretize                                   \
         mid-atlantic. Should be divisible by 20 */

#define NQMAX                                                                  \
  48 /* number of points used in Gaussian                                      \
        kernel                             */

#define STENCIL                                                                \
  4 /* kernel runs from -STENCIL to +STENCIL                                   \
        std deviations                     */

/* variables needed to communicate with subroutines */
static long jst, nst; /* start (j) and end(n) of active reference swaption  */
static double Ract, phiact;   /* rate & phi of active ref. swaption   */
static Date s[MAX_REF_DATES]; /* dates s[j] of reference swaptions */
static double
    sy[MAX_REF_DATES]; /* time from tnow (in years) to dates of ref. swptns  */
static double B[MAX_REF_DATES]; /* discount factors from tnow to s[j] */
static double
    aB[MAX_REF_DATES]; /* products of discount factors and coverages    */
static double
    Ract_sh[MAX_REF_DATES]; /* fixed rates used for short ref. swaptions     */
static double
    Ract_lg[MAX_REF_DATES]; /* fixed rates used for long ref. swaptions      */
static double R_theor_sh[MAX_REF_DATES]; /* rates to be used if they aren't */
static double R_theor_lg[MAX_REF_DATES]; /*       too far from par */
static double Rmin_sh[MAX_REF_DATES];    /* minimum rates to be used for    */
static double Rmin_lg[MAX_REF_DATES];    /*     reference swaptions    */
static double Rmax_sh[MAX_REF_DATES];    /* maximum rates to be used for    */
static double Rmax_lg[MAX_REF_DATES];    /*     reference swaptions    */
static double
    prem[MAX_NUM_DATES]; /* premiums [j] paid if mid_atl exercised at sj */
static double
    origprem[MAX_NUM_DATES]; /* premiums [j] paid if mid_atl exercised at sj */
static int pr; /* 1 for payer  , 0 for reciever                   */
static double
    v_bs_sh[MAX_REF_DATES]; /* BS prices for short reference swaptions       */
static double
    v_bs_lg[MAX_REF_DATES]; /* BS prices for long reference swaptions        */

/* model parameters */
static double
    g[MAX_REF_DATES]; /* g[j] is integral of exp(-lambda*t) from tj to tn   */
static double zeta[MAX_REF_DATES]; /* zeta[j] is integral from zero to tj of
                                      sigma stuff */
static double kappa; /* FRA std dev. (short ref. swptns in annual markets  */
static double overall_tau; /* default tau value (essentially irrelevent) */

/* variables used for convenience */
static double Gj1, Gj2, sq2pi; /* miscellaneous constants */
static String error; /* error control                                      */

/* definitions to transfer sig's and tau's to market structure */
static long num_sig, num_tau; /* number of sigma and tau values               */
static double sig_val[MAX_NUM_DATES];  /* LGM model parameter sigma[j]  */
static double sig_date[MAX_NUM_DATES]; /*     at date sig_date[j] */
static double tau_val[MAX_NUM_DATES];  /* LGM model parameter tau[j]  */
static double tau_date[MAX_NUM_DATES]; /*     at date tau_date[j] */

/* static ptr to yc */
static String s_ycname = NULL;

/* Bug Fix A (BFA): Slightly changes the value of the floating leg if we are
currently in a floating period in which the interest rate has already been set.
The current rate is a new input parameter rflt_current; the program ignores this
input if it is zero or negative. To implement  , need to uncomment and delete
some lines */

/* Bug Fix B (BFB): Correctly accounts for the accrual of the floating rate
 * branch. Implemented */
/* MAIN PROGRAM to evaluate american swaptions */
/* NOTE: Notional must be 1   */
String LGMamerican(
    long nfix, /* number of fixed periods                                  */
    Date tfix_start[], /* [0  ,1  ,...  ,nfix-1]  start dates for fixed coupons
                        */
    Date tfix_end[],   /* [0  ,1  ,...  ,nfix-1]  end dates for fixed coupons   */
    Date tfix_pay[],   /* [0  ,1  ,...  ,nfix-1]  pay dates for fixed coupons   */
    double full_pay_fix[], /* [0  ,1  ,...  ,nfix-1]  fixed coupons (with
                              notional at end)   */
    double Inprem[],    /* [0  ,1  ,...  ,nfix-1]  premium for exer. at j (with
                           notional) */
    char *fix_basis,    /* basis for fixed coupons    */
    int early_flag_fix, /* 0=subtract accrual from frst pymnt; 1=add accrual to
                           fee */
    long nflt,    /* number of floating dates                                 */
    Date t_flt[], /* [0  ,1  ,...  ,nflt-1] all flting dtes (all strt dtes plus
                     final end dte) */
    /* BFA: remove this line and uncomment the next line to implement */
    /* double   rflt_current */ /* current floating rate (if fixed); ignored if
                                   non-positive   */
    char *flt_basis,            /* basis for floating coupons            */
    int early_flag_flt, /* 0 = subtrct accrual from frst pymnt; 1 = add accrual
                           to prem */
    Date t_first_ex,    /* first exerercise date    */
    int lag_exer_start, /* days between exercise and start */
    int cal_or_bus,     /* 0 = cal. days  , 1 = bus. days for lag_exer_start     */
    BusDayConv conv_start, /* business day convention for start */
    char *pay_rec_str,     /* SRT_RECEIVER or SRT_PAYER     */
    String yc_name,        /* Yield curve name    */
    int end_of_day_flag,   /*      */
    int fix_tau_flag,      /* 1 means calibrate with fixed tau */
    int usecaps,   /* if fixed tau  , then 1=use caplets to calibrate  , 0=use
                      long swptns */
    double In_tau, /* if fixed tau  , use this value for tau (in years) */
    String (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to get swaption vols */
    char *char_vol_type,        /*       */
    int update_ts_flag,   /* 1=cmpute new sigmas and taus & store in yldcrv;
                             0=don't bother */
    int find_exer_bdry,   /* 1=find swap rates at exercise bocrvary; 0=don't
                             bother       */
    String outfile,       /* output log file (unused)       */
    double *LGMamervalue, /* answer */
    SrtLgmExerBdryData
        *lgmExerBdryData) /* ptr to exercise boundary data structure         */
{
  /* declarations */
  Date t0, t1, tnow, t0_ex, t0_start;
  long i, j, k, mdiv, mdivst, count, countmax;
  long j0, In_ex, npay;
  long mdivar[2];
  double accrual, rflt, rdivst;
  Date It_ex[MAX_NUM_DATES], It_start[MAX_NUM_DATES];
  double Iprem[MAX_NUM_DATES], Iredpay[MAX_NUM_DATES], delta_tex[2];
  double LGMvalue[2]; /* mid-atlantic values */
  int end_of_day, update_flag, find_exer_flag;
  SrtReceiverType srt_p_r_flag;
  SrtDiffusionType srt_vol_type;
  SrtBasisCode srt_fix_basis;
  SrtBasisCode srt_float_basis;
  SrtErr srt_err;

  /* BFA: delete this line and the next to implement BFA */
  double rflt_current = 0.0;

  /* initialize */
  SrtCurvePtr yldcrv = NULL;
  yldcrv = lookup_curve(yc_name);

  error = NULL;
  full_pay_fix[nfix - 1] =
      full_pay_fix[nfix - 1] - 1; /* rem. notional from fixed leg */

  /* Transform char_vol_type into SrtDiffusionType */
  interp_diffusion_type(char_vol_type, &srt_vol_type);

  /* Transform fixed and floating into basis code */
  srt_err = interp_basis(fix_basis, &srt_fix_basis);
  if (srt_err)
    return (srt_err);

  srt_err = interp_basis(flt_basis, &srt_float_basis);
  if (srt_err)
    return (srt_err);

  /* Pay-Receive flag */
  srt_err = interp_rec_pay(pay_rec_str, &srt_p_r_flag);
  if (srt_err)
    return (srt_err);

  /* get earliest exercise date t0_ex and earliest start (as of) date t0_start
   */
  tnow = get_clcndate_from_yldcrv(yldcrv);
  end_of_day = (int)(end_of_day_flag ? 1 : 0);
  t0_ex = tnow;
  if (end_of_day == 1)
    t0_ex =
        add_unit(tnow, 1, SRT_DAY, SUCCEEDING); /* PH 2/28/98 trivial bug fix */
  t0_ex = (t0_ex > t_first_ex ? t0_ex : t_first_ex);
  if (cal_or_bus == 0)
    t0_start = add_unit(t0_ex + lag_exer_start - 1, 1, SRT_DAY, conv_start);
  else
    t0_start = add_unit(t0_ex, lag_exer_start, SRT_BDAY, conv_start);
  if (t0_start < tnow)
    t0_start = tnow;
  if (t0_start < t0_ex)
    t0_start = t0_ex; /* PH 2/28/98 potential bug fix */

  /* check inputs */
  /* PH 2/28/98 potential bug fix */
  if ((nfix < 1) || (nfix > MAX_NUM_DATES) || (nflt < 1) ||
      (nflt > MAX_NUM_DATES)) {
    error = "too many or too few dates";
    return (error);
  }

  if (t0_start >= (tfix_end[nfix - 1] - 2)) {
    error = "no exer. dates left";
    return (error);
  }

  for (i = 0; i < (nfix - 1) && tfix_start[i] < tfix_end[i] &&
              tfix_start[i] < tfix_start[i + 1] &&
              tfix_end[i] < tfix_end[i + 1] && tfix_pay[i] < tfix_pay[i + 1];
       i++)
    ;
  if (tfix_start[nfix - 1] >= tfix_end[nfix - 1])
    i = nfix - 2; /* PH 2/28/98 potential bug fix */
  for (j = 0; j < (nflt - 1) && t_flt[j] < t_flt[j + 1]; j++)
    ;

  if (i < (nfix - 1) || j < (nflt - 1)) {
    error = "dates are wrong";
    return (error);
  } else
    ; /* input looks OK */

  /* Get midatlantic value for countmax different sequences of dates */
  countmax = 2;
  if (t0_start >= (tfix_end[nfix - 1] - 30))
    countmax = 1; /* only do every month */
  else
    ; /* do every month & every other month */

  /* subdivide pay periods into 1 mo and 2mo subdivisions  , but ensure that the
     first exercise start date and all fixed start dates are included */

  rdivst = (tfix_end[nfix - 1] - tfix_start[0]) /* months per fixed period */
           * 12. / 365.25 / ((double)nfix);
  mdivst = DTOL(rdivst + .5);

  if (mdivst <= 1) /* only do every month if fixed dates are */
    countmax = 1;  /* less than two months apart             */

  mdivar[0] = mdivst;
  if (mdivar[0] < 1)
    mdivar[0] = 1;
  delta_tex[0] = rdivst / ((double)mdivar[0]);

  if (countmax >= 2) {
    mdivar[1] = DTOL(((double)mdivst) / 2. + .1);
    delta_tex[1] = rdivst / ((double)mdivar[1]);
  }

  /* evaluate the deal as a midatlantic with mdivar[count-1]
     subdivisions per fixed pay period */

  for (count = 1; count <= countmax; count++) {
    mdiv = mdivar[count - 1];
    In_ex = nfix * mdiv;
    if (In_ex > MAX_NUM_DATES) /* PH 2/28/98 potential bug fix */
    {
      error = "too many dates";
      return (error);
    }

    /* subdivide fixed pay periods into mdiv subintervals */
    for (i = 0; i < nfix; i++) {
      t0 = tfix_start[i];
      t1 = tfix_end[i];
      for (k = 0; k < mdiv; k++) {
        j = i * mdiv + k;
        It_start[j] =
            (Date)(t0 + ((double)(t1 - t0)) * ((double)k) / ((double)mdiv) +
                   0.5);
        if (cal_or_bus == 0)
          It_ex[j] = add_unit(It_start[j] - lag_exer_start + 1, -1, SRT_DAY,
                              conv_start);
        else
          It_ex[j] =
              add_unit(It_start[j], -lag_exer_start, SRT_BDAY, conv_start);
        if (It_ex[j] > It_start[j])
          It_ex[j] = It_start[j];
      } /* end k loop */
    }   /* end i loop */

    /* slip first exercise date in */
    j = 0;
    while (j < (In_ex - 1) && It_start[j] <= t0_start)
      j++;
    if (j == 0) /* first "as of" date is before first coupon start date */
    {
      j0 = 0;
      It_ex[0] = ((t0_ex > It_ex[0]) ? t0_ex : It_ex[0]);
    } else if (It_start[j] <= t0_start) {
      j0 = 0; /* first "as of" date is within last coupon period */
      It_start[0] = t0_start;
      It_ex[0] = t0_ex;
      In_ex = 1;
    } else {
      j0 = j - 1;
      if (It_start[j] < (t0_start + 14))
        j0 = j; /* use first "as of" date in place of coupon start date */
      else
        ; /* else append "as of" date before coupon start date */
      It_start[j0] = t0_start;
      It_ex[j0] = t0_ex;
    }

    /* exercise dates are now It_ex[j0  ,j0+1  ,...  ,In_ex-1]
       start dates are now It_start[j0  ,j0+1  ,...  ,In_ex-1]

       compute premiums and reduction in first fixed coupon payment */

    i = 0;
    for (j = j0; j < In_ex; j++) {
      while (i < nfix && tfix_pay[i] <= It_start[j])
        i++;
      /* "as of" date for exer. j occurs in period i */
      Iprem[j] = Inprem[i]; /* floating branch + fee                       */
      if (tfix_start[i] < It_start[j]) {
        accrual = full_pay_fix[i] *
                  coverage(tfix_start[i], It_start[j], srt_fix_basis) /
                  coverage(tfix_start[i], tfix_end[i], srt_fix_basis);
        accrual = ((accrual < full_pay_fix[i]) ? accrual : full_pay_fix[i]);
      } /* PH 2/28/98 potential bug fix */
      else
        accrual = 0.;
      if (early_flag_fix == 0)
        Iredpay[j] =
            accrual; /* reduce first payment by accrued fixed interest */
      else {
        Iredpay[j] = 0.; /* else add accrual to the fee */
        Iprem[j] = Iprem[j] + accrual;
      }
    }
    /* adjust floating leg */
    i = 1;
    for (j = j0; j < In_ex; j++) {
      while (i < (nflt - 1) && t_flt[i] <= It_start[j])
        i++; /* "as of" date occurs in floating period i */
      /* BFB: correct the floating leg to properly account for accruals */
      /* Get floating rate */
      if (t_flt[i - 1] >=
          tnow) /* rate has not been fixed; estimate from df's */
        rflt = (1.0 / swp_f_df(t_flt[i - 1], t_flt[i], yldcrv) - 1.0) /
               coverage(t_flt[i - 1], t_flt[i], srt_float_basis);

      else if (rflt_current > 0) /* rate has been fixed; use input rate */
        rflt = rflt_current;

      else /* rate has been fixed  , but is unknown; est. from df's */
      {
        t0 = tnow + t_flt[i] - t_flt[i - 1];
        rflt = (1.0 / swp_f_df(tnow, t0, yc_name) - 1.0) /
               coverage(tnow, t0, srt_float_basis);
      }
      /* Adjust floating leg */
      t0 = ((It_start[j] < t_flt[i - 1]) ? t_flt[i - 1] : It_start[j]);
      if (early_flag_flt == 0)
        Iprem[j] = Iprem[j] - 1. +
                   (1. + rflt * coverage(t0, t_flt[i], srt_float_basis)) *
                       swp_f_df(It_start[j], t_flt[i], yc_name);
      else
        /* Premium = (Premium minus the complete floating leg) - accrued
        interest + *
        /* (aggregation of following floating legs + first floating leg) */
        /* We do premium minus the complete floating leg because we insert a
         * more */
        /* accurate floating leg */
        Iprem[j] =
            (Iprem[j] - 1.) -
            rflt * coverage(t_flt[i - 1], t0, srt_float_basis) +
            (1. + rflt * coverage(t_flt[i - 1], t_flt[i], srt_float_basis)) *
                swp_f_df(It_start[j], t_flt[i], yc_name);
    }

    /* feed autocal */
    In_ex = In_ex - j0; /* effective number of exercise dates */
    i = 0;
    while (i < nfix &&
           tfix_pay[i] <= It_start[j0]) /* PH 2/18/98 potential bug fix */
      i++;           /* i is first coupon after first "as of" date */
    npay = nfix - i; /* effective number of fixed coupons */

    /* evaluate mid-atlantic */
    update_flag = 0;
    if ((update_ts_flag != 0) && (count == 1))
      update_flag = 1;
    find_exer_flag = 0;
    if ((find_exer_bdry != 0) && (count == 1))
      find_exer_flag = 1;
    full_pay_fix[nfix - 1] = full_pay_fix[nfix - 1] + 1.;

    error = LGMautocal(In_ex, &It_ex[j0], &It_start[j0], &Iprem[j0], npay,
                       &tfix_start[i], &tfix_end[i], &tfix_pay[i],
                       &full_pay_fix[i], &Iredpay[j0], pay_rec_str, yc_name,
                       end_of_day_flag, fix_tau_flag, usecaps, In_tau, GetVol,
                       char_vol_type, update_flag, find_exer_flag, outfile,
                       &LGMvalue[count - 1], lgmExerBdryData,
                       /*
                       Not interested in LGM reference swaption data  ,
                       tau/sigma data or intrinsic value
                       */
                       NULL, NULL, NULL, 0);

    full_pay_fix[nfix - 1] = full_pay_fix[nfix - 1] - 1.;

    if (error != NULL)
      return (error);
  } /* end count loop */

  /* extrapolate to get value of the american */
  if (countmax == 1)
    *LGMamervalue = LGMvalue[0]; /* use monthly exercise for american value */
  else                           /* extrapolate (error goes like h*h        */
    *LGMamervalue = (delta_tex[1] * delta_tex[1] * LGMvalue[0] -
                     delta_tex[0] * delta_tex[0] * LGMvalue[1]) /
                    (delta_tex[1] * delta_tex[1] - delta_tex[0] * delta_tex[0]);
  return (error);
}

/* ========================================================================= */
/* ========================================================================= */

/* MAIN PROGRAM to find zeta(t) and G(t) to match market and evaluate
   mid_atlantic */
/* ========================================================================= */

String LGMautocal(
    long In_ex,      /* J = In_ex-1; In_ex is # of exercises      */
    Date It_ex[],    /* notification (exercise) dates   T[0  ,1  ,...  ,J]    */
    Date It_start[], /* start dates for exercise Tj                     s[0  ,1
                        ,...  ,J] */
    double Iprem[],  /* premium (inc. notional) paid to exercise at Tj  f[0  ,1
                        ,...  ,J] */
    long npay,       /* I = npay-1; npay is # of fixed leg coupons       */
    Date tfix_start[], /* start dates for coupon i                   tst[0  ,1
                          , ...  , I]  */
    Date tfix_end[], /* end dates for coupon i                     tend[0  ,1  ,
                        ...  , I] */
    Date tpay[], /* pay dates for coupon i                     tpay[0  ,1  , ...
                    , I] */
    double fixpymnt[], /* coupons (incl. notional in last)           fixpymnt[0
                          ,...  ,I] */
    double Iredpay[],  /* reduction in 1rst payment after exercise   Iredpay[0
                          ,..J]    */
    char *pay_rec_str, /* SRT_RECEIVER or SRT_PAYER */
    String yc_name,    /* pointer to market structures    */
    int end_of_day_flag, /* 0: intra day 1: end of day */
    int fix_tau_flag,    /* 1 means calibrate with fixed tau */
    int usecaps,   /* if fixed tau  , then 1=use caplets to calibrate  , 0=use
                      long swptns */
    double In_tau, /* if fixed tau  , use this value for tau (in years) */
    String (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *), /* function to get swaption vols */
    char *char_vol_type,
    int update_ts_flag, /* 1=compte new sigs and taus & store in yldcrv; 0=don't
                           bother */
    int find_exer_bdry, /* 1 = find swap rates at exercise bocrvary; 0 = don't
                           bother  */
    String outfile,     /* output file name for log file (unused)     */
    double *LGMvalue,   /* LGM value of mid-atlantic             */
    SrtLgmExerBdryData *lgmExerBdryData, /* ptr to exercise boundary data
                                            structure (NULL => not req'd) */
    SrtLgmRefSwptnData *lgmRefSwptnData, /* ptr to reference swaption data
                                            structure (NULL => not req'd) */
    SrtLgmTSData *lgmTSData,  /* ptr to tau/sigma data (NULL => not req'd) */
    double *intrinsicValue,   /* ptr to intrinsic value (NULL => not req'd) */
    int skip_deal_evaluation) /* 0 for a noraml valuation  , 1 for a simple
                                 volatility points extraction */
{
  /* Declarations */
  Date *t_ex = NULL, *t_start = NULL;
  double *pvpay = NULL, *phi_cr = NULL, *redpay = NULL;
  long *ifirst = NULL;
  long i, j, n, nlast, nex;
  long num_months; /* delta_t for term structure. Depends on calibration
                      instruments:
                                   /* 6 for semiannual swaptions  , 12 for
                      annual swaptions */
  /* 3 or 6 if usecaps==1 depending on most liquid caplets; */
  Date tnow, tfirst, nextdate;
  int icmpnd;       /* 0 = annual ref  , swaptions  ,
                      1 = semiannual ref. swaptions */
  double smoothing; /* strike smoothing parameter; 0.0=no smoothing; 1.0=full
                       smoothing */
  BasisCode natural_basis;       /* parameters defining the  */
  SrtCompounding natural_compnd; /* most liquid swaptions in */
  BusDayConv natural_conv;       /* the currency             */
  int natural_spot_lag;          /* S.R. may-07-96 */
  SrtReceiverType srt_p_r_flag;
  SrtDiffusionType srt_vol_type;
  SrtErr srt_err;
  String error_str;
  int end_of_day;
  SrtLgmTS *TS = NULL;
  SrtCurvePtr yldcrv = NULL;
  double intrinsic = 0;
  Date spotdate_temp;

  /* Transform char_vol_type into SrtDiffusionType */

  interp_diffusion_type(char_vol_type, &srt_vol_type);

  /* Pay-Receive flag */
  srt_err = interp_rec_pay(pay_rec_str, &srt_p_r_flag);

  if (srt_err)
    return (srt_err);

  /* initializations */

  s_ycname = yc_name;
  yldcrv = lookup_curve(yc_name);
  error = NULL;
  sq2pi = sqrt(2. * SRT_PI);
  tnow = get_clcndate_from_yldcrv(yldcrv);
  end_of_day = (int)(end_of_day_flag ? 1 : 0);

  pr = 1; /* payer */
  if (srt_p_r_flag == SRT_RECEIVER)
    pr = 0; /* reciever */

  /* define most liquid swaption or caplet market for currency */
  /* sets natural_compnd to SRT_SEMIANNUAL or SRT_ANNUAL  , sets natural_basis;
     sets natural_conv  , and sets natural_spot_lag  */

  error_str =
      get_ccy_defaults(        /* S.R. may-07-96 */
                       yldcrv, /* Market parameters */
                       fix_tau_flag, usecaps, &natural_basis, &natural_compnd,
                       &natural_conv, &natural_spot_lag);

  if (error_str)
    return (error_str);
  else
    ;

  if ((usecaps != 1) || (fix_tau_flag != 1)) /* PH 2/18/98 ModC */
  {
    icmpnd = 0;      /*PH 2/18/98 */
    num_months = 12; /* parameter values when calibrating to annual swaptions */

    if (natural_compnd == SRT_SEMIANNUAL) {
      icmpnd = 1; /*PH 2/18/98 */
      num_months = 6;
    } /* parameter values when calibrating to semiannual swaptions */
    else
      ;
  } else /* PH 2/18/98 ModC */
  {
    icmpnd = 0; /*PH 2/18/98 */
    num_months =
        6; /* parameter values when calibrating to semiannual caplets */

    if (natural_compnd == SRT_QUARTERLY)
      num_months =
          3; /* parameter values when calibrating to quarterly caplets */
    else
      ;
  }
  /* PH 2/18/98 eliminated unused code */
  /* set default tau */
  overall_tau = 50.; /* default tau */
  if (fix_tau_flag == 1)
    overall_tau = In_tau;

  smoothing = 1.0; /* strike smoothing parameter: 0.0 = no smoothing  , 1.0 =
                      full smoothing */

  /* -------------------------*/
  /* STEP ONE: Check inputs   */
  /* -------------------------*/
  /* find the effective number of exercise dates */
  tfirst = tnow;

  /* Check end of day: if true  , move to next business day */
  if (end_of_day == 1)
    tfirst =
        add_unit(tnow, 1, SRT_DAY, SUCCEEDING); /* PH 2/28/98 trivial bug fix */

  nex = -1;
  nextdate = tfirst;
  for (j = 0; j < In_ex; j++) {
    if (It_ex[j] >= nextdate) {
      nex++;
      nextdate = It_ex[j] + MIN_EX_INT;
    }
  }

  if (nex > (MAX_NUM_DATES - 1) || npay > MAX_NUM_DATES) {
    error = "too many dates";
    return (error);
  }

  if (nex < 0 || npay < 1) {
    error = "# of dates <= 0";
    return (error);
  }

  /* Checks the dates for the fixed payments */
  /* Increment i until finds dates in the wrong order: error returned below */
  for (i = 0;
       i < (npay - 1) && tpay[i] <= tpay[i + 1] &&
       tfix_start[i] < tfix_start[i + 1] && tfix_end[i] < tfix_end[i + 1] &&
       tfix_start[i] < tfix_end[i] && tfix_start[i] < tpay[i];
       i++)
    ;

  if (tfix_start[npay - 1] >= tfix_end[npay - 1] ||
      tfix_start[npay - 1] >= tpay[npay - 1])
    i = npay - 2; /* PH 2/28/98 trivial bug fix */

  if (i < (npay - 1)) {
    return serror("Fixed payment dates are incorrect");
  }

  /* Checks exercise dates  */
  for (j = 0; j < (In_ex - 1) && It_ex[j] < It_ex[j + 1] &&
              It_start[j] >= It_ex[j] && It_start[j] <= It_start[j + 1];
       j++)
    ;

  if (It_start[In_ex - 1] < It_ex[In_ex - 1])
    j = In_ex - 2; /* PH 2/28/98 trivial bug fix */

  if (j < (In_ex - 1)) {
    return serror("Exercise dates are incorrect");
  }

  /* Removes all fixed coupons (+ associated dates) before
   * tfirst(today+eofday)*/
  while (npay > 0 && tpay[0] < tfirst) {
    tpay++;       /* Next payment date */
    tfix_start++; /* Next coupon start date */
    tfix_end++;   /* Next coupon end date */
    fixpymnt++;   /* Next coupon */
    npay--;       /* One less coupon to consider */
  }
  if (npay == 0) {
    return serror("The deal has expired: no more fixed coupon!");
  }

  /* -------------------------*/
  /* STEP TWO: Allocate space */
  /* -------------------------*/
  t_ex = NULL;
  t_start = NULL;
  ifirst = NULL;
  pvpay = NULL;
  redpay = NULL;
  phi_cr = NULL;

  t_ex =
      (Date *)srt_calloc(nex + 1, sizeof(Date)); /* t_ex[0  ,1  ,...  ,nex]   */
  t_start =
      (Date *)srt_calloc(nex + 1, sizeof(Date)); /* t_start[0  ,..  ,nex]   */
  ifirst =
      (long *)srt_calloc(nex + 1, sizeof(long)); /* ifirst[0  ,..  ,nex]    */
  pvpay =
      (double *)srt_calloc(npay, sizeof(double)); /* pvpay[0  ,...  ,npay-1] */
  redpay = (double *)srt_calloc(nex + 1,
                                sizeof(double)); /* redpay[0  ,...  ,nex]   */
  phi_cr = (double *)srt_calloc(nex + 1,
                                sizeof(double)); /* phi_cr[0  ,...  ,nex]   */

  if (t_ex == NULL || t_start == NULL || ifirst == NULL || pvpay == NULL ||
      redpay == NULL || phi_cr == NULL) {
    error = "alloc. failed in LGMautocal";
    return (error);
  }

  /* ------------------------------- */
  /* STEP THREE: Copy exercise dates */
  /* ------------------------------- */
  /* ... omit exercise dates closer than MIN_EX_INT days apart  , and omit
     exercise dates before today ... */
  nex = -1;
  nextdate = tfirst;
  for (j = 0; j < In_ex; j++) {
    if (It_ex[j] >= nextdate) {
      nex++;
      t_ex[nex] = It_ex[j];       /* t_ex[0  ,1  ,2  , ...  , nex]    */
      t_start[nex] = It_start[j]; /* t_start[0  ,1  ,2  , ...  , nex] */
      origprem[nex] = Iprem[j];   /* prem[0  ,1  ,2  , ...  , nex]    */
      redpay[nex] = Iredpay[j];   /* redpay[0  ,1  , ...  , nex]    */
      nextdate = t_ex[nex] + MIN_EX_INT;
    }
  }

  /* eliminate any exercise dates that have no underlying left */

  for (; nex >= 0 && t_start[nex] > (tpay[npay - 1] - 10); nex--)
    ; /* PH 2/28/98 trivial bug fix */
  if (nex < 0) {
    error = "no valid exercise dates";
    ERROR_CHECK;
  }

  /* This is the deal that will be evaluated:
          exercise dates t_ex[j]  , j=0  ,1  ,...  ,nex
          if exercised at t_ex[j]  , then fixed leg will receive the payments
                  fixpymnt[i]-redpay[i] paid on tpay[i] for i=ifirst[j]
                  fixpymnt[i] paid on tpay[i] for i=ifirst[j]+1  , ifirst[j]+2
     , ...  , npay-1 the floating leg will receive a single payment of
     origprem[j] paid on  t_start[j]
          ifirst[j] (calculated next) is first i for which tpay[i]>t_start[j] */

  /* Miscellaneous bookkeeping:
    find first coupon paid for each exercise date ... ifirst[j]
    compute pv of the interest payments at each pay date ... pvpay[i]
    compute pv of the premium payments at each start date ... prem[j]
    compute pv of the 1rst payment reduction at each start ...redpay[j] */

  error = pvcashflow(tnow, yc_name, nex, t_ex, t_start, npay, tfix_start,
                     tfix_end, tpay, fixpymnt, ifirst, pvpay, redpay);
  ERROR_CHECK;

  /* compute the intrinsic value of the mid_atlantic (in case its
     wanted at some future date) */
  intrinsic = comp_intrin_val(nex, ifirst, redpay, npay, pvpay);
  if (intrinsicValue)
    *intrinsicValue = intrinsic;

  /* Check the price if we are at the last exercice date */
  if ((nex == 0) && (tfirst == t_ex[0])) {
    /* we have a swap */
    *LGMvalue = intrinsic;
    FreeMost(t_ex, t_start, ifirst, pvpay, phi_cr, redpay);
    return NULL;
  }

  /* ------------------------------------- */
  /* STEP FOUR: Choose reference swaptions */ /* This next part has been
                                                 simplified */
  /* ------------------------------------- */
  /* choose dates s[0  , ...  , n] for reference swaptions */

  choose_dates(tnow, t_start[0], &n, npay, tpay, nex, t_ex, icmpnd, num_months,
               natural_conv);
  nlast = ((icmpnd == 1) ? n - 2 : n - 1);

  for (j = 0; j <= n; j++) {
    sy[j] = coverage(tnow, s[j],
                     BASIS_ACT_365); /* sy[0  ,..  ,n] years from tnow    */
    B[j] = swp_f_df(tnow, s[j], yc_name);
  } /* B[0  ,..  ,n] today's dis.factors */
  aB[0] = 0.;
  for (j = 1; j <= n; j++)
    aB[j] = B[j] *
            coverage(s[j - 1], s[j], natural_basis); /* discounted coverages */

  /* We can have some problems if there is only one exercise */
  /* and the exercise notification date is before spot lag */
  /* Solution we put the spot of the yield curve equal to t_ex[0] */
  if ((nex == 0) && (t_ex[0] < get_spotdate_from_yldcrv(yldcrv))) {
    spotdate_temp = get_spotdate_from_yldcrv(yldcrv);
    set_spotdate_from_yldcrv(yldcrv, t_ex[0]);
  }

  /* Reference instruments:
          If fix_tau_flag != 1  , code uses swaptions with start dates are s[j]
     , j=0  ,1  ,...  ,nlast. End dates for short swaptions are s[j+1] (annual
     markets) or s[j+2] (semiannual markets). End dates for long swaptions are
     s[n]. Pay dates are s[i]  , i=j+1  ,...  ,n. Exercise dates precede s[j] by
     natural_spot_lag. If fix_tau_flag == 1  , and usecaps !=1 code uses long
     swaptions but not the short ones If fix_tau_flag == 1 and usecaps ==1  ,
     code uses caps with start dates s[j] and end dates
                          s[j+1]  , j=0  ,1  ,...  ,nlast*/

  /* Choose fixed rates (strikes) for reference swaptions */
  /* Find the minimum & maximum strikes allowed
        (only allow R to be FAR_OUT Black-Scholes std deviations from par) */

  error = R_restrictions(GetVol, srt_vol_type, n, nlast, icmpnd, natural_compnd,
                         natural_basis, fix_tau_flag, usecaps);
  /*sets Rmin_sh[j]  , Rmax_sh[j]  , Rmin_lg[j]  , and Rmax_lg[j] */
  ERROR_CHECK;

  if ((fix_tau_flag != 1) || (usecaps != 1))
    error =
        Olivier_lg(n, nlast, icmpnd, nex, npay, ifirst, t_start, pvpay, redpay);
  ERROR_CHECK;

  if ((fix_tau_flag != 1) || (usecaps == 1))
    error = Olivier_sh(n, nlast, nex, npay, ifirst, t_start, pvpay, redpay,
                       icmpnd, fix_tau_flag);
  ERROR_CHECK;

  /* To simplify the coding  , when we don't need both long and short rates  ,
     we copy the needed one into the unneeded one */
  if ((fix_tau_flag == 1) && (usecaps != 1)) {
    for (j = 0; j <= nlast; j++) {
      R_theor_sh[j] = R_theor_lg[j];
      Rmin_sh[j] = Rmin_lg[j];
      Rmax_sh[j] = Rmax_lg[j];
    }
  }
  if ((fix_tau_flag == 1) && (usecaps == 1)) {
    for (j = 0; j <= nlast; j++) {
      R_theor_lg[j] = R_theor_sh[j];
      Rmin_lg[j] = Rmin_sh[j];
      Rmax_lg[j] = Rmax_sh[j];
    }
  }

  /* Smooth strikes & restrict strikes to permitted range */
  smooth_strikes(nlast, smoothing);

  R_theor_sh[nlast] = R_theor_lg[nlast];

  for (j = 0; j <= nlast; j++) {
    Ract_sh[j] = (R_theor_sh[j] > Rmin_sh[j] ? R_theor_sh[j] : Rmin_sh[j]);
    Ract_sh[j] = (Ract_sh[j] < Rmax_sh[j] ? Ract_sh[j] : Rmax_sh[j]);
    Ract_lg[j] = (R_theor_lg[j] > Rmin_lg[j] ? R_theor_lg[j] : Rmin_lg[j]);
    Ract_lg[j] = (Ract_lg[j] < Rmax_lg[j] ? Ract_lg[j] : Rmax_lg[j]);
  }

  /* --------- */
  /* STEP FIVE */
  /* --------- */
  /*   Fit zeta[t] and g[t] so that the LGM value of the reference
     swaptions or caplets match market value */
  if (fix_tau_flag != 1)
    error = mid_atl_fit_market(GetVol, srt_vol_type, n, nlast, icmpnd,
                               natural_compnd, natural_basis, natural_spot_lag,
                               lgmRefSwptnData);
  else
    error = mid_atl_fit_market_fix_tau(GetVol, srt_vol_type, n, nlast, icmpnd,
                                       natural_compnd, natural_basis, usecaps,
                                       natural_spot_lag, lgmRefSwptnData);
  ERROR_CHECK;

  /* ------------ */
  /* STEP SIX: Find the value of the midatlantic using the newly fitted
     G(t) and zeta(t) */
  /* ------------ */

  if (!skip_deal_evaluation) {
    error = mid_atl_eval(tnow, n, yldcrv,                   /* misc. info   */
                         nex, t_ex, t_start,                /* option info  */
                         npay, tpay, pvpay, redpay, ifirst, /* payment info */
                         LGMvalue, phi_cr);                 /* output       */
    ERROR_CHECK;
  }
  /* ------------ */
  /* STEP SEVEN: Find critical strikes (strikes at the optimal exercise point
                 in case they are wanted in the future                        */
  /* ------------ */

  if ((!skip_deal_evaluation) && (find_exer_bdry != 0)) {
    error = Findexerbdry(tnow, nex, t_ex, tpay[npay - 1], phi_cr, n, yc_name,
                         icmpnd, natural_basis, natural_conv, natural_spot_lag,
                         lgmExerBdryData);
    if (error != NULL) {
      error = "exer. bdry not found";
      ERROR_CHECK;
    }
  }

  /* We had some problems if there is only one exercise */
  /* and the exercise notification date is before spot lag */
  /* Solution we put back the spot of the yield curve */
  if ((nex == 0) && (t_ex[0] == get_spotdate_from_yldcrv(yldcrv))) {
    set_spotdate_from_yldcrv(yldcrv, spotdate_temp);
  }

  /* ------------ */
  /* STEP EIGHT: Find sigma and tau values equivalent to g[j] and zeta[j] */
  /* ------------ */

  if ((update_ts_flag != 0) && (lgmTSData)) {
    error = Findsigtau(n, tnow, natural_spot_lag);
    if (error != 0) {
      error = "term structure not reset";
      ERROR_CHECK;
    }

    lgmTSData->NTS = num_sig; /* = num_tau = common no. of tau/sig points */

    if (!lgmTSData->TSArr)
      lgmTSData->TSArr = srt_calloc(num_sig, sizeof(SrtLgmTS));

    for (i = 0; i < num_sig; ++i) {
      TS = lgmTSData->TSArr + i;

      TS->tauDt = (long)tau_date[i];
      TS->tau = tau_val[i];
      TS->sigmaDt = (long)sig_date[i];
      TS->sigma = sig_val[i];
    }
  }

  /* Free everything and return  */
  FreeMost(t_ex, t_start, ifirst, pvpay, phi_cr, redpay);
  return (error);
}

/* ======================================================================== */
/* ======================================================================== */
/* ======================================================================== */

/*******************************************************************************
Static functions
*******************************************************************************/

/* fit g[j] and zeta[j] by matching LGM values of vanilla swaptions to
   market (BS) values. Swaptions used:
      start s[j]  , end s[n]  , strike Ract_lg[j]  ,  (fit exactly)
      start s[j]  , end one year later  , strike Ract_sh[j] (fit whenever
   possible) for j = 0  ,1  , ...  ,nlast
   Uses global Newton routine (semiannual or annual) for fit                */

static String mid_atl_fit_market(String (*GetVol)(),
                                 SrtDiffusionType srt_vol_type, long n,
                                 long nlast, int icmpnd,
                                 SrtCompounding compounding, SrtBasisCode basis,
                                 int spot_lag,
                                 SrtLgmRefSwptnData *lgmRefSwptnData) {
  double sqznew, sqzmin, sqzmax;
  double gnew, gmin, gmax;
  double phistsh, phistlg;
  double v1, sig;
  long j, dj;
  int check;
  Err err;
  SwapDP sdp;
  SrtLgmRefSwptn *refSwptn = NULL;
  double fwd;

  error = NULL;

  if (lgmRefSwptnData) {
    lgmRefSwptnData->NrefShortSwptn = nlast + 1;
    lgmRefSwptnData->NrefLongSwptn = nlast + 1;
    lgmRefSwptnData->isSemi = icmpnd;

    if (!lgmRefSwptnData->refShortSwptnArr)
      lgmRefSwptnData->refShortSwptnArr =
          srt_calloc(nlast + 1, sizeof(SrtLgmRefSwptn));
    if (!lgmRefSwptnData->refLongSwptnArr)
      lgmRefSwptnData->refLongSwptnArr =
          srt_calloc(nlast + 1, sizeof(SrtLgmRefSwptn));
  }

  /* Get Black-Scholes values */
  dj = 2;
  if (icmpnd == 0)
    dj = 1;

  for (j = 0; j <= nlast; j++) {
    /* SHORT Swaptions */
    if (err = swp_f_setSwapDP(s[j], s[j + dj], compounding, basis, &sdp))
      return err;

    err = swp_f_ForwardRate_SwapDP(&sdp, s_ycname, "CASH", &fwd);
    if (err)
      return err;

    /* The Vol is used ATM
            error = (*GetVol)(s[j]  , s[j+dj]  , fwd  , &sig);
 */
    error = (*GetVol)(s[j], s[j + dj], Ract_sh[j], 0, &sig);
    if (error != NULL) {
      error = "failed to get vol in fit market";
      return (error);
    }
    /*BS value(short swaptns)*/
    jst = j;
    nst = j + dj;
    v_bs_sh[j] = Vbs_jn(Ract_sh[j], sig, srt_vol_type, spot_lag);

    /* Save short term swaption details */
    if (lgmRefSwptnData) {
      refSwptn = lgmRefSwptnData->refShortSwptnArr + j;

      refSwptn->bgnDt = s[j];
      refSwptn->endDt = s[j + dj];
      refSwptn->strike = Ract_sh[j];
      refSwptn->bsVol = sig;
      refSwptn->bsPv = v_bs_sh[j];
    }

    /* LONG Swaptions */
    if (err = swp_f_setSwapDP(s[j], s[n], compounding, basis, &sdp))
      return err;

    err = swp_f_ForwardRate_SwapDP(&sdp, s_ycname, "CASH", &fwd);
    if (err)
      return err;

    /* The Vol is used ATM
     error = (*GetVol)(s[j]  , s[n]  , fwd  , &sig);
*/
    error = (*GetVol)(s[j], s[n], Ract_lg[j], 0, &sig);
    if (error != NULL) {
      error = "failed to get vol in fit market";
      return (error);
    }
    /*BS value(long swaptns)*/
    jst = j;
    nst = n;
    v_bs_lg[j] = Vbs_jn(Ract_lg[j], sig, srt_vol_type, spot_lag);

    /* Save long term swaption details */
    if (lgmRefSwptnData) {
      refSwptn = lgmRefSwptnData->refLongSwptnArr + j;

      refSwptn->bgnDt = s[j];
      refSwptn->endDt = s[n];
      refSwptn->strike = Ract_lg[j];
      refSwptn->bsVol = sig;
      refSwptn->bsPv = v_bs_lg[j];
    }
  }

  /* Set g[nlast]  , ...  , g[n] to default values  */
  phistsh = 0.;
  phistlg = 0.;
  g[n] = 0.;
  g[n - 1] = overall_tau * (exp(-sy[n - 1] / overall_tau) -
                            exp(-sy[n] / overall_tau)); /* wolog ! */
  if (icmpnd != 0)
    g[n - 2] = overall_tau *
               (exp(-sy[n - 2] / overall_tau) - exp(-sy[n] / overall_tau));

  /* Find zeta[nlast] by matching market price */
  sqzmax = 2. * sig * sqrt(sy[nlast]) /
           g[nlast]; /* fit LGM value computed by Vjn  */
  sqznew =
      sqzmax * (1. - B[n] / B[nlast]) / 2.; /* to BS value in v_bs_sh[nlast]  */
  sqzmin = .2 * sqznew;
  Ract = Ract_sh[nlast];
  jst = nlast;
  nst = n;
  phiact = phistsh;
  /* Solves for the last swaption : only zeta to be fit */
  ScSlowNewt(&sqznew, sqzmin, sqzmax, /* fit zeta */
             &check, v_bs_sh[nlast], 0, Vjn, Vjn_der);
  if (check != 0) {
    error = "Newton failed at zeta[nlast]";
    return (error);
  }
  phistsh = phiact;
  phistlg = phiact;

  /*  update zeta array */
  zeta[nlast] = sqznew * sqznew;
  zeta[n - 1] = zeta[nlast] * sy[n - 1] / sy[nlast];
  zeta[n] = zeta[nlast] * sy[n] / sy[nlast];

  /* Step backwards in j  , finding g(t) and zeta(t) at t=s[j] by
     fitting short and long swaptions with exercise date sj to market */
  for (j = (nlast - 1); j >= 0; j--) {
    /* set maximum zeta[j] and minimum g[j] values */
    g[j] = g[j + 1] +
           0.05 * (g[j + 1] - g[n]) * (sy[j + 1] - sy[j]) / (sy[n] - sy[j + 1]);
    zeta[j] = zeta[j + 1];

    /* find LGM value of long swaption */
    jst = j;
    nst = n;
    Ract = Ract_lg[j];
    phiact = phistlg;
    v1 = Vjn(g[j], 1); /* Values the swaption */
    phistlg = phiact;

    /* fit long swaption by changing either zeta or g  , keeping the other at
     * its extreme value */
    if (v1 > v_bs_lg[j]) /* decrease zeta[j] to fit long swaption to market */
    {
      sqzmax = sqrt(zeta[j]);
      sqznew = sqzmax * sqrt(sy[j] / sy[j + 1]);
      if (sqznew > sqzmax)
        sqznew = sqzmax;
      sqzmin = .5 * sqznew;
      ScSlowNewt(&sqznew, sqzmin, sqzmax, /* fit long swaption */
                 &check, v_bs_lg[j], 0, Vjn, Vjn_der);
      if (check != 0) {
        error = "Newton failed at zeta[j]";
        return (error);
      } else
        zeta[j] = sqznew * sqznew;
    }

    else if (v1 < v_bs_lg[j]) /* increase g[j] to fit long swaption to market */
    {
      gmin = g[j];
      gnew = g[j + 1] +
             (g[j + 1] - g[n]) * (sy[j + 1] - sy[j]) / (sy[n] - sy[j + 1]);
      gmax = 2.0 * gnew;
      ScSlowNewt(&gnew, gmin, gmax, /* fit long swaption */
                 &check, v_bs_lg[j], 1, Vjn, Vjn_der);
      if (check != 0) {
        error = "Newton failed at zeta[j]";
        return (error);
      } else
        g[j] = gnew;
    } else
      ;
    phistlg = phiact;

    /* if the LGM value of the short swaption is less  than the BS value  , then
       a unique pair zeta[j] and  g[j] exist which fit LGM values exactly ...
       else we cannot do better than the bdry values we already have */

    /* find LGM value of short swaption */
    jst = j;
    nst = j + dj;
    Ract = Ract_sh[j];
    phiact = phistsh;
    v1 = Vjn(g[j], 1);
    phistsh = phiact;
    if (v1 >= v_bs_sh[j])
      continue; /* can't improve; go to new j */

    /* can improve: use Newton to fit both g[j] and zeta[j]  */
    if (icmpnd != 0)
      error = Newt_Semi(j, n, &phistlg, &phistsh);
    else
      error = Newt_Annual(j, n, &phistlg, &phistsh);
    if (error != NULL)
      return (error); /* didn't converge ... return bad news */
    else
      ; /* converged; get new j */
  }     /* end j loop */

  return (error);
}

/* fit zeta[j] by matching LGM values of vanilla long swaptions (if usecaps!=1)
        or caplets (if usecaps==1) to market (BS) values.
        Reference swaptions:
      start s[j]  , end s[n]  , strike Ract_lg[j]	 for j = 0  ,1  , ...
   ,nlast Reference caplets start s[j]  , end s[j+1]  , strike Ract_sh[j]  for j
   = 0  ,1  , ...  ,nlast
   Uses global Newton routine (semiannual or annual) for fit                */

static String mid_atl_fit_market_fix_tau(
    String (*GetVol)(), SrtDiffusionType srt_vol_type, long n, long nlast,
    int icmpnd, SrtCompounding compounding, SrtBasisCode basis, int usecaps,
    int spot_lag, SrtLgmRefSwptnData *lgmRefSwptnData) {
  double sqznew, sqzmin, sqzmax, target;
  double v1, sig;
  long j, dj;
  int check;
  Err err;
  SwapDP sdp;
  SrtLgmRefSwptn *refSwptn = NULL;
  double fwd;

  error = NULL;
  dj = 1;

  if (lgmRefSwptnData) {
    lgmRefSwptnData->NrefShortSwptn = nlast + 1;
    lgmRefSwptnData->NrefLongSwptn = nlast + 1;
    lgmRefSwptnData->isSemi = icmpnd;

    if (!lgmRefSwptnData->refShortSwptnArr)
      lgmRefSwptnData->refShortSwptnArr =
          srt_calloc(nlast + 1, sizeof(SrtLgmRefSwptn));
    if (!lgmRefSwptnData->refLongSwptnArr)
      lgmRefSwptnData->refLongSwptnArr =
          srt_calloc(nlast + 1, sizeof(SrtLgmRefSwptn));
  }

  /* Get Black-Scholes values */
  for (j = 0; j <= nlast; j++) {
    if (usecaps == 1) /* Caps */
    {
      if (err = swp_f_setSwapDP(s[j], s[j + dj], compounding, basis, &sdp))
        return err;

      err = swp_f_ForwardRate_SwapDP(&sdp, s_ycname, "CASH", &fwd);
      if (err)
        return err;

      /* The Vol is used ATM
      error = (*GetVol)(s[j]  , s[j+dj]  , fwd  , &sig);
      */
      error = (*GetVol)(s[j], s[j + dj], Ract_sh[j], 1, &sig);
      if (error != NULL) {
        error = "failed to get vol in fit market";
        return (error);
      }
      /*BS value(caps)*/
      jst = j;
      nst = j + dj;
      v_bs_sh[j] = Vbs_jn(Ract_sh[j], sig, srt_vol_type, spot_lag);

      /* Save short term swaption details */
      if (lgmRefSwptnData) {
        refSwptn = lgmRefSwptnData->refShortSwptnArr + j;

        refSwptn->bgnDt = s[j];
        refSwptn->endDt = s[j + dj];
        refSwptn->strike = Ract_sh[j];
        refSwptn->bsVol = sig;
        refSwptn->bsPv = v_bs_sh[j];
      }
    }

    if (usecaps != 1) /* LONG Swaptions */
    {
      if (err = swp_f_setSwapDP(s[j], s[n], compounding, basis, &sdp))
        return err;

      err = swp_f_ForwardRate_SwapDP(&sdp, s_ycname, "CASH", &fwd);
      if (err)
        return err;

      /* The Vol is used ATM
           error = (*GetVol)(s[j]  , s[n]  , fwd  , &sig);
  */
      error = (*GetVol)(s[j], s[n], Ract_lg[j], 0, &sig);
      if (error != NULL) {
        error = "failed to get vol in fit market";
        return (error);
      }

      /*BS value(long swaptns)*/
      jst = j;
      nst = n;
      v_bs_lg[j] = Vbs_jn(Ract_lg[j], sig, srt_vol_type, spot_lag);

      /* Save long term swaption details */
      if (lgmRefSwptnData) {
        refSwptn = lgmRefSwptnData->refLongSwptnArr + j;

        refSwptn->bgnDt = s[j];
        refSwptn->endDt = s[n];
        refSwptn->strike = Ract_lg[j];
        refSwptn->bsVol = sig;
        refSwptn->bsPv = v_bs_lg[j];
      }
    }
  }

  /* Define g[0  ,1  ,...  ,n] from tau  */
  g[n] = 0.;
  for (j = 0; j < n; j++)
    g[j] =
        overall_tau * (exp(-sy[j] / overall_tau) - exp(-sy[n] / overall_tau));

  /* Find zeta[nlast] by matching market price */
  sqzmax = 2. * sig * sqrt(sy[nlast]) /
           g[nlast]; /* fit LGM value computed by Vjn  */
  sqznew = sqzmax * (1. - B[n] / B[nlast]) / 2.; /* to BS value  */
  sqzmin = .2 * sqznew;
  if (usecaps == 1) {
    jst = nlast;
    nst = jst + 1;
    Ract = Ract_sh[jst];
    target = v_bs_sh[nlast];
  } else {
    jst = nlast;
    nst = n;
    Ract = Ract_lg[nlast];
    target = v_bs_lg[nlast];
  }
  phiact = 0.;
  ScSlowNewt(&sqznew, sqzmin,
             sqzmax, /* Find zeta to match last swaption or cap */
             &check, target, 0, Vjn, Vjn_der);
  if (check != 0) {
    error = "Newton failed at zeta[nlast]";
    return (error);
  }

  /*  update zeta array */
  zeta[nlast] = sqznew * sqznew;
  zeta[n - 1] = zeta[nlast] * sy[n - 1] / sy[nlast];
  zeta[n] = zeta[nlast] * sy[n] / sy[nlast];

  /* Step backwards in j  , finding g(t) and zeta(t) at t=s[j] by
     fitting caplet or long swaption with exercise date sj to market */
  for (j = (nlast - 1); j >= 0; j--) {
    zeta[j] = zeta[j + 1]; /* set maximum zeta[j] value */

    if (usecaps == 1) /* find LGM value of swaption or cap */
    {
      jst = j;
      nst = jst + 1;
      Ract = Ract_sh[jst];
      target = v_bs_sh[jst];
    } else {
      jst = j;
      nst = n;
      Ract = Ract_lg[jst];
      target = v_bs_lg[jst];
    }

    v1 = Vjn(g[j], 1); /* Values the swaption */
                       /* fit swaption or cap by changing  zeta */
    if (v1 > target)   /* decrease zeta[j] to fit long swaption or caplet */
    {
      sqzmax = sqrt(zeta[j]);
      sqznew = sqzmax * sqrt(sy[j] / sy[j + 1]);
      if (sqznew > sqzmax)
        sqznew = sqzmax;
      sqzmin = .5 * sqznew;
      ScSlowNewt(&sqznew, sqzmin, sqzmax, /* fit swaption or caplet */
                 &check, target, 0, Vjn, Vjn_der);
      if (check != 0) {
        error = "Newton failed at zeta[j]";
        return (error);
      } else
        zeta[j] = sqznew * sqznew;
    } else
      zeta[j] = zeta[j + 1]; /* CANNOT fit to market; use best value */
  }                          /* end j loop */

  /* the array zeta[j] has actually been fitted at the exercise dates of the
     reference swaptions
     instead of s[j] ... use linear interpolation to correct zeta[j] for
     spot_lag */

  /*    for (j=0; j<=nlast; j++)
      {   sly = s[j] - add_unit(s[j]  , -spot_lag  , SRT_BDAY  , SUCCEEDING);
          sly = sly/365.0;
          zeta[j] = zeta[j] + (zeta[j+1]-zeta[j])*sly/(sy[j+1]-sy[j]);
      }
      slope = (zeta[nlast]-zeta[nlast-1])/(sy[nlast]-sy[nlast-1]);
      if (slope > 1.0e-8)
      {   for (j=(nlast+1); j<=n; j++)
              zeta[j] = zeta[nlast] + slope*(sy[j]-sy[nlast]);
      }
  */
  return (error);
}

/***********************************************************************/
static String Newt_Semi(long j, long n, double *philgptr, double *phishptr) {
  double sqznew, sqzold, sqzmin, sqzmax;
  double gnew, gold, gmin, gmax;
  double zetasave, gsave;
  double stepz, stepg;
  double v1, v2, v11, v12, v21, v22;
  double lambda, missby, crit1, crit2;
  long i, k;

  /* set extreme values & initial guess for zeta[j] */
  sqzmax = zeta[j];
  zeta[j] = zeta[j + 1] * sy[j] / sy[j + 1];
  if (zeta[j] > sqzmax)
    zeta[j] = sqzmax;
  zetasave = zeta[j];
  sqzmax = sqrt(sqzmax);
  sqznew = sqrt(zeta[j]);
  sqzmin = sqznew / 2.0;

  /* set extreme values & initial guess for g[j] */
  gmin = g[j];
  gnew = 2. * g[j + 1] - g[j + 2];
  if (gnew < gmin)
    gnew = gmin;
  g[j] =
      g[j + 1] + (g[j + 1] - g[n]) * (sy[j + 1] - sy[j]) / (sy[n] - sy[j + 1]);
  if (g[j] > gnew)
    gnew = g[j];
  g[j] = gnew;
  gsave = gnew;
  gmax = 2.0 * gnew;

  stepz = .15 * sqzmax; /* set maximum step lengths */
                        /*    stepg = .15*gmin; */
  stepg = .15 * gmax;

  crit1 = 1.0e-8 * (1. + fabs(gnew)); /* set convergence criterion */
  crit2 = 1.0e-8 * sqznew + 1.0e-12;

  /* get error and its derivatives */
  jst = j;
  nst = j + 2;
  Ract = Ract_sh[j];
  phiact = *phishptr;
  Vjn_all(&v1, &v11, &v12);
  v1 = v1 / v_bs_sh[j] - 1.0; /* error in short swaption */
  v11 = v11 / v_bs_sh[j];     /* deriv. wrt sqrt(zeta)   */
  v12 = v12 / v_bs_sh[j];     /* deriv. wrt g[j]         */
  *phishptr = phiact;

  jst = j;
  nst = n;
  Ract = Ract_lg[j];
  phiact = *philgptr;
  Vjn_all(&v2, &v21, &v22);
  v2 = v2 / v_bs_lg[j] - 1.0; /* error in long swaption  */
  v21 = v21 / v_bs_lg[j];     /* deriv. wrt sqrt(zeta)   */
  v22 = v22 / v_bs_lg[j];     /* deriv. wrt g[j]         */
  *philgptr = phiact;

  /* main Newton iteration */
  for (i = 1; i <= 40; i++) {
    sqzold = sqznew; /* store old values */
    gold = gnew;
    missby = v1 * v1 + v2 * v2;

    /* new Newton step */
    if (fabs(v11) < 0.3e-3 && fabs(v12) < 0.3e-3) {
      sqznew =
          sqzold - v2 * v21 * sqzold * sqzold /
                       (sqzold * sqzold * v21 * v21 + gold * gold * v22 * v22);
      gnew = gold - (v2 + v21 * (sqznew - sqzold)) / v22;
    } else {
      sqznew = sqzold - (v22 * v1 - v12 * v2) / (v11 * v22 - v12 * v21);
      gnew = gold - (-v21 * v1 + v11 * v2) / (v11 * v22 - v12 * v21);
    }

    /* reduce step if out of bounds */
    if (fabs(gnew - gold) > stepg) {
      sqznew = sqzold + (sqznew - sqzold) * stepg / fabs(gnew - gold);
      gnew = gold + (gnew - gold) * stepg / fabs(gnew - gold);
    }
    if (fabs(sqznew - sqzold) > stepz) {
      gnew = gold + (gnew - gold) * stepz / fabs(sqznew - sqzold);
      sqznew = sqzold + (sqznew - sqzold) * stepz / fabs(sqznew - sqzold);
    }
    if (gnew <= gmin) {
      sqznew = sqzold + (sqznew - sqzold) * (gmin - gold) / (gnew - gold);
      gnew = gmin;
    } else if (gnew >= gmax) {
      sqznew = sqzold + (sqznew - sqzold) * (gmax - gold) / (gnew - gold);
      gnew = gmax;
    }
    if (sqznew <= sqzmin) {
      gnew = gold + (gnew - gold) * (sqzmin - sqzold) / (sqznew - sqzold);
      sqznew = sqzmin;
    } else if (sqznew >= sqzmax) {
      gnew = gold + (gnew - gold) * (sqzmax - sqzold) / (sqznew - sqzold);
      sqznew = sqzmax;
    }
    /* update arrays */
    g[j] = gnew;
    zeta[j] = sqznew * sqznew;

    /* calculate values at new point */
    jst = j;
    nst = j + 2;
    Ract = Ract_sh[j];
    phiact = *phishptr;
    Vjn_all(&v1, &v11, &v12);
    v1 = v1 / v_bs_sh[j] - 1.0; /* error in short swaption */
    v11 = v11 / v_bs_sh[j];     /* deriv. wrt sqrt(zeta)   */
    v12 = v12 / v_bs_sh[j];     /* deriv. wrt g[j]         */
    *phishptr = phiact;

    jst = j;
    nst = n;
    Ract = Ract_lg[j];
    phiact = *philgptr;
    Vjn_all(&v2, &v21, &v22);
    v2 = v2 / v_bs_lg[j] - 1.0; /* error in long swaption  */
    v21 = v21 / v_bs_lg[j];     /* deriv. wrt sqrt(zeta)   */
    v22 = v22 / v_bs_lg[j];     /* deriv. wrt g[j]         */
    *philgptr = phiact;

    /* if needed  , reduce Newton step length until error decreases */
    for (k = 1; k <= 5 && (v2 * v2 + v1 * v1 >= missby); k++) {
      lambda = missby / (missby + v2 * v2 + v1 * v1);
      if (lambda < .1)
        lambda = .1;
      sqznew = sqzold + lambda * (sqznew - sqzold);
      gnew = gold + lambda * (gnew - gold);
      zeta[j] = sqznew * sqznew;
      g[j] = gnew;

      /* compute values at reduced step */
      jst = j;
      nst = j + 2;
      Ract = Ract_sh[j];
      phiact = *phishptr;
      Vjn_all(&v1, &v11, &v12);
      v1 = v1 / v_bs_sh[j] - 1.0; /* error in short swaption */
      v11 = v11 / v_bs_sh[j];     /* deriv. wrt sqrt(zeta)   */
      v12 = v12 / v_bs_sh[j];     /* deriv. wrt g[j]         */
      *phishptr = phiact;

      jst = j;
      nst = n;
      Ract = Ract_lg[j];
      phiact = *philgptr;
      Vjn_all(&v2, &v21, &v22);
      v2 = v2 / v_bs_lg[j] - 1.0; /* error in long swaption  */
      v21 = v21 / v_bs_lg[j];     /* deriv. wrt sqrt(zeta)   */
      v22 = v22 / v_bs_lg[j];     /* deriv. wrt g[j]         */
      *philgptr = phiact;
    }

    /* check convergence */
    if (fabs(gnew - gold) < crit1 && fabs(sqznew - sqzold) < crit2)
      break;
    if ((v1 * v1 + v2 * v2) < 1.0e-22)
      break;

    /* no convergence yet  , take another Newton step */
  } /* end i loop */

  /*  Check for convergence using loose criterion */
  error = NULL;
  crit1 = 10.0 * crit1;
  crit2 = 10.0 * crit2;
  if (fabs(gnew - gold) < crit1 && fabs(sqznew - sqzold) < crit2)
    return (error);
  else if ((v1 * v1 + v2 * v2) < 1.0e-18)
    return (error);

  else if (j == 0) /* conv. failed  , but not disastrous */
  {
    zeta[j] = zetasave;
    g[j] = gsave;
    return (error);
  }

  else {
    error = "Newton failed in fit market"; /* disaster */
    return (error);
  }
}

/***********************************************************************/
/* Here we exploit the fact that a short swaption is an option on an FRA:
   the link between g and zeta is explicit...
*/
static String Newt_Annual(long j, long n, double *philgptr, double *phishptr) {
  double sqznew, sqzmin, sqzmax;
  double gnew, gmin, gmax;
  int check;

  error = NULL;
  /* set extreme values and initial guess for g[j] */
  gmin = g[j];
  gnew = 2. * g[j + 1] - g[j + 2];
  if (gnew < gmin)
    gnew = gmin;
  g[j] =
      g[j + 1] + (g[j + 1] - g[n]) * (sy[j + 1] - sy[j]) / (sy[n] - sy[j + 1]);
  if (g[j] > gnew)
    gnew = g[j];
  g[j] = gnew;
  gmax = 2.0 * gnew;

  /* fit short swaption */
  jst = j;
  nst = j + 1;
  Ract = Ract_sh[j];
  phiact = *phishptr;
  ScSlowNewt(&gnew, gmin, gmax, &check, v_bs_sh[j], 1, Vjn, Vjn_der);
  if (check != 0) {
    error = "Newton failed at zeta[j]";
    return (error);
  } else
    g[j] = gnew;
  *phishptr = phiact;

  kappa = (g[j] - g[j + 1]) * sqrt(zeta[j]);

  /* Fit large swaption on curve (g[j]-g[j+1])*sqrt(zeta[j]) = kappa   */
  sqzmax = sqrt(zeta[j]);
  sqznew = sqrt(zeta[j] * sy[j] / sy[j + 1]);
  if (sqznew > sqzmax)
    sqznew = sqzmax;
  sqzmin = .5 * sqznew;

  jst = j;
  nst = n;
  Ract = Ract_lg[j];
  phiact = *philgptr;
  ScSlowNewt(&sqznew, sqzmin, sqzmax, &check, v_bs_lg[j], 1, Vjn1, Vjn1_der);
  if (check != 0) {
    error = "Newton failed at zeta[j]";
    return (error);
  }
  *philgptr = phiact;

  /* update array with new values */
  zeta[j] = sqznew * sqznew;
  g[j] = g[j + 1] + kappa / sqznew;
  return (error);
}

/***********************************************************************/
static String pvcashflow(Date tnow, String ycName, long nex, Date t_ex[],
                         Date t_start[], /* exercise info */
                         long npay, Date tfix_start[], Date tfix_end[],
                         Date tpay[], double fixpymnt[], /* payment info  */
                         long ifirst[], double pvpay[], double redpay[])
/* output        */

/* find first cpn whose end date & pay date are after each start
        date ... ifirst[j]
  compute pv of the interest payments at each pay date ... pvpay[i]
  compute pv of the premium payments at each start date ... prem[j]
  compute pv of the reduction in 1rst payment at each start date redpay[j] */
{
  long i, j;
  double df;

  error = "no disc. factors";

  /* for each exercise  , find first coupon whose pay date > settlement date */
  i = 0;
  for (j = 0; j <= nex; j++) {
    while (tpay[i] <= t_start[j])
      i++;
    ifirst[j] = i;
  }

  /* calculate pv of payments */
  for (i = 0; i < npay; i++) {
    df = swp_f_df(0.0, tpay[i], ycName);
    if (df == SRT_DF_ERROR)
      return (error);
    pvpay[i] = fixpymnt[i] * df;
    /* Slight fix: coupons can prove negative (ZC option - funding) */
  }

  /* calculate pv of premium */
  for (j = 0; j <= nex; j++) {
    df = swp_f_df(tnow, t_start[j], ycName);
    prem[j] = origprem[j] * df;
    if (prem[j] < 0)
      return (error);
    df = swp_f_df(tnow, tpay[ifirst[j]], ycName);
    redpay[j] = redpay[j] * df;
    if (redpay[j] < 0)
      return (error);
  }

  error = NULL;
  return (error);
}

/**************************************************************************/
/* compute intrinsic value of mid-atlantic */

static double comp_intrin_val(long nex, long ifirst[], double redpay[],
                              long npay, double pvpay[]) {
  long ilast, j, i;
  double sum, temp, intrin_val;

  intrin_val = 0.;
  sum = 0.;

  ilast = npay;
  for (j = nex; j >= 0; j--) {
    for (i = ifirst[j]; i < ilast; i++)
      sum = sum + pvpay[i];
    ilast = ifirst[j];
    temp = sum - prem[j] - redpay[j];
    if (pr != 0)
      temp = -temp;
    if (intrin_val < temp)
      intrin_val = temp;
  }
  return (intrin_val);
}

/**************************************************************************/
/* choose dates for reference swaptions or caplets */
static void choose_dates(Date tnow, Date tstart, long *nptr, long npay,
                         Date tpay[], long nex, Date t_ex[], int icmpnd,
                         long num_months, BusDayConv natural_conv) {
  long i, n, nmin, start_ex;
  long num_fudge_minus, num_fudge_plus;
  /*	long  j  , nmax  , num_fudge;	*/
  Date tmin;
  /*	Date  tmaybe  , tpaymin; */

  nmin = ((icmpnd == 1) ? 2 : 1); /* if icmpnd=1  , min tenor is two periods */
  num_fudge_minus = 7;            /* num_fudge are the number of days we can */
  num_fudge_plus = 7;  /*	alter the schedule (after 1rst exer) to match */
  if (num_months >= 6) /* schedule to pay dates */
    num_fudge_plus = 35;
  if (num_months >= 12)
    num_fudge_plus = 35;

  /*	num_fudge = ((num_months>6)	? 2:1);
      if (tpay[npay-1] <= (tnow+20))
      {   for (j=0; j<=nmin; j++)
              s[j] = add_unit(tnow+20  , ((int) (num_months*j))  ,SRT_MONTH  ,
     natural_conv); n = nmin;
      }
      else
      {
  */
  /* first reference swaption start date will be before tmin  ,
          but no more than one period before tmin */

  start_ex = 0;
  tmin = add_unit(tnow, 7, SRT_DAY, natural_conv);
  if ((t_ex[0] < tmin) && (nex > 0)) {
    start_ex = 1;
  }
  /*        tmaybe = add_unit(tnow  , num_months  , SRT_MONTH  , natural_conv) +
     14; if (tmin < tmaybe)	*/	  /* Avoids calib. on instrmnts with less than 14 days to exer */
  /*           tmin = tmaybe;

          n = 0;
          tpaymin = tpay[npay-1];
          while (tpaymin >= tmin)
          {   n++;
              tpaymin = add_unit(tpay[npay-1]  , -((int)(num_months*n))  ,
     SRT_MONTH  , natural_conv);
          } */                                       /* tpaymin is first "paydate" <  tmin  */
  /*
                  nmax = n;
          if (nmax > (MAX_REF_DATES - 1))
              nmax = MAX_REF_DATES-1;
          for (j=0; j<=nmax; j++)
              s[j] = add_unit(tpay[npay-1]  , -((int)(num_months*(n-j)))  ,
     SRT_MONTH  , natural_conv); n = nmax;
  */
  i = 0;
  n = 0;
  s[0] = t_ex[start_ex];
  while (((s[n] < tpay[npay - 1]) || (n < nmin)) && (n < MAX_REF_DATES)) {
    n++;
    s[n] = add_unit(s[n - 1], num_months, SRT_MONTH, natural_conv);
    /* make date fall on pay date if possible */
    while (i < (npay - 1) &&
           tpay[i] <= s[n] - num_fudge_minus) /* PH 2/20/98 trivial change */
      i++;
    if ((tpay[i] >= s[n] - num_fudge_minus) &&
        (tpay[i] <= s[n] + num_fudge_plus))
      s[n] = tpay[i];
  }

  /*    }
   */
  /* OVE fix: we do not need to add dates:
          we want the last European to be priced properly	 */
  /* PH: 2/18/98: ensure that we have one complete swaption or caplet */

  /*    if (n<nmin)
      {   for (j=n+1; j<=nmin && j < MAX_REF_DATES; j++)
              s[j] = add_unit(t_ex[start_ex]  , num_months*n  , SRT_MONTH  ,
     natural_conv); n=nmin;
      }
  */
  /* make s[j] for j>=1 fall on a pay date if possible */
  /*   i=0;
      for (j=1; j<=n; j++)
      {   while (i<(npay-1) && tpay[i]<=s[j])
              i++;
          if (tpay[i]<=(s[j]+7) && tpay[i]>=(s[j]-7))
              s[j]=tpay[i];
      }
  */

  *nptr = n;
  return;
}

/**************************************************************************/
/* use three point smoothing to prevent weird deals from defeating code */

static void smooth_strikes(long nlast, double smth) {
  double rlast_lg, rlast_sh, rtempsh0, rtemplg0, rtemp;
  long j;

  if ((smth < 0) || (smth > 1))
    smth = 1.; /* smth=1 is full smoothing; smth=0 is no smoothing */

  /* OVE fix: since we do not add dates make a check if last one */
  /* PH 2/18/98 If nlast < 1 (one period)  , there is no smoothing to be done */
  if (nlast < 1)
    return;

  else if (nlast == 1) {
    rlast_lg = (.75 * R_theor_lg[1] + .25 * R_theor_lg[0]) * smth +
               (1. - smth) * R_theor_lg[1];
    rlast_sh = (.75 * R_theor_sh[1] + .25 * R_theor_sh[0]) * smth +
               (1. - smth) * R_theor_sh[1];
    R_theor_lg[0] = (.75 * R_theor_lg[0] + .25 * R_theor_lg[1]) * smth +
                    (1. - smth) * R_theor_lg[0];
    R_theor_sh[0] = (.75 * R_theor_sh[0] + .25 * R_theor_sh[1]) * smth +
                    (1. - smth) * R_theor_sh[0];
    R_theor_lg[1] = rlast_lg;
    R_theor_sh[1] = rlast_sh;
    return;
  } else
    ;
  /* Main case: nlast>1 */
  rlast_lg = (.75 * R_theor_lg[nlast] + .50 * R_theor_lg[nlast - 1] -
              .25 * R_theor_lg[nlast - 2]) *
                 smth +
             (1. - smth) * R_theor_lg[nlast];
  rlast_sh = (.75 * R_theor_sh[nlast] + .50 * R_theor_sh[nlast - 1] -
              .25 * R_theor_sh[nlast - 2]) *
                 smth +
             (1. - smth) * R_theor_sh[nlast];

  rtempsh0 = R_theor_sh[0];
  rtemplg0 = R_theor_lg[0];
  for (j = 0; j <= (nlast - 1); j++) {
    rtemp = R_theor_sh[j];
    R_theor_sh[j] =
        smth * 0.25 * (rtempsh0 + 2.0 * R_theor_sh[j] + R_theor_sh[j + 1]) +
        (1. - smth) * R_theor_sh[j];
    rtempsh0 = rtemp;

    rtemp = R_theor_lg[j];
    R_theor_lg[j] =
        smth * 0.25 * (rtemplg0 + 2.0 * R_theor_lg[j] + R_theor_lg[j + 1]) +
        (1. - smth) * R_theor_lg[j];
    rtemplg0 = rtemp;
  }
  R_theor_sh[nlast] = rlast_sh;
  R_theor_lg[nlast] = rlast_lg;

  return;
}

/*************************************************************************/
static String Olivier_lg(long n, long nlast, int icmpnd, long nex, long npay,
                         long ifirst[], Date t_start[], double pvpay[],
                         double redpay[])

/* Strikes for long reference swaptions are set to be at the money where
   the swaps underlying the midatlantic are at the money. */
{
  double ratio[MAX_NUM_DATES];
  double sum, rat;
  long i, j, ilast;
  error = NULL;

  /* compute ratio = pv of fixed leg / pv of premium */
  ilast = npay;
  sum = 0.;
  for (j = nex; j >= 0; j--) {
    for (i = ifirst[j]; i < ilast; i++)
      sum = sum + pvpay[i];
    ilast = ifirst[j];
    ratio[j] = (sum - redpay[j]) / prem[j];
  }

  /* interpolate on the ratios to find value at s[j]  , and choose the
     interest rate of the long swaption so that it has the same ratio  */

  sum = 0.;
  for (j = nlast + 1; j < n; j++)
    sum = sum + aB[j + 1];

  i = nex;
  for (j = nlast; j >= 0; j--) {
    sum = sum + aB[j + 1];
    while (i > 0 && t_start[i] > s[j])
      i--;
    if (i == nex)
      rat = ratio[nex];
    else if (t_start[i] > s[j])
      rat = ratio[0];
    else
      rat = (ratio[i] * ((double)(t_start[i + 1] - s[j])) +
             ratio[i + 1] * ((double)(s[j] - t_start[i]))) /
            ((double)(t_start[i + 1] - t_start[i]));
    R_theor_lg[j] = (rat * B[j] - B[n]) / sum;
  }
  /* quick fix PH Feb. 10 1997 */
  /* eliminates slight arbitrage for in the money deals
          when the paydates extend more than a couple years past the exercise
     dates */

  rat = R_theor_lg[0];
  for (j = 0; j <= nlast; j++) {
    if (s[j] <= t_start[nex])
      rat = R_theor_lg[j]; /* Good: R was found by interpolation */
    else                   /* Bad: R was found by extrapolation */
      R_theor_lg[j] = rat; /* Replace with last interpolated value */
  }

  return (error);
}

static String Olivier_sh(long n, long nlast, long nex, long npay, long ifirst[],
                         Date t_start[], double pvpay[], double redpay[],
                         int icmpnd, int fix_tau_flag)

/* Strikes for reference swaptions are set to be at the money with
   the delta irr  of the underlying swaps */
{
  double ratio[MAX_NUM_DATES];
  double sum, rat, rat1, rat2, prem1, prem2;
  long i, i1, i2, j, ilast, jlast, jend, delay;
  error = NULL;

  delay = icmpnd + 1;
  /* compute ratio = pv of fixed leg / pv of premium */
  ilast = npay;
  sum = 0.;
  for (j = nex; j >= 0; j--) {
    for (i = ifirst[j]; i < ilast; i++)
      sum = sum + pvpay[i];
    ilast = ifirst[j];
    ratio[j] = (sum - redpay[j]) / prem[j];
  }

  /* interpolate on ratio[j] to find its value at s[j] and s[j+delay]  , and
     choose the interest rate of the short swaption so that it is equivalent to
     buying a long swaption at s[j] with the first ratio and selling one at
     s[j+delay] at the second ratio  */

  i1 = nex;
  i2 = nex;
  for (j = nlast; j >= 0; j--) {
    while (i1 > 0 && t_start[i1] > s[j])
      i1--;
    if (i1 == nex) {
      rat1 = ratio[nex];
      prem1 = origprem[nex];
    } /* interpolation */
    else if (t_start[i1] > s[j]) {
      rat1 = ratio[0];
      prem1 = origprem[0];
    } else {
      rat1 = (ratio[i1] * ((double)(t_start[i1 + 1] - s[j])) +
              ratio[i1 + 1] * ((double)(s[j] - t_start[i1]))) /
             ((double)(t_start[i1 + 1] - t_start[i1]));
      prem1 = (origprem[i1] * ((double)(t_start[i1 + 1] - s[j])) +
               origprem[i1 + 1] * ((double)(s[j] - t_start[i1]))) /
              ((double)(t_start[i1 + 1] - t_start[i1]));
    }

    jlast = j + delay;
    if (jlast > n)
      jlast = n;
    while (i2 > 0 && t_start[i2] > s[jlast])
      i2--;
    if (i2 == nex) {
      rat2 = ratio[nex];
      prem2 = origprem[nex];
    } /* interpolation */
    else if (t_start[i2] > s[jlast]) {
      rat2 = ratio[0];
      prem2 = origprem[0];
    } else {
      rat2 = (ratio[i2] * ((double)(t_start[i2 + 1] - s[jlast])) +
              ratio[i2 + 1] * ((double)(s[jlast] - t_start[i2]))) /
             ((double)(t_start[i2 + 1] - t_start[i2]));
      prem2 = (origprem[i2] * ((double)(t_start[i2 + 1] - s[jlast])) +
               origprem[i2 + 1] * ((double)(s[jlast] - t_start[i2]))) /
              ((double)(t_start[i2 + 1] - t_start[i2]));
    }

    /* compute strike of swaption */
    sum = 0.;
    jend = jlast;
    for (i = (j + 1); i <= jend; i++)
      sum = sum + aB[i];
    R_theor_sh[j] =
        (rat1 * B[j] - B[jlast] * (1.0 + prem2 * (rat2 - 1.0) / prem1)) / sum;
  }
  /* quick fix PH Feb. 10 1997 */
  /* eliminates slight arbitrage for in the money deals
          when the paydates extend more than a couple years past the exercise
     dates */

  /*	PH 2/18/98 this doesn't seem to be necessary */

  if (fix_tau_flag != 1) {
    rat = R_theor_lg[0]; /* if not enough exercise dates  , use R_lg for R_sh */
    if (s[0] <= t_start[nex] && t_start[nex] <= s[delay]) {
      rat = ((s[delay] - t_start[nex]) *
                 R_theor_lg[0] /* feather shorts if more exer. dates */
             + (t_start[nex] - s[0]) * R_theor_sh[0]) /
            (s[delay] - s[0]);
    }
  } else {
    rat = R_theor_sh[0];
  }

  for (j = 0; j <= nlast; j++) {
    if (s[j + delay] <= t_start[nex])
      rat = R_theor_sh[j]; /* Good: R was found by interpolation */
    else                   /* Bad: R was found by extrapolation */
      R_theor_sh[j] = rat; /* Replace with last interpolated value */
  }
  return (error);
}

/**************************************************************************/
/* calculate the value of the reference swaption */

static double Vjn(double par, int whichpar) {
  double ans, sum, sqz, target1, phi;
  long i;
  int check;

  /* update zeta[jst]  */
  if (whichpar == 0)
    zeta[jst] = par * par;
  else if (whichpar == 1)
    g[jst] = par;

  /* find phi at the money */
  target1 = 0.;
  phi = phiact;
  ScNewt(&phi, &check, target1, pay_payer, pay_payer_der);
  if (check != 0) {
    error = "no phi critical found";
    phi = 0.;
  }
  phiact = phi;
  sqz = sqrt(zeta[jst]);

  /* find the LGM value */
  sum = 0.;
  for (i = jst + 1; i < nst; i++)
    sum = sum + aB[i] * norm_accurate(phi / sqz - g[i] * sqz);
  ans = sum * Ract +
        (B[nst] + aB[nst] * Ract) * norm_accurate(phi / sqz - g[nst] * sqz) -
        B[jst] * norm_accurate(phi / sqz - g[jst] * sqz);
  if (pr == 1) {
    sum = 0.;
    for (i = jst + 1; i < nst; i++)
      sum = sum + aB[i];
    ans = ans - Ract * sum - B[nst] - Ract * aB[nst] + B[jst];
  }
  return (ans);
}

/* Calculate the derivative of the above function with respect to
   sqz=sqrt(zeta[j]). phiact must be set to phi at the money.*/

static double Vjn_der(double par, int whichpar) {
  double sum, sqz, tmp, ans;
  long i;

  sqz = sqrt(zeta[jst]);
  if (whichpar == 1) {
    tmp = phiact / sqz - g[jst] * sqz;
    ans = B[jst] * sqz / sq2pi * exp(-tmp * tmp / 2.);
    return (ans);
  }

  tmp = phiact / sqz - g[nst] * sqz;
  sum = (B[nst] + Ract * aB[nst]) * (g[jst] - g[nst]) * exp(-tmp * tmp / 2.);

  for (i = jst + 1; i < nst; i++) {
    tmp = phiact / sqz - g[i] * sqz;
    sum = sum + Ract * aB[i] * (g[jst] - g[i]) * exp(-tmp * tmp / 2.);
  }

  ans = sum / sq2pi;
  return (ans);
}

/* calculate value of reference swaption on kappa = constant curve */

static double Vjn1(double par, int whichpar) {
  double ans, sum, sqz, target1, phi;
  long i;
  int check;

  /* update zeta[jst]  */
  zeta[jst] = par * par;
  g[jst] = g[jst + 1] + kappa / par;

  /* find phi at the money */
  target1 = 0.;
  phi = phiact;
  ScNewt(&phi, &check, target1, pay_payer, pay_payer_der);
  if (check != 0) {
    error = "no phi critical found";
    phi = 0.;
  }
  phiact = phi;
  sqz = sqrt(zeta[jst]);
  /* find the LGM value */
  sum = 0.;
  for (i = jst + 1; i < nst; i++)
    sum = sum + aB[i] * norm_accurate(phi / sqz - g[i] * sqz);
  ans = sum * Ract +
        (B[nst] + aB[nst] * Ract) * norm_accurate(phi / sqz - g[nst] * sqz) -
        B[jst] * norm_accurate(phi / sqz - g[jst] * sqz);

  if (pr == 1) {
    sum = 0.;
    for (i = jst + 1; i < nst; i++)
      sum = sum + aB[i];
    ans = ans - Ract * sum - B[nst] - Ract * aB[nst] + B[jst];
  }

  return (ans);
}

/* Calculate the derivative of the above function with respect to
   sqz=sqrt(zeta[j]). phiact must be set to phi at the money.*/

static double Vjn1_der(double par, int whichpar) {
  double sum, sqz, tmp, ans;
  long i;

  sqz = sqrt(zeta[jst]);
  tmp = phiact / sqz - g[nst] * sqz;

  sum =
      (B[nst] + Ract * aB[nst]) * (g[jst + 1] - g[nst]) * exp(-tmp * tmp / 2.);
  for (i = jst + 1; i < nst; i++) {
    tmp = phiact / sqz - g[i] * sqz;
    sum = sum + Ract * aB[i] * (g[jst + 1] - g[i]) * exp(-tmp * tmp / 2.);
  }
  ans = sum / sq2pi;

  return (ans);
}

/*************************************************************************/
/* par swap rates and variants */
/* calculate the par swap rate in the future (at s[jst] for swap
   starting at s[jst] and ending at s[nst]  , given phi and zeta at s[jst] */

static double Rpar_jn(double ph) {
  double sum, ans;
  long i;
  sum = 0.;
  for (i = jst + 1; i <= nst; i++)
    sum = sum + aB[i] * exp(g[i] * (ph - .5 * g[i] * zeta[jst]));

  ans = (B[jst] * exp(g[jst] * (ph - .5 * g[jst] * zeta[jst])) -
         B[nst] * exp(g[nst] * (ph - .5 * g[nst] * zeta[jst]))) /
        sum;

  return (ans);
}

/* calculate the par swap rate (today) for swap beginning at s[jst]
   and ending at s[nst]        */

static double Rbs_par_jn() {
  double sum, ans;
  long i;

  sum = 0.;
  for (i = jst + 1; i <= nst; i++)
    sum = sum + aB[i];

  ans = (B[jst] - B[nst]) / sum;

  return (ans);
}

/************************************************************************/
/* Calculate the Black-Scholes price of swaption beginning at s[jst]
   and ending at s[nst]. Fixed rate r  , BS volatility sig  , expiry is
   s[jst]  , and pr=1 for a payer and pr=0 for reciever.                */

static double Vbs_jn(double r, double sig, SrtDiffusionType srt_vol_type,
                     int spot_lag) {
  double ans, rswp, chi, var;
  double years_to_expiry, spot_lag_years;

  rswp = Rbs_par_jn();
  spot_lag_years = s[jst] - add_unit(s[jst], -spot_lag, SRT_BDAY, SUCCEEDING);
  years_to_expiry = sy[jst] - spot_lag_years / 365.0;

  if (years_to_expiry < 1. / 365.)
    years_to_expiry = 1. / 365.;

  var = sig * sqrt(years_to_expiry);
  /* Lognormal Black-Scholes */
  if (srt_vol_type == SRT_LOGNORMAL) {
    chi = (log(r / rswp)) / var + .5 * var;
    ans = (B[jst] - B[nst]) *
          ((r / rswp) * norm_accurate(chi) - norm_accurate(chi - var));
  } else
      /* Normal Black-Scholes */
      if (srt_vol_type == SRT_NORMAL) {
    chi = (r - rswp) / var;
    ans = (B[jst] - B[nst]) / rswp * var *
          (chi * norm_accurate(chi) + gauss(chi));
  }
  if (pr == 1)
    ans = ans + (B[jst] - B[nst]) * (1. - r / rswp);

  return (ans);
}

/************************************************************************/
/* Values of swaps in the future within the LGM framework  , and
   similar functions */
/* Calculate the value of vanilla swap (payer) at time s[jst]  ,
   given phi. Swap starts at s[jst]  , ends at s[nst]  , and has
   interest rate Ract   */

static double pay_payer(double ph) {
  double sum, ans;
  long i;

  sum = 0.;
  for (i = jst + 1; i < nst; i++)
    sum = sum + aB[i] * exp(-(g[jst] - g[i]) *
                            (ph - .5 * (g[i] + g[jst]) * zeta[jst]));

  ans = sum * Ract - B[jst] +
        (B[nst] + Ract * aB[nst]) *
            exp(-(g[jst] - g[nst]) * (ph - .5 * (g[jst] + g[nst]) * zeta[jst]));

  return (ans);
}

/* Calculate the derivative of the above function with respect to phi */

static double pay_payer_der(double ph) {
  double sum, ans;
  long i;

  sum = 0.;
  for (i = jst + 1; i < nst; i++)
    sum = sum -
          aB[i] * (g[jst] - g[i]) *
              exp(-(g[jst] - g[i]) * (ph - .5 * (g[i] + g[jst]) * zeta[jst]));

  ans = sum * Ract -
        (B[nst] + Ract * aB[nst]) * (g[jst] - g[nst]) *
            exp(-(g[jst] - g[nst]) * (ph - .5 * (g[jst] + g[nst]) * zeta[jst]));

  return (ans);
}

/*********************************************************************/
/* Calculate the LGM value (val) of vanilla swaption starting at s[jst]  ,
   ending at s[nst]  , given zeta[jst] and all g[i] for i>=jst.
   Also calculate the derivatives wrt sqrt(zeta[jst])  , g[j]  , and g[j+1]
     (v1  , v2  , and v3).
    phiact = initial guess for "at the money" phi needs to be set
    (at zero if nothing else)  , and the fixed rate Ract needs to be set */

static void Vjn_all(double *val, double *v1, double *v2) {
  double sum, sum1, sum2, sqz, most, phi, target1, tmp;
  long i;
  int check;

  sqz = sqrt(zeta[jst]);
  check = 0;

  /* find phi at the money */
  target1 = 0.;
  phi = phiact;
  ScNewt(&phi, &check, target1, pay_payer, pay_payer_der);
  if (check != 0) {
    error = "no phi critical found";
    phi = phiact;
  }
  phiact = phi;

  /* find the LGM values and derivatives */
  sum = sum1 = 0.;
  for (i = jst + 1; i < nst; i++) {
    tmp = phi / sqz - g[i] * sqz;
    sum = sum + aB[i] * norm_accurate(tmp);
    sum1 = sum1 + aB[i] * (g[jst] - g[i]) * exp(-tmp * tmp / 2.);
  }

  tmp = phi / sqz - g[nst] * sqz;
  most = (B[nst] + Ract * aB[nst]) * exp(-tmp * tmp / 2.);
  sum1 = (Ract * sum1 + most * (g[jst] - g[nst])) / sq2pi;
  sum = Ract * sum + (B[nst] + Ract * aB[nst]) * norm_accurate(tmp);

  tmp = phi / sqz - g[jst] * sqz;
  sum = sum - B[jst] * norm_accurate(tmp);
  sum2 = B[jst] * sqz / sq2pi * exp(-tmp * tmp / 2.);

  if (pr == 1) /* use parity for payer */
  {
    sum = sum + B[jst] - B[nst] - Ract * aB[nst];
    for (i = jst + 1; i < nst; i++)
      sum = sum - Ract * aB[i];
  }

  *val = sum;
  *v1 = sum1;
  *v2 = sum2;
  return;
}

/**************************************************************************/
/* Calculate the swap rates that are FAR_OUT std deviations from par
   for swaptions that start at s[jst] and end at s[jst+2] and s[nst] */

static String R_restrictions(String (*GetVol)(), SrtDiffusionType srt_vol_type,
                             long n, long nlast, int icmpnd,
                             SrtCompounding compounding, SrtBasisCode basis,
                             int fix_tau_flag, int usecaps) {
  double rbs, sig;
  long j, dj;
  SwapDP sdp;
  double fwd;
  Err err;

  error = NULL;
  dj = 2;
  if (icmpnd == 0)
    dj = 1;

  for (j = 0; j <= nlast; j++) {
    if ((fix_tau_flag != 1) ||
        (usecaps == 1)) /*get par swap rate for short swaption */
    {
      jst = j;
      nst = j + dj;
      rbs = Rbs_par_jn();
      /* PH: 2/18/98 Will this work for caplets? If so  , uncomment
                      if (err = srt_f_setSwapDP(s[j]  , s[j+dj]  , compounding
         , basis  , &sdp)) return err;

                      err = swp_f_ForwardRate_SwapDP( &sdp  , s_yldcrv  , "CASH"
         , &fwd); if (err) return err;
      */
      /*   error = (*GetVol)(s[j]  ,s[j+dj]-s[j]  ,fwd  ,&sig); */ /* Get BS
                                                                      volatility
                                                                    */
      error = (*GetVol)(s[j], s[j + dj], rbs, 0, &sig); /* Get BS volatility */
      if (error != NULL) {
        error = "vol failed in R_restriction";
        return (error);
      }

      /* extreme strikes for short swaptions */
      if (srt_vol_type == SRT_LOGNORMAL) {
        Rmin_sh[j] =
            rbs * exp(-sig * sqrt(sy[j]) * FAR_OUT - .5 * sig * sig * sy[j]);
        Rmax_sh[j] =
            rbs * exp(+sig * sqrt(sy[j]) * FAR_OUT - .5 * sig * sig * sy[j]);
      } else {
        Rmin_sh[j] = rbs - sig * sqrt(sy[j]) * FAR_OUT;
        Rmax_sh[j] = rbs + sig * sqrt(sy[j]) * FAR_OUT;
      }
    }
    if ((fix_tau_flag != 1) || (usecaps != 1)) {
      jst = j;
      nst = n; /*get par swap rate for long swaption */
      rbs = Rbs_par_jn();

      if (err = swp_f_setSwapDP(s[j], s[n], compounding, basis, &sdp))
        return err;

      err = swp_f_ForwardRate_SwapDP(&sdp, s_ycname, "CASH", &fwd);
      if (err)
        return err;

      error = (*GetVol)(s[j], s[n], fwd, 0, &sig); /* Get BS volatility*/
      if (error != NULL) {
        error = "vol failed in R_restriction";
        return (error);
      }

      /* extreme strikes for long swaptions */
      if (srt_vol_type == SRT_LOGNORMAL) {
        Rmin_lg[j] =
            rbs * exp(-sig * sqrt(sy[j]) * FAR_OUT - .5 * sig * sig * sy[j]);
        Rmax_lg[j] =
            rbs * exp(+sig * sqrt(sy[j]) * FAR_OUT - .5 * sig * sig * sy[j]);
      } else {
        Rmin_lg[j] = rbs - sig * sqrt(sy[j]) * FAR_OUT;
        Rmax_lg[j] = rbs + sig * sqrt(sy[j]) * FAR_OUT;
      }
    }
  }
  return (error);
}

/**************************************************************************/
/*MID ATLANTIC EVALUATOR
1. Calculate the value of a mid-atlantic using the convolution method
2. Find critical phi values where optimal exercise occurs                 */

static String
mid_atl_eval(Date tnow, /* calculation date */
                        /* model parameters*/
             long n,    /* zeta[0  ,1  ,...  ,n]  , g[0  ,1  ,...  ,n]
                             at s[0  ,1  ,...  ,n]               */
                        /* market place */
             SrtCurvePtr m,
             /* exercise info */
             long nex, Date t_ex[], /* exercise dates t_ex[0  ,...  ,nex]   */
             Date t_start[],        /* start dates t_start[0  ,...  ,nex]   */
                                    /*payment info */
             long npay,
             Date tpay[],     /* pay dates tpay[0  ,1  ,...  ,npay-1]    */
             double pvpay[],  /* pv of pymnts pvpay[0  ,1  ,...  ,npay-1] */
             double redpay[], /* pv of 1rst payment reduction [0  ,...nex] */
             long ifirst[],   /* ifirst[0  ,...  ,nex] */
                              /*answer*/
             double *LGMvalue,
             /*strikes*/
             double phi_cr[]) /* critical phi values for strikes    */

{
  /* Declarations */
  double *wgt, *wgt2, *Q0, *Q1;
  double *zeta_ex, *G_start, *G_pay, *eff_pay, *postfac, *effred;
  double phi, phiwidth, phimin, dphi;
  double qmin, dq, offset, scale;
  double var, alfa, beta, vact;
  double Vs1, Vs2;
  double exk, rik, qcr;
  double v_cr, exk_cr, rik_cr, rik2_cr, exi_cr, ri_cr;
  long i, j, k, nq, nq2, nphi, j0, k0;
  long kk, klast, kklast, ik, ik_cr, ik2_cr, ii_cr;
  int check;
  double phi_neg, phi_pos, phi_width_narrow;
  double exk_ratio;

  double coarse_pv;
  double fine_pv;

  /*Initialize */
  error = NULL;
  nq = NQMAX;
  nphi = NDIV;
  nq = (nq + 3) / 4;
  nq2 = 2 * nq;
  nq = 4 * nq;
  nphi = (nphi + 3) / 4;
  nphi = 4 * nphi;
  Vs1 = -33.3333e+257;
  Vs2 = +33.3333e+257;

  /*STEP ONE: Allocate space */
  wgt = NULL;
  wgt2 = NULL;
  Q0 = NULL;
  Q1 = NULL;
  zeta_ex = NULL;
  G_start = NULL;
  G_pay = NULL;
  eff_pay = NULL;
  postfac = NULL;
  effred = NULL;

  zeta_ex = (double *)srt_calloc(nex + 1,
                                 sizeof(double)); /* zeta_ex[0  ,...  ,nex]   */
  G_start = (double *)srt_calloc(nex + 1,
                                 sizeof(double)); /* G_start[0  ,...  ,nex]   */
  G_pay =
      (double *)srt_calloc(npay, sizeof(double)); /* G_pay[0  ,...  ,npay-1]  */
  wgt = (double *)srt_calloc(nq + 1, sizeof(double)); /* wgt[0  ,...  ,nq] */
  wgt2 =
      (double *)srt_calloc(nq2 + 1, sizeof(double)); /* wgt[0  ,...  ,nq2] */
  Q0 =
      (double *)srt_calloc(nphi + 1, sizeof(double)); /* Q0[0  ,...  ,nphi] */
  Q1 =
      (double *)srt_calloc(nphi + 1, sizeof(double)); /* Q1[0  ,...  ,nphi] */
  eff_pay = (double *)srt_calloc(npay * (nphi + 1), sizeof(double));
  /* eff_pay[k][i]  , 0<=k<=nphi  , 0<=i<=npay-1 */
  effred = (double *)srt_calloc((nex + 1) * (nphi + 1), sizeof(double));
  /* effred[k][j]  , 0<=k<=nphi  , 0<=j<=nex     */
  postfac = (double *)srt_calloc(npay * (nex + 1), sizeof(double));
  /* postfac[j][i]  , 0<=j<=nex  , 0<=i<=npay-1  */

  if (wgt == NULL || wgt2 == NULL || Q0 == NULL || Q1 == NULL ||
      zeta_ex == NULL || G_start == NULL || G_pay == NULL || eff_pay == NULL ||
      postfac == NULL || effred == NULL) {
    error = "alloc. failed in evaluator";
    return (error);
  }

  /* STEP TWO: Set up numerical integration */
  /* Interpolate G and zeta values at needed dates */

  interpolate_par(tnow, n, nex, t_ex, zeta_ex, t_start, G_start, npay, tpay,
                  G_pay);

  genweights(nq, STENCIL, wgt);   /* generate integration weights */
  genweights(nq2, STENCIL, wgt2); /* generate integration weights */

  /* Find initial strike to center grid on */
  j = nex;
  init_point(&phi, &check, prem[j], ifirst[j], npay, pvpay, redpay[j], G_pay,
             G_start[j], zeta_ex[j]);
  if (check != 0) {
    error = "can't find initial strike";
    ERROR_CHECK1;
  }

  phi_cr[nex] = phi;
  phi_width_narrow = WIDTH * sqrt(zeta_ex[nex]);
  phi_neg = phi < 0.0 ? phi : 0;
  phi_pos = phi > 0.0 ? phi : 0;
  phimin = phi_neg - phi_width_narrow;
  phiwidth = (phi_pos - phi_neg) / 2.0 + phi_width_narrow;
  dphi = 2. * phiwidth / ((double)nphi); /* phi = phimin + k*dphi */

  qmin = -STENCIL;
  dq = -2. * qmin / ((double)nq); /* q = qmin + i*dq */

  /* Setup payoff matrix (trade space for time!) */
  for (k = 0; k <= nphi; k++) {
    phi = (2 * k - nphi) * phi_width_narrow;
    phi += (nphi - k) * phi_neg + k * phi_pos;
    phi /= nphi;
    k0 = k * npay;
    for (i = 0; i < npay; i++)
      eff_pay[k0 + i] = pvpay[i] * exp(G_pay[i] * phi);
    k0 = k * (nex + 1);
    for (j = 0; j <= nex; j++) {
      i = ifirst[j];
      effred[k0 + j] =
          redpay[j] * exp(G_pay[i] * (phi - .5 * G_pay[i] * zeta_ex[j]));
    }
  }

  for (j = 0; j <= nex; j++) {
    j0 = j * npay;
    for (i = 0; i < npay; i++)
      postfac[j0 + i] = exp(-.5 * G_pay[i] * G_pay[i] * zeta_ex[j]);
  }

  /* Initialize Q0  , Q1 */
  for (kk = 0; kk <= nphi; kk++) {
    Q0[kk] = Q1[kk] = 0.;
  }

  /* Step backwards from the last exercise date  , first replacing
     Q0 with the payoff when the payoff is larger  , then convolving
     against a Gaussian to get the values at the next previous
     exercise date                                              */

  for (j = nex; j >= 0; j--) {
    j0 = j * npay;

    /* replace with payoff for all phi if payoff is larger */
    klast = 0;
    kklast = 0;
    for (kk = 0; kk <= nphi; kk++) {
      k = (pr == 0 ? kk : (nphi - kk));
      phi = (2 * k - nphi) * phi_width_narrow;
      phi += (nphi - k) * phi_neg + k * phi_pos;
      phi /= nphi;
      /* compute payoff */
      vact = -effred[k * (nex + 1) + j];
      k0 = k * npay;
      for (i = ifirst[j]; i < npay; i++)
        vact = vact + eff_pay[k0 + i] * postfac[j0 + i];

      vact = vact -
             prem[j] * exp(G_start[j] * (phi - .5 * G_start[j] * zeta_ex[j]));
      if (pr != 0)
        vact = -vact;

      if (vact > Q0[k]) {
        klast = k; /* replace with payoff */
        kklast = kk;
        Vs1 = Q0[k];
        Q0[k] = vact;
      } else if (kk == kklast + 1)
        Vs2 = vact;
      else
        ;
      if (vact < 0.0)
        kk = nphi + 1; /*no more replacements */
      else
        ;
    }

    /* find exercise boundary */
    if (klast == 0 || klast == nphi) {
      if (j < nex)
        exk_cr = (phi_cr[j + 1] - phimin) / dphi;
      else
        exk_cr = -phimin / dphi;
    } else if (pr == 0)
      exk_cr = ((double)klast) +
               (Q0[klast] - Vs1) / (Q0[klast + 1] - Vs2 + Q0[klast] - Vs1);
    else
      exk_cr = ((double)klast) -
               (Q0[klast] - Vs1) / (Q0[klast - 1] - Vs2 + Q0[klast] - Vs1);
    phi_cr[j] = phimin + dphi * exk_cr;

    /* Evaluate mid-atlantic  , Q0  , at previous time step by integration */
    if (j > 0)
      var = sqrt(zeta_ex[j] - zeta_ex[j - 1]);
    else
      var = sqrt(zeta_ex[j]);
    if (var < 2.0e-5)
      var = 2.0e-5;

    /* Transfer values */
    for (k = 0; k <= nphi; k++) {
      Q1[k] = Q0[k];
      Q0[k] = 0.;
    }

    offset = var * qmin / dphi; /* phi_prime = phimin + dphi*k_prime */
    scale = var * dq / dphi;    /* k_prime = k + offset + scale*i      */

    /* At exercise boundary  , a break in slope occurs */
    /* Find value of payoff and k where the break occurs */

    v_cr = EXPAY(phi_cr[j], npay, ifirst[j], pvpay, redpay[j], prem[j], G_pay,
                 G_start[j], zeta_ex[j]);
    ik_cr = DTOL(exk_cr);
    rik_cr = ik_cr;
    ik2_cr = 2 * (DTOL(exk_cr / 2.));
    rik2_cr = ik2_cr;

    /* convolve against Gaussian with half the resolution  , and subtract
       1/3 of this value from 4/3 of the value at full resolution to get
       result with the O(h*h) error extrapolated away */

    exk_ratio = var / dphi * (double)qmin / (double)nq;
    for (k = 0; k <= nphi; k++)

    /* convolve with discretization 2dq and 2dphi */
    {
      coarse_pv = 0.0;
      fine_pv = 0.0;

      exi_cr = (exk_cr - (double)k) / (2. * scale) + ((double)nq) / 4.;
      qcr = qmin + 2.0 * dq * exi_cr;
      ii_cr = DTOL(exi_cr);
      ri_cr = (double)ii_cr;
      if (fabs(qcr) <= -qmin)
        beta = exp(-qcr * qcr / 2.) / (2. * sq2pi);
      else
        beta = 0.;

      for (i = 0; i <= nq2; i++) {
        exk = k + exk_ratio * (nq - 4 * i);
        ik = 2 * (DTOL(exk / 2.));
        rik = ik;
        alfa =
            iterp_integrand(exk, rik, ik, Q1, nphi, v_cr, exk_cr, rik2_cr, 2);

        coarse_pv = coarse_pv + wgt2[i] * alfa;
        /* Q0[k] = Q0[k] -  wgt2[i]*alfa/3.0; */ /* sum is -1/3 of the integral
                                                    evaluated */
        /*   on the coarse grid                  */

        if (i == ii_cr) /* correction for break in slope */
          coarse_pv = coarse_pv +
                      beta * 2. * dq * (v_cr - alfa) * (1. + ri_cr - exi_cr);
        /*    Q0[k] = Q0[k] - beta*2.*dq*(v_cr-alfa)*(1. + ri_cr - exi_cr)/3.;
         */
        if (i == (ii_cr + 1))
          coarse_pv =
              coarse_pv + beta * 2. * dq * (v_cr - alfa) * (exi_cr - ri_cr);
        /*    Q0[k] = Q0[k] - beta*2.*dq*(v_cr-alfa)*(exi_cr - ri_cr)/3.; */
      }

      /* convolve with discretization dq and dphi */
      exi_cr = 2.0 * exi_cr;
      ii_cr = DTOL(exi_cr);
      ri_cr = (double)ii_cr;

      for (i = 0; i <= nq; i++) {
        exk = k + exk_ratio * (nq - 2 * i);
        ik = DTOL(exk);
        rik = ik;

        alfa = iterp_integrand(exk, rik, ik, Q1, nphi, v_cr, exk_cr, rik_cr, 1);
        fine_pv = fine_pv + wgt[i] * alfa;
        /*Q0[k] = Q0[k] + 4.0* wgt[i]*alfa/3.0 ;* /* adds 4/3 the integral
         * evaluated on  */
        /* the fine grid                       */

        if (i == ii_cr) /* correction for break in slope */
          fine_pv = fine_pv + beta * dq * (v_cr - alfa) * (1. + ri_cr - exi_cr);
        /*   Q0[k] = Q0[k] + 4.0*beta*dq*(v_cr-alfa)*(1. + ri_cr - exi_cr)/3.;
         */
        if (i == (ii_cr + 1))

          if (i == (ii_cr + 1))
            fine_pv = fine_pv + beta * dq * (v_cr - alfa) * (exi_cr - ri_cr);
        /*    Q0[k] = Q0[k] + 4.0*beta*dq*(v_cr-alfa)*(exi_cr - ri_cr)/3.; */

      } /* end i loop */

      Q0[k] = -1.0 / 3.0 * coarse_pv + 4.0 / 3.0 * fine_pv;
    } /* end k loop */

  } /* end j step */

  /* Find value at tnow (phi=0) */
  exk = -phimin / dphi;
  ik = DTOL(exk);
  rik = ik;
  *LGMvalue = (1. + rik - exk) * Q0[ik] + (exk - rik) * Q0[ik + 1];
  ik = 2 * (DTOL(exk / 2.));
  rik = ik;
  *LGMvalue = 4. / 3. * (*LGMvalue) - (2. + rik - exk) * Q0[ik] / 6. -
              (exk - rik) * Q0[ik + 2] / 6.;
  /* Free everything and return */
  FreeMore(zeta_ex, G_start, G_pay, wgt, wgt2, Q0, Q1, eff_pay, postfac,
           effred);
  return (error);
}

/*************************************************************************/
/* interpolate the between Q1[ik] and Q1[ik+di] to get the integrand at exk */

static double iterp_integrand(double exk, double rik, long ik, double Q1[],
                              long nphi, double v_cr, double exk_cr,
                              double rik_cr, long di) {
  double alfa, edi;

  /* OVE fix: if interpolation gives negative prices  , set to 0 */
  edi = (double)di;
  if (exk <= 0.) {
    alfa = Q1[0] + (Q1[di] - Q1[0]) * exk / edi;
  } else if (exk >= nphi) {
    alfa = Q1[nphi] + (Q1[nphi] - Q1[nphi - di]) * (exk - (double)nphi) / edi;
  } else if (exk <= rik_cr || exk >= (rik_cr + edi)) {
    alfa = (Q1[ik] * (rik + edi - exk) + Q1[ik + di] * (exk - rik)) / edi;
  } else if (exk < exk_cr) {
    alfa = (v_cr * (exk - rik + .01) + Q1[ik] * (exk_cr - exk + .01)) /
           (exk_cr - rik + .02);
  } else {
    alfa =
        (v_cr * (edi + .01 + rik - exk) + Q1[ik + di] * (exk - exk_cr + .01)) /
        (rik + edi - exk_cr + .02);
  }

  if (alfa < 0.0)
    alfa = 0.0;

  return (alfa);
}

/*************************************************************************/
/* interpolate to zeta & G values at needed dates;
   store G(t_start[j]) in G_start[0  ,...  ,nex];
   store G(tpay[i]) in G_pay[0  ,1  ,...  ,npay-1];
   store zeta(t_ex[j]) in zeta_ex[0  ,1  ,...  ,nex]   */

static void interpolate_par(Date tnow, long n, long nex, Date t_ex[],
                            double zeta_ex[], Date t_start[], double G_start[],
                            long npay, Date tpay[], double G_pay[]) {
  long i, j;

  j = 0;
  for (i = 0; i <= npay - 1; i++) {
    while (s[j] < tpay[i] && j < n)
      j++;
    if (s[j] < tpay[i])
      G_pay[i] = g[n] + (g[n] - g[n - 1]) * (tpay[i] - s[n]) /
                            ((double)(s[n] - s[n - 1]));
    else if (j == 0)
      G_pay[i] =
          g[0] + (g[1] - g[0]) * (tpay[i] - s[0]) / ((double)(s[1] - s[0]));
    else
      G_pay[i] = (g[j] * (tpay[i] - s[j - 1]) + g[j - 1] * (s[j] - tpay[i])) /
                 ((double)(s[j] - s[j - 1]));
  }

  j = 0;
  for (i = 0; i <= nex; i++) {
    while (s[j] < t_ex[i] && j < n)
      j++;
    if (s[j] < t_ex[i])
      zeta_ex[i] = zeta[n] + (zeta[n] - zeta[n - 1]) * (t_ex[i] - s[n]) /
                                 ((double)(s[n] - s[n - 1]));
    else if (j == 0)
      zeta_ex[i] = zeta[0] * (t_ex[i] - tnow) / ((double)(s[0] - tnow));
    else
      zeta_ex[i] =
          (zeta[j] * (t_ex[i] - s[j - 1]) + zeta[j - 1] * (s[j] - t_ex[i])) /
          ((double)(s[j] - s[j - 1]));
  }

  j = 0;
  for (i = 0; i <= nex; i++) {
    while (s[j] < t_start[i] && j < n)
      j++;
    if (s[j] < t_start[i])
      G_start[i] = g[n] + (g[n] - g[n - 1]) * (t_start[i] - s[n]) /
                              ((double)(s[n] - s[n - 1]));
    else if (j == 0)
      G_start[i] =
          g[0] + (g[1] - g[0]) * (t_start[i] - s[0]) / ((double)(s[1] - s[0]));
    else
      G_start[i] =
          (g[j] * (t_start[i] - s[j - 1]) + g[j - 1] * (s[j] - t_start[i])) /
          ((double)(s[j] - s[j - 1]));
  }
  return;
}
/**************************************************************************/
/* generate integration weights */

static void genweights(long nq, double wid, double wgt[]) {
  double dq, qmin, qi, extra;
  long i;

  qmin = -wid;
  dq = 2. * wid / ((double)nq);

  for (i = 1; i < nq; i++) {
    qi = qmin + dq * ((double)i);
    wgt[i] = (qi + dq) * norm_accurate(qi + dq) - 2. * qi * norm_accurate(qi) +
             (qi - dq) * norm_accurate(qi - dq);
    wgt[i] = wgt[i] + (1. / sq2pi) * (exp(-(qi + dq) * (qi + dq) / 2.) -
                                      2. * exp(-qi * qi / 2.) +
                                      exp(-(qi - dq) * (qi - dq) / 2.));
  }

  wgt[0] = wgt[nq] =
      -qmin * norm_accurate(qmin) + (qmin + dq) * norm_accurate(qmin + dq) +
      (1. / sq2pi) *
          (exp(-(qmin + dq) * (qmin + dq) / 2.) - exp(-qmin * qmin / 2.));
  extra = qmin * norm_accurate(qmin) + 1. / sq2pi * exp(-qmin * qmin / 2.);

  if (pr == 0) {
    wgt[0] = wgt[0] + extra;
    wgt[1] = wgt[1] - extra;
  } else {
    wgt[nq] = wgt[nq] + extra;
    wgt[nq - 1] = wgt[nq - 1] - extra;
  }

  for (i = 0; i <= nq; i++)
    wgt[i] = wgt[i] / dq;

  return;
}

/**************************************************************************/
/* Calculate the value of the midatlantic's payoff at exercise date t_ex[j]
   given the value of phi                                                */

static double EXPAY(double phi, long npay, long i0, double pvpay[], double redj,
                    double strike, double G_pay[], double Gst, double zetaj) {
  double sum, ans;
  long i;

  sum = 0.;
  for (i = i0; i <= npay - 1; i++) {
    sum =
        sum + (pvpay[i] - redj) * exp(G_pay[i] * (phi - .5 * G_pay[i] * zetaj));
    redj = 0.;
  }

  ans = sum - strike * exp(Gst * (phi - .5 * Gst * zetaj));
  if (pr != 0)
    ans = -ans;

  /* OVE fix: prevent negative interpolated prices */
  if (ans < 0.0)
    ans = 0.0;

  return (ans);
}

/**************************************************************************/
/* Find sigma and tau values equivalent to g[j] and zeta[j].
   This routine calculates g(t) and zeta(t) at all event
   dates of the mid-atlantic  , determines the sigmas and taus that
   yield these values  , and stores them in the market (m) term
   structure.    */

static String Findsigtau(long n, Date tnow, long spot_lag)

{
  double lambda, target, glc;
  double dt, th;
  long j, k;
  int check;

  error = NULL;
  check = 0;

  /* set dates */
  num_sig = num_tau = n + 1;
  for (j = 0; j <= n; j++) {
    tau_date[j] = s[j];
    sig_date[j] = add_unit(s[j], -spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
  }
  glc = exp(-sy[n] / overall_tau);

  /* step backwards  , computing tau_val[k] from the integral g[j] */
  for (k = n; k > 0; k--) {
    dt = sy[k] - sy[k - 1];

    target = (g[k - 1] - g[k]) / (dt * glc);
    if (target <= Gnew(log(.0000001 / glc))) {
      lambda = (log(.0000001 / glc)) / dt;
      glc = .0000001;
    } else {
      if (target <= 1.0)
        th = 0.; /* initial guess */
      else
        th = log(target) + log(1. + log(target)) + log(2.);

      ScNewt(&th, &check, target, Gnew, Gnew_der);
      if (check != 0) {
        error = "can't convert tau";
        return (error);
      }

      lambda = th / dt; /* correct lambda */

      glc = glc * exp(th); /* update glc     */
    }

    if (fabs(lambda) < 0.005) /* convert to tau = 1/lambda */
      tau_val[k] = 200.0;
    else
      tau_val[k] = 1.0 / lambda;
  }

  /* get value at k=0 */
  tau_val[0] = tau_val[1];

  /* Correct the maturities by the spot_lag for the proper sigma computation */
  for (j = 0; j <= n; j++) {
    sy[j] -= (s[j] - add_unit(s[j], spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING)) /
             365.0;
  }

  /* Calculate sigma values to match zeta curve */
  /* note that glc is value at s[0] */
  th = sy[0] / tau_val[0];
  sig_val[0] = glc * sqrt(2. * zeta[0] / tau_val[0] / (1.0 - exp(-2. * th)));

  /* construct sigma at rest of the dates */
  for (k = 1; k <= n; k++) {
    th = (sy[k] - sy[k - 1]) / tau_val[k];

    sig_val[k] = glc * sqrt(2. / tau_val[k] * (zeta[k] - zeta[k - 1]) /
                            (exp(2. * th) - 1.));

    glc = glc * exp(-th);
  }

  return (error);
}
/*Wasn't that fun*/

/***********************************************************************/
/* Find the swap rates r_exer[j] at the exercise bocrvary for vanilla swaptions
   with exercise dates t_ex[j] and ending at tend */

static String Findexerbdry(Date tnow, long nex, Date t_ex[], Date tend,
                           double phi_cr[], long n, String yc_name, int icmpnd,
                           BasisCode natural_basis, BusDayConv natural_conv,
                           int natural_spot_lag,
                           SrtLgmExerBdryData *lgmExerBdryData) {
  long i, j, ncoup, mult;
  Date tfirstst, coupondate, trealend;
  double Bend, aB, sum;
  Date *tnatst = NULL, *paydate = NULL;
  double *zeta_ex = NULL, *G_pay = NULL, *G_st = NULL;
  double *B_st = NULL, *aB_pay = NULL, *aB_first = NULL;
  long *ifirstpay = NULL;

  /* This holds the start date & par rate for the exercise boundary swaptions */
  SrtLgmExerBdry *exerBdryArr = NULL;

  error = NULL;

  /* find number of coupon dates */
  if (natural_conv == NO_BUSDAY_CONVENTION)
    tfirstst = add_unit(t_ex[0], natural_spot_lag, SRT_DAY, natural_conv);
  else
    tfirstst = add_unit(t_ex[0], natural_spot_lag, SRT_BDAY, natural_conv);

  trealend = bus_date_method(tend, natural_conv);
  mult = 6;
  if (icmpnd == 0)
    mult = 12;

  ncoup = 0;
  coupondate = bus_date_method(tend, natural_conv);
  while (coupondate > tfirstst) {
    ncoup++;
    coupondate =
        add_unit(tend, -((int)(mult * ncoup)), SRT_MONTH, natural_conv);
  }

  /* allocate space */
  tnatst = (Date *)srt_calloc(
      nex + 1, sizeof(Date)); /* tnatst[0  , 1  , ...  , nex]     */
  exerBdryArr = srt_calloc(nex + 1, sizeof(SrtLgmExerBdry));
  zeta_ex = (double *)srt_calloc(
      nex + 1, sizeof(double)); /* zeta_ex[0  , 1  , ...  , nex]    */
  G_st = (double *)srt_calloc(
      nex + 1, sizeof(double)); /* G_st[0  , 1  , ...  , nex]       */
  B_st = (double *)srt_calloc(
      nex + 1, sizeof(double)); /* B_st[0  , 1  , ...  , nex]       */
  aB_first = (double *)srt_calloc(
      nex + 1, sizeof(double)); /* aB_first[0  , 1  , ...  , nex]   */
  ifirstpay = (long *)srt_calloc(
      nex + 1, sizeof(long)); /* ifirstpay[0  , 1  , ...  , nex]  */

  paydate = (Date *)srt_calloc(
      ncoup + 1, sizeof(Date)); /* paydate[0  , 1  , ...  , ncoup]  */
  G_pay = (double *)srt_calloc(
      ncoup + 1, sizeof(double)); /* G_pay[0  , 1  , ...  , ncoup]    */
  aB_pay = (double *)srt_calloc(
      ncoup + 1, sizeof(double)); /* aB_pay[0  , 1  , ...  , ncoup]   */

  if (tnatst == NULL || exerBdryArr == NULL || paydate == NULL ||
      zeta_ex == NULL || G_pay == NULL || G_st == NULL || G_st == NULL ||
      aB_pay == NULL || aB_first == NULL || ifirstpay == NULL) {
    error = "alloc. failed in Findexerbdry";
    return (error);
  }

  /* generate start_upon_exercise dates */

  if (natural_conv == NO_BUSDAY_CONVENTION)
    for (j = 0; j <= nex; j++)
      tnatst[j] = add_unit(t_ex[j], natural_spot_lag, SRT_DAY, natural_conv);
  else
    for (j = 0; j <= nex; j++)
      tnatst[j] = add_unit(t_ex[j], natural_spot_lag, SRT_BDAY, natural_conv);

  /* trim unused exercise dates */
  while ((nex >= 0) && (tnatst[nex] >= trealend))
    nex--;
  if (nex < 0) {
    error = "no exercise dates";
    ERROR_CHECK2;
  }

  lgmExerBdryData->NexerBdry = nex + 1;

  /* generate coupon dates (end & pay dates) */
  for (i = 1; i <= ncoup; i++)
    paydate[i] =
        add_unit(tend, -((int)(mult * (ncoup - i))), SRT_MONTH, natural_conv);
  paydate[0] = tfirstst;

  /* find first coupon date after each exercise date */
  i = 1;
  for (j = 0; j <= nex; j++) {
    while (paydate[i] <= tnatst[j])
      i++;
    ifirstpay[j] = i;
  }

  /* generate discount factors and discounted payments */
  for (i = 1; i <= ncoup; i++) {
    aB_pay[i] = swp_f_df(tnow, paydate[i], yc_name) *
                coverage(paydate[i - 1], paydate[i], natural_basis);
    if (aB_pay[i] <= 0) {
      error = "no disc. factor";
      ERROR_CHECK2;
    }
  }

  for (j = 0; j <= nex; j++) {
    B_st[j] = swp_f_df(tnow, tnatst[j], yc_name);
    aB_first[j] = swp_f_df(tnow, paydate[ifirstpay[j]], yc_name) *
                  coverage(tnatst[j], paydate[ifirstpay[j]], natural_basis);
  }
  Bend = swp_f_df(tnow, paydate[ncoup], yc_name);

  /* interpolate to get zeta and G at the needed dates */

  interpolate_par(tnow, n, nex, t_ex, zeta_ex, tnatst, G_st, ncoup + 1, paydate,
                  G_pay);

  /* compute par swap rate at the exercise boundary */
  for (j = 0; j <= nex; j++) {
    sum = 0.;
    for (i = ifirstpay[j]; i <= ncoup; i++) {
      if (i == ifirstpay[j])
        aB = aB_first[j];
      else
        aB = aB_pay[i];
      sum = sum + aB * exp(G_pay[i] * (phi_cr[j] - .5 * G_pay[i] * zeta_ex[j]));
    }

    exerBdryArr[j].parRate =
        (B_st[j] * exp(G_st[j] * (phi_cr[j] - .5 * G_st[j] * zeta_ex[j])) -
         Bend *
             exp(G_pay[ncoup] * (phi_cr[j] - .5 * G_pay[ncoup] * zeta_ex[j]))) /
        sum;

    exerBdryArr[j].bgnDate = t_ex[j];
  }

  /* output */
  lgmExerBdryData->exerBdryArr = exerBdryArr;

  /* free unneeded arrays & return */
  FreeLots(tnatst, paydate, zeta_ex, G_pay, G_st, B_st, aB_pay, aB_first,
           ifirstpay);
  return (error);
}

/***********************************************************************/
/* calculate Gnew(th) and its derivative */

static double Gnew(double th) {
  double ans;
  /* protect against th=0 */
  if (fabs(th) < 1.e-7)
    ans = 1 + .5 * th + th * th / 6.;
  else
    ans = (exp(th) - 1.) / th;
  return (ans);
}

static double Gnew_der(double th) {
  double ans;
  /* protect against th=0 */
  if (fabs(th) < 1.e-7)
    ans = .5 + th / 3.;
  else
    ans = (1. + (th - 1.) * exp(th)) / th / th;
  return (ans);
}

/************************************************************************/
/* This routine calculates the phi (called x) where EXPAY is zero */

static void init_point(double *x, int *check, double target, long i0, long npay,
                       double pvpay[], double redj, double G_pay[], double Gst,
                       double zetaj) {
  double sum, sum1, temp, xold, paymnt;
  long i, j;

  *check = 0;
  xold = 1.e30;
  *x = 0.;

  for (j = 0; j <= 10; j++) /* Newton step */
  {
    sum = sum1 = 0.;
    for (i = i0; i < npay; i++) {
      if (i == i0)
        paymnt = pvpay[i] - redj;
      else
        paymnt = pvpay[i];
      temp = paymnt *
             exp(-(Gst - G_pay[i]) * (*x - (Gst + G_pay[i]) * zetaj / 2.));
      sum = sum + temp;
      sum1 = sum1 - temp * (Gst - G_pay[i]);
    }
    *x = *x + (target - sum) / sum1;
    if (fabs(*x - xold) < 1.e-09)
      return;
    if (fabs(target / sum - 1.) < 1.e-09)
      return;
    xold = *x;
  }

  *check = 1; /* no convergence */
  return;
}

/************************************************************************/
/* Newton: This version of Newton is non-public since it is only
guaranteed to converge if f'(x) has one sign and f"(x) has one sig_valn */

static void ScNewt(double *x, int *check, double target, double (*func1)(),
                   double (*deriv1)()) {
  double value, crit, xold, xx, df;
  long j;

  *check = 0;
  xold = 1.e30;

  for (j = 0; j <= 10; j++) {
    xx = *x;
    value = (*func1)(xx);
    df = (*deriv1)(xx);
    *x = xx = xx + (target - value) / df; /* take Newton step */
    crit = .2e-12 * (1. + fabs(xx));

    if (fabs(xx - xold) < crit) /* check convergence */
      return;
    crit = .2e-14 * (1. + fabs(value));
    if (fabs(target - value) < crit)
      return;

    xold = xx; /* no convergence yet; take new Newton step */
  }

  *check = 1; /* no convergence*/
  return;
}

/* Newton: This version of Newton is non-public since it is only
guaranteed to converge if f'(x) has one sign.
   Note: x needs to be positive                                  */

static void ScSlowNewt(double *x, double xmin, double xmax, int *check,
                       double target, int gzfind, double (*func1)(),
                       double (*deriv1)()) {
  double xold, xnew, fold, fnew, crit1, crit2, dfdx, lambda;
  long j, k;

  /* initiailize  */
  *check = 0;
  xmin = xmin / 1.5;
  xmax = xmax * 1.5;
  xnew = *x;
  fnew = (*func1)(xnew, gzfind) - target;

  /* newton steps */
  for (j = 1; j <= 15; j++) {
    xold = xnew;
    fold = fnew;
    dfdx = (*deriv1)(xold, gzfind);
    if (fabs(dfdx) < 1.0e-14 * xold)
      /* OVE it has to be >= 0.0 otherwise blows up on out of the money (go the
       * wrong way) */
      dfdx = xold * 1.0e-14 * (dfdx >= 0 ? 1.0 : -1.0);
    xnew = xold - fold / dfdx; /* full newton step */

    /* reduce step if out of bounds */
    if (xnew < xmin)
      xnew = xmin;
    if (xnew > xmax)
      xnew = xmax;

    /* new function value */
    fnew = (*func1)(xnew, gzfind) - target;

    /* if needed  , reduce Newton step until error decreases */
    for (k = 1; k <= 5 && (fnew * fnew - fold * fold) > DBL_EPSILON; k++) {
      lambda = fold / (fold - fnew);
      if (lambda < .1)
        lambda = .1;
      xnew = xold + lambda * (xnew - xold);
      fnew = (*func1)(xnew, gzfind) - target;
    }

    /* check convergence */
    crit1 = .2e-12 * (1. + xnew);
    if (fabs(xnew - xold) < crit1) {
      *x = xnew;
      return;
    }
    crit2 = .2e-14 * (1. + fabs(target));
    if (fabs(fnew) < crit2) {
      *x = xnew;
      return;
    }

    /* update bounds */
    if (fold * fnew < 0.) {
      xmin = (xold < xnew ? xold : xnew);
      xmax = (xold > xnew ? xold : xnew);
    } else {
      if (xnew < xold)
        xmax = xnew;
      else
        xmin = xnew;
    }
  } /* end loop ... go back and take new Newton step */

  /* fell through loop ... check for convergence using looser bounds */
  crit1 = .5e-7;
  crit2 = .5e-7;
  if (fabs(xnew - xold) < crit1) {
    *x = xnew;
    return;
  }
  if (fabs(fnew) < crit2) {
    *x = xnew;
    return;
  }

  /* admit defeat */
  *check = 1;
  return;
}

/* ================================================================================
 */

/** OVE: gets the market convention:
                - either from the YC passed (if initialised)
                - or     from the default CCY parameters
**/

static Err get_ccy_defaults(SrtCurvePtr yldcrv, int fix_tau_flag, int usecaps,
                            BasisCode *basis, SrtCompounding *compnd,
                            BusDayConv *conv, int *spot_lag) {
  String ccy_str;
  SrtCcyParam *ccy_param;
  Err err;

  /* OVE: this will work as long as SrtInitYc is called  , that sets default
   * ccy_param */

  ccy_param = get_ccyparam_from_yldcrv(yldcrv);
  if (!ccy_param) {
    ccy_str = get_curve_ccy(yldcrv);
    err = swp_f_get_CcyParam_from_CcyStr(ccy_str, &ccy_param);
    if (err)
      return err;
  }

  if ((usecaps != 1) || (fix_tau_flag != 1))
  /* calibration with swaptions */
  {
    *basis = ccy_param->swap_basis_code;
    *compnd = ccy_param->compd;
    *conv = ccy_param->swap_bus_day_conv;
  } else
  /* calibration with caplets */
  {
    *basis = ccy_param->cash_basis_code;
    *compnd = SRT_QUARTERLY;
    *conv = ccy_param->cash_bus_day_conv;
  }

  *spot_lag = ccy_param->spot_lag;

  return NULL;
}

/***************************************************************************************/
static void FreeMost(void *a1, void *a2, void *a3, void *a4, void *a5,
                     void *a6) {
  srt_free(a1);
  srt_free(a2);
  srt_free(a3);
  srt_free(a4);
  srt_free(a5);
  srt_free(a6);

  return;
}

static void FreeMore(void *a1, void *a2, void *a3, void *a4, void *a5, void *a6,
                     void *a7, void *a8, void *a9, void *a10) {
  srt_free(a1);
  srt_free(a2);
  srt_free(a3);
  srt_free(a4);
  srt_free(a5);
  srt_free(a6);
  srt_free(a7);
  srt_free(a8);
  srt_free(a9);
  srt_free(a10);
  return;
}

static void FreeLots(void *a1, void *a2, void *a3, void *a4, void *a5, void *a6,
                     void *a7, void *a8, void *a9) {
  srt_free(a1);
  srt_free(a2);
  srt_free(a3);
  srt_free(a4);
  srt_free(a5);
  srt_free(a6);
  srt_free(a7);
  srt_free(a8);
  srt_free(a9);
  return;
}

static void FreeRest(void *a1, void *a2) {
  srt_free(a1);
  srt_free(a2);
  return;
}

/**********************************************************************/

#undef ERROR_CHECK
#undef ERROR_CHECK1
#undef MAX_NUM_DATES
#undef MIN_EX_INT
#undef MAX_REF_DATES
#undef FAR_OUT
#undef WIDTH
#undef NDIV
#undef NQMAX
#undef STENCIL
