#ifndef TBAUTIL_H
#define TBAUTIL_H

#include "MBSPPFuncs.h"
#include "MBSPTUtil.h"
#include "stdio.h"
/***macros: constants****/
// switches to decide whether to use certain corrections:

/////

#define NbSFS 200

#define MBS_ERROR 1E-7 /* Error tolerance */
#define MAXMAT 30      /* Maximum maturity for bechmark swaps */
#define JUMPCOEFF 1.5  /* Coefficient used to determine the size of the jump */
#define JUMPCOEFF2                                                             \
  2. /* Coefficient used to determine the size of the second jump */
#define NBEXER 50 /* Maximum number of call dates */
#define NBFWD 50  /* Maximum number of forward prices */
#define NBSWAPTION                                                             \
  5 /* Number of benchmark swaptions used to compute spot volatility */

#define MAXBUFF 128
#define FEBEFFECT 0.93442623
// directories
#define DATADIR "O:/Options/Tba/datfiles/"
#define PPFDATADIR "O:/Options/Tba/Ppf/datfiles/"
#define OUTPUTDIR "O:/Options/Tba/outfiles/"

/* Maximun number of period-dependent prepayments */
#define ACURA 0.03
/* Acuracy of most estimation of implied ppf or vol etc. */
#define NTRIALS 12
/* max number of iterations to get implied stuff */
#define Nsigma 4
//#define THRESHOLD 0.9

/*dates related:*/

#define OKAY 0
#define ERR -1

/*macro function*/

/*** Structs***/
typedef struct {
  long Today,                  /* Today's date */
      Spotdays;                /* Number of business days until the spot date */
  int MMB,                     /* Money market basis (360 or 365) */
      IndexMat,                /* Index maturity in years */
      SwaptionExp[NBSWAPTION], /* Expiration (in years) of the benchmark
                              swaptions */
      TStree,                  /* 0 for treasury tree      , 1 for swap tree */
      TSdiscount;    /* 0 for treasury discount      , 1 for swap discount */
  char AoS,          /* Yield curve type (annual or semi-annual) */
      YearBase[4],   /* Year basis convention for the swap curve */
      IndexFreq;     /* Index frequency */
  double MMYield[6], /* Money market rates (O/N      ,1m      ,2m      ,3m ,6m
                    ,1yr) for tree building */
      MSYield[6],   /* Libor rates if above is treasury      , vice versa      ,
                   for   discounting. */
      SwapYield[8], /* Swap rates (2      ,3      ,4      ,5      ,7      ,10
                       ,20      ,30yr) */
      YieldS[8],    /* Yield rates for discount */
      Beta,         /* Mean reversion coefficient */
      SwaptionVol[NBSWAPTION]; /* Volatility of the benchmark swaptions */
  long ValueDate,              /* Output: value date */
      Period[2 * MAXMAT + 4], /* Output: length of money market periods followed
                             by 6m increments */
      Holidays[MAXHOLIDAYS];
  int numHolidays;
  double OAS,               /* Option adjusted spread (cf Lattice() in dev2.c */
      Zero[2 * MAXMAT + 4], /* Output: Zeros at money market points followed by
                           6m interval points */
      ZeroS[2 * MAXMAT + 4]; /* Output: Zeros at mm points followed by 6m
                                interval points for YieldS*/

} TERM_DATA;

typedef struct {
  double *ZeroCoupon, /* Zero coupon price maturing at each node */
      *ZeroCouponS,   /* Zero coupon price of discounting bond at each node */
      *FwdRate,       /* Forward rate from one period to the following */
      *FwdRateS,  /* Discounting forward rates from one period to the following
                   */
      *FwdVol,    /* Instantaneous volatility of forward rate */
      *Drift,     /* Drift of the short term rate at each node */
      *Delta,     /* convexity adjustment for discounting rates when tree rates
                 differ */
      *FwdPrepay, /* Forward prepayment at each node in the tree */
      *PrepayVol, /* Instantaneous volatility of forward prepayment */
      *Length,    /* Length of time steps */
      *Accrued,   /* Accrued interest at each node */
      *FwdStrike, /* Strike value at each node */
      Jump,       /* Size of interest rate jump (in log space) */
      Jump2;      /* Size of the prepayment jump (in log space) */
  int Ppy,        /* Number of period per year */
      NbNodes,    /* Total number of nodes */
      ExEnd,      /* Last exercise date expressed in number of nodes */
      NodeMax,    /* Maximum node reached in the interest rate tree */
      NodeMax2,   /* Maximum node reached in the prepayment tree */
      ExerPosition[NBEXER], /* Position (node number) of exercise dates */
      FwdPosition[NBFWD],   /* Position (node number) of forward prepayments */
      *ExerFlag; /* 1 if node is an exercise date      , 0 otherwise */
} TREE_DATA;

typedef struct {
  double fast;          /* fast triggered speed */
  double slow;          /* slow triggered speed */
  double trigrs[NbSFS]; /* trigger REFINANCE    rates   */
  double slowpo[NbSFS]; /* slow prepayment price -- po*/
  double slowio[NbSFS]; /* slow prepayment price --io */
  double fastpo[NbSFS]; /* fast price -- po */
  double fastio[NbSFS]; /* fast prepay price--io*/
  int ngroup;           /* number of trigger rates in each group */
} TRIGGER;

typedef struct {
  double IndexLevel[NBPOINTS]; /* base ppf index level- refirate */
  double AmortLevel[NBPOINTS]; /* base ppf */
  double SeasoningLevel[3];    /* seasoning refirate level */
  int SeasoningMat[3];         /* seasoning months */
  double Seasonality[12];      /* seasonality factor */
  int Delay;                   /* delay in pp */
  double Speedup;              /* speed up factor */
  double rational;
} PREPAY;

typedef struct {
  int NbFwd;            /* Number of forward prepayments */
  long FwdMat[NBFWD];   /* Maturity of the forward prepayments */
  double Spot,          /* Spot value of the prepayment */
      FwdPrepay[NBFWD], /* Forward values of prepayment */
      SpotVol,          /* Spot volatility of prepayment */
      Beta;             /* Mean reversion of prepayment */
} PREPAY_DATA;

typedef struct {
  long ptime;      /* present time yyyymm */
  long otime;      /* origination time  */
  long wam;        /* wam */
  int term;        /* original term of mortgage in month */
  double Gwac;     /* Wac gross	*/
  double Refirate; /* today's refirate */
  int Delay;       /* delay of prepayment in number of months*/
  double Change;   /* forecast ppf with a change in refirate */
  char FoG;        /* Fanny or Ginnie */
} DEAL;

typedef struct {
  int SFSNb,     /* Number of SFS per group */
      Nbweights, /* number of weights used in prepayment function */
      NbDeals,   /* number of deals to be priced */
      Delay, /* Number of coupon delay between index fixing and prepayment */
      DelayPmt,
      /* Number of days of payment delay */
      NbFixedPrepay,        /* Number of fixed (deterministic) prepayments */
      NbExer,               /* Number of exercise dates */
      Term;                 /* Term of the mortgage in months */
  long Maturity,            /* Maturity of the deal */
      Exer[NBEXER];         /* Exercise dates */
  double Strike[NBEXER],    /* Strikes */
      Coupon,               /* Underlying coupon (annualized) */
      Gwac,                 /* Gross Coupn (service fee included */
      Refirate,             /* Refinance rate of the trading day */
      IndexLevel[NBPOINTS], /* Index level that define the prepayment function
                             */
      AmortLevel[MAXNUMP]
                [NBPOINTS], /* Corresponding annual CPR level or Avg Life*/
      AmortLevel1[MAXNUMP][NBPOINTS],
      Amort[MAXSFSNB], /* Additional amortization once the trigger is reached
                      for each SFS */
      Speedup,         /* Speed up the PPF */
      rational,        /* rationality speed up */
      TriggerRate[MAXSFSNB],    /* Trigger rate for each SFS */
      SeasoningLevel[NUMRAMPS], /* 3 Levels of Refi/Coupon that define seasoning
                               (cf SeasoningMat) */
      SeasoningMat[NUMRAMPS],   /* 3 Maturities that define seasoning */
      Seasonality[12],          /* Seasonality level for each month */
      InitialSMM[MAXFIXEDPREPAY],
      FixedPrepay[MAXFIXEDPREPAY]; /* Deterministic prepayment level */
  char Freq,                       /* Frequency of mortgage payments */
      Hedge, /* Calculate the hedge results ('Y' or 'N') */
      Ppform;
  char Scenario;  /* 'Y' means to do parallel shift of YC */
  double Shift;   /* shift amount */
  char Schedule;  /* 'N' means no scheduled amortization */
  char AlorCPR;   /* 'A' means average life      , 'C' means CPR */
  char Balloon;   /* 'Y' means ballon */
  int BalloonSch; /* original amortization schedule of    balloon  */
  char Always;    /* if yes      , always build tree 10 yr beyond */
  double Inpio;   /* input io */
  double Inppo;   /* input po */
  double stripratio;
  double collatoral;
  int casenumber;  /* I built in the program of 10 cases */
  double implvol;  /* implied volatility multiplier  */
  int pay_delay;   /* payment delay */
  char Swaption;   /* swaption: No      , Call      , Put */
  int ExerT;       /* exercise time: periods from now */
  char invfloater; /* inverse floater or not */
  double output[10];
} DEAL_DATA;
/***Functions***/
char *Manager(char *data_dir,
              TERM_DATA *term_data, /* Structure of term structure data  */
              TREE_DATA *tree_data, /* Structure of tree data */
              DEAL_DATA *deal_data, /* Structure of deal structure data  */
              PREPAY_DATA *prepay_data,
              HISTORY *history,     // historical index values
              MBSDENSITY *density); // density function for prepay propensity

char *
Input_Check(TERM_DATA *term_data,     /* Structure of term structure data  */
            DEAL_DATA *deal_data,     /* Structure of deal structure data  */
            PREPAY_DATA *prepay_data, /* Structure of prepayment data */
            TREE_DATA *tree_data,
            HISTORY *history,     // historical index values
            MBSDENSITY *density); // density function for prepay propensity

int Conv_Freq(char Freq);

double Conv_AoS(char AoS);

int Holiday(long date_i, long *holidays, int numHolidays);
long Nxtbusday(long date, int len, long *holidays, int num);

//#ifdef USESRTDATES

long Daysact(long date1_i, long date2_i);

//#endif//USESRTDATES

void pksort3(int n, double *arr, int *brr, int *crr);

/*interpolation and/or integration*/
double qinterpp(int price_n, double *price_x, double *price, int dens_nb,
                double *dens_x, double *dens);

double qinterpi(double *price_x, double *price, double *dens_x, double *dens);

double qinterps(double *x, double *y, double x1, double *coeff);

void linterp(double x, double *y, double x0, double x1, double y0, double y1);

double linterpFlat(int n, double *xs, double *ys,
                   double x); // interpolate or do flat extrapolation

void merror_register(char *err_text);

int *mbs_ivector(int nl, int nh);

double *mbs_dvector(int nl, int nh);

double **mbs_dmatrix(int nrl, int nrh, int ncl, int nch);

void mbs_free_ivector(int *v, int nl, int nh);

void mbs_free_dvector(double *v, int nl, int nh);

void mbs_free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

void Tree_Free(int NbNodes, TREE_DATA *tree_data);

void rs(double *io, double *po, double *rtn, double *spd, double I, double P,
        double *rr, double *ss, double *IO_spd, double *IO_rtn, double *PO_spd,
        double *PO_rtn);

double smmal1(double smm, int N);

//////////

#endif