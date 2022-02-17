#ifndef SWP_H_CMSOPT_H
#define SWP_H_CMSOPT_H

#define CMS_TOLERANCE 1e-10
#define CMS_EPSILON 1.0e-8
#define DEFAULT_CMS_DELTA_STRIKE 0.0005
#define DEFAULT_CMS_DELTA 0.0005
#define MAX_CMS_SWAPS 8000
#define MAX_STDEV 8

/* This the New version including Smile and swaption replication for GRFN Use
   only
        -> Flat Vol*/
Err swp_f_cmsoption(double, double, SrtCompounding, double, double, double,
                    SrtReceiverType, double, double, int, SrtDiffusionType,
                    double *);

Err swp_f_Cms_Option(
    double dFwdSwapRate,                /*Forward Swap Rate */
    double dMaturity,                   /* CMS Maturity as double */
    double dNumPeriods, double dStrike, /* Option Strike */
    double dFrequency,                  /* CMS Frequency */
    SrtReceiverType PayRec,             /* Pay or Rec */
    double dDelay,                      /* Delay */
    double dRateConv,                   /* date adjustment */
    SrtDiffusionType VolType,           /* Vol type Lognormal or Normal */
    double dFlatVol,                    /* Flat Vol if used */
    int iMethod, /* 0: Use Flat Vol      , 1: linear interpolation      , 2:
                    FullSmile*/
    Date dStart, /* The following parameters are used in the GetVol function */
    Date dEnd,   /* and are useless in the rest of the code */
    SRT_Boolean bAdjForSpread, double dSpread, char *szVolCurveName,
    long lNumStrikesInVol, double *pdStrikesVol, double *dCMSOptionValue);

/* Calculates the CMS option taking into account the fact that in NY      ,
market data are physical-delivery swaptions and not cash-settle swaptions */
Err swp_f_Cms_OptionNY(
    double dFwdSwapRate,                /*Forward Swap Rate */
    double dMaturity,                   /* CMS Maturity as double */
    double dNumPeriods, double dStrike, /* Option Strike */
    double dFrequency,                  /* CMS Frequency */
    SrtReceiverType PayRec,             /* Pay or Rec */
    double dDelay,                      /* Delay */
    double dRateConv,                   /* date adjustment */
    SrtDiffusionType VolType,           /* Vol type Lognormal or Normal */
    double dFlatVol,                    /* Flat Vol if used */
    int iMethod, /* 0: Use Flat Vol      , 1: linear interpolation      , 2:
                    FullSmile*/
    Date dStart, /* The following parameters are used in the GetVol function */
    Date dEnd,   /* and are useless in the rest of the code */
    SRT_Boolean bAdjForSpread, double dSpread, char *szVolCurveName,
    long lNumStrikesInVol, double *pdStrikesVol, double *dCMSOptionValue,
    String ycName, String refRateCode);

Err swp_f_Tec_Option(
    double dFwdSwapRate,                /*Forward Swap Rate */
    double dMaturity,                   /* TEC Maturity as double */
    double dNumPeriods, double dStrike, /* Option Strike */
    double dMargin,                     /* the Tec Margin */
    double dFrequency,                  /* TEC Frequency */
    double dPaymentPower,               /* TEC Power */
    SrtReceiverType PayRec,             /* Pay or Rec */
    double dDelay,                      /* Delay */
    double dRateConv,                   /* date adjustment */
    SrtDiffusionType VolType,           /* Vol type Lognormal or Normal */
    double dFlatVol,                    /* Flat Vol if used */
    int iMethod, /* 0: Use Flat Vol      , 1: linear interpolation      , 2:
                    FullSmile*/
    Date dStart, /* The following parameters are used in the GetVol function */
    Date dEnd,   /* and are useless in the rest of the code */
    SRT_Boolean bAdjForSpread, double dSpread, char *szVolCurveName,
    long lNumStrikesInVol, double *pdStrikesVol, double *dTECOptionValue);

/* Flat version using direct integration of the Payoff */
Err swp_f_FlatTecCmsOption(
    double dFwdSwapRate, Date StartDate, Date TheoEndDate, double dMaturity,
    double dNumPeriods, double dStrike, double dFrequency,
    SrtReceiverType PayRec,   /* Pay or Rec */
    double dDelay,            /* Delay */
    SrtDiffusionType VolType, /* Vol type Lognormal or Normal */
    double dFlatVol,          /* Flat Vol if used */
    double dRateConv,         /* Basis adjustment */
    double dMargin,           /* Margin for Tec */
    double dPower,            /* Power for Tec */
    double *dTecCmsGenericValue);

Err swp_f_BnpCmsRate(double dFwdSwapRate, /*Forward Swap Rate */
                     double dMaturity, double dNumPeriods, double dFrequency,
                     double dDelay, double dVol, SrtDiffusionType VolType,
                     int iMethod, double dlambda, double dCmsSpread,
                     double *dBnpCmsRateValue);

void Implied_SpreadLambda(double dFwdSwapRate,  /*Forward Swap Rate */
                          double dCmsFullSmile, /*Forward Cms rate */
                          double dCmsVegaFS,    /* Vega full smile ???*/
                          double dMaturity, double dNumPeriods,
                          double dFrequency,

                          double dVol, SrtDiffusionType VolType,
                          double *dlambda, double *dCmsSpread);

/* TEST FUNCTIONS */

void set_numberpoints(double numpoints, int type);
void set_upperboundvol(double strike);
void set_upperboundIntegral(double bound);
double get_upperbound(void);
double get_upperstrike(void);
double get_uppervol(void);

Err swp_f_truncvol(char *vol_id, Ddate start, Ddate end, double strike,
                   double *volatility, double *power);

typedef struct {
  double vol;  /* Vol or BetaVol */
  double beta; /* SABR or HESTON params */
  double alpha;
  double rho;
  double zeta;
  double lambda;     /* mean reversion of the vol */
  double volinfty;   /* mean reversion level of the vol */
  double NstD;       /* mean reversion of the vol */
  double Pi;         /* mean reversion level of the vol */
  double bump_vol;   /* for smile method */
  int bump_type;     /* 0: additive      , 1: mult */
  int is_vol_bumped; /* 0 no : 1 yes */
  int atm;           /* ATM or strike */
  int input_method;  /* 0: Use Flat Vol      , 1: linear interpolation      ,
                                2: FullSmile      ,    3: ??      , 4: Atm or
                        strike*/
  int comp_method;
  int input_SABR;
  int input_SABRAF;
  char *szVolCurveName; /* getvol params */
  SrtDiffusionType voltype_input;
  SrtDiffusionType voltype_used;
  Date start;
  Date end;
} CMSVol;

typedef struct {
  double CAP_VolUD;
  double CAP_IntUD;
  double FLOOR_VolUD;
  double FLOOR_IntUD;
  long Nx;
  int CAP_Mesh;
  int FLOOR_Mesh;
} CMSParams;

typedef struct {
  double vol_bump;
  int vol_bump_type; /* 0: add      , 1: mult */
  double alpha_bump;
  double beta_bump;
  double rho_bump;
  double yc_bump;
  long num_bump;
  int greek_type; /* 0 : none      , 1: all      , 2: vol      ,
                     3: alpha      , 4: beta      , 5: rho */
} CMSGreeks;

void CMSTECfreevol(CMSVol *Volstruct);
void CMSTECinitvol(CMSVol *Volstruct);
void CMSTECinitParams(CMSParams *Paramstruct);

/* if one wants to use it
/*
/* First step  : init the struct : CMSTECinitvol(&myVol);
/* Second step : define the method + (optional params szVolCurveName      ,
numstrikes      , strikes      , vol      , beta      , etc....if necessary
/* Third step  : call the function
/* Last step   : free the structure :  CMSTECfreevol(&myVol);
/*
/* Example of second step :
/* 1. flat log vol case	: myVol.method = 0      , myVol.vol = 10%      ,
myVol.voltype = SRT_LOGNORMAL */
/* 2. approx smile log	: myVol.method = 1      , myVol.voltype = SRT_LOGNORMAL
, + fill myVol.start      , myVol.end      , myVol.szVolCurveName      ,
myVol.numstrikes      , myVol.strikes
/* 3. full sabr			: myVol.method = 2      , myVol.szVolCurveName =
vol_name(mkt)      , myVol.start      , myVol.end       , myVol.voltype =
Inputtype atm vols (norm      , log      , beta);
/* 4. strike vol		: myVol.method = 4; myVol.atm = 0;
myVol.szVolCurveName = vol_name(mkt); myVol.voltype = type of vol in the
mkt(norm      , log);
*/

Err swp_f_CMSTEC_Option(                     /* underlying */
                        double dFwdSwapRate, /*Forward Swap Rate */
                        double dnfp,  /* number of full period of the default
                                     swap rate */
                        double dfreq, /* frequency of the default swap rate 1:A
                                         , 2:S      , 4:Q */
                        double dRateConv, /* adjustment swap default/swap und */
                        /* Option parameters */
                        double dMaturity,       /* CMS Maturity as double */
                        double dStrike,         /* Option Strike */
                        SrtReceiverType PayRec, /* Pay or Rec */
                        double dPower, /* 1: CMS      , 4: tec      , etc.... */
                        double dMargin,
                        /* delay parameters */
                        double dCvgZC, /* cvg of the delay */
                        double dFwdZC, /* B(0      ,Tpay) / B(0      ,Tstart) */
                        /* vol structure */
                        CMSVol *Volstruct, /* Flat Vol*/
                        /* Greek struct */
                        CMSGreeks *Greekstruct,
                        /* numerical param */
                        CMSParams *params,
                        /* Output */
                        double *dCMSOptionValue);

Err swp_f_CMSTEC_rate(                     /* underlying */
                      double dFwdSwapRate, /*Forward Swap Rate */
                      double dnfp,  /* number of full period of the default swap
                                   rate */
                      double dfreq, /* frequency of the default swap rate 1:A ,
                                   2:S      , 4:Q */
                      double dRateConv, /* adjustment swap default/swap und */
                      /* Option parameters */
                      double dMaturity, /* CMS Maturity as double */
                      double dPower, /* 1: CMS      , 4: tec      , etc.... */
                      double dMargin,
                      /* delay parameters */
                      double dCvgZC, /* cvg of the delay */
                      double dFwdZC, /* B(0      ,Tpay) / B(0      ,Tstart) */
                      /* vol structure */
                      CMSVol *Volstruct, /* Flat Vol*/
                      /* Greek struct */
                      CMSGreeks *Greekstruct,
                      /* numerical param */
                      CMSParams *params,
                      /* Output */
                      double *dCMSTECrate);

Err swp_f_CMSTEC(                     /* underlying */
                 double dFwdSwapRate, /*Forward Swap Rate */
                 double
                     dnfp, /* number of full period of the default swap rate */
                 double dfreq, /* frequency of the default swap rate 1:A      ,
                                  2:S , 4:Q */
                 double dRateConv, /* adjustment swap default/swap und */
                 /* Option parameters */
                 double dMaturity,       /* CMS Maturity as double */
                 double dStrike,         /* Option Strike */
                 SrtReceiverType PayRec, /* Pay or Rec */
                 double dPower, /* 1: CMS      , 4: tec      , etc.... */
                 double dMargin,
                 /* delay parameters */
                 double dCvgZC, /* cvg of the delay */
                 double dFwdZC, /* B(0      ,Tpay) / B(0      ,Tstart) */
                 /* vol structure */
                 CMSVol *Volstruct, /* Flat Vol*/
                 /* numerical param */
                 CMSParams *params,
                 /* Output */
                 double *dCMSTEC);

Err swp_f_CMSTEC_GREEKS(                     /* underlying */
                        double dFwdSwapRate, /*Forward Swap Rate */
                        double dnfp,  /* number of full period of the default
                                     swap rate */
                        double dfreq, /* frequency of the default swap rate 1:A
                                         , 2:S      , 4:Q */
                        double dRateConv, /* adjustment swap default/swap und */
                        /* Option parameters */
                        double dMaturity,       /* CMS Maturity as double */
                        double dStrike,         /* Option Strike */
                        SrtReceiverType PayRec, /* Pay or Rec */
                        double dPower, /* 1: CMS      , 4: tec      , etc.... */
                        double dMargin,
                        /* delay parameters */
                        double dCvgZC, /* cvg of the delay */
                        double dFwdZC, /* B(0      ,Tpay) / B(0      ,Tstart) */
                        /* vol structure */
                        CMSVol *Volstruct, /* Flat Vol*/
                        /* Greek struct */
                        CMSGreeks *Greekstruct,
                        /* numerical param */
                        CMSParams *params,
                        /* Output */
                        double *dCMSTEC);

/* Get the Weights and the points for the intregration */
/* The numeraire is not taken into account */
Err swp_f_CMSTEC_Weights(              /* underlying */
                         double dnfp,  /* number of full period of the default
                                      swap rate */
                         double dfreq, /* frequency of the default swap rate 1:A
                                          , 2:S      , 4:Q */
                         double
                             dRateConv, /* adjustment swap default/swap und */
                         /* Option parameters */
                         double dStrike,         /* Option Strike */
                         SrtReceiverType PayRec, /* Pay or Rec */
                         double
                             dPower, /* 1: CMS      , 4: tec      , etc.... */
                         double dMargin,
                         /* delay parameters */
                         double dCvgZC, /* cvg of the delay */
                         double
                             dSpdZC, /* B(0      ,Tpay) / B(0      ,Tstart) */
                         /* numerical param */
                         CMSParams *params,
                         /* Output */
                         double *pdW, double *pdX);

#endif
