
#ifndef _SRT_H_LGM_US_TYPES_H_
#define _SRT_H_LGM_US_TYPES_H_

#include "opsabr.h"
#include "srt_h_lgmtypes.h"

/**********************************************************************************/
/* Structure for Callable Time Swap */
typedef struct SrtCallTimeSwap_ {
  /* Exercise information */
  long nEx;     /* number of exercise dates */
  long FirstEx; /* Reserved for LGMautocal */
  Date *tEx;    /* [0      ,1      ,...      ,nEx-1] notification dates; deal
                   starts on next    coupon date */
  Date *tSet;   /* [0      ,1      ,...      ,nEx-1] settlement dates */
  long *iSet; /* [0      ,1      ,...      ,nEx-1] on exer      , settlement is
                 tCpnStart[iSet[j]] */
  double *strike; /* [0      ,1      ,...      ,nEx-1] fee paid by option holder
                     to exercise */
  SrtReceiverType PayRec; /* RECEIVER or PAYER */
                          /*	Coupons: each coupon = marg*cvg + (marg_accrued +
                           * gear*LiborCash)*accrued*cvg_accrued */
  long nCpn;              /* number of coupon periods */
  Date *tCpnStart; /* [0      ,1      ,...      ,nCpn-1] coupon start date */
  Date *tCpnPay;   /* [0      ,1      ,...      ,nCpn-1] coupon pay date */
  double *tCpn;
  double *tFundingPayment;
  double *gear; /* in general      , 0 (fixed coupon) or 1 (floating digital) -
                   [0 ,1      ,...      ,nCpn-1]  margins accrued*/
  double *tCvgCpn; /* [0      ,1      ,...      ,nCpn-1]  margins accrued  */
                   /* Barriers details */
  double call_spread_up;
  double call_spread_low;
  double **barriers; /* [i      ,0]: lower barrier      , [i      ,1]: upper
                        barrier */
                     /* Observation details */
  long observation_freq;
  Date ***observationdays; /* [0      ,.. ,nCpn-1]*[0..observation_freq-1][0..1]
                              , 0 for the Start      , 1 for the End  */
  double **ratiodays_for_subperiod; /* [0      ,..
                                       ,nCpn-1]*[0..observation_freq-1] */
                                    /* Index details */
  long undIndex;
  long undAccrued;
  double CorrelStart; /* for floating digital      , correlation between the
                         cash libor and the first fixing of the period; */
  double CorrelEnd; /* for floating digital      , correlation between the cash
                       libor and the last exercise of the period */
                    /* Time swap forward pricing */
  double **tForwardTS;   /* [0..1]*[0..nEx-1] for each exercise date. Filled in
                            LGMExtractInfoFromDeal  */
  int timeswapvolmethod; /*  2 = Sliding      , 3 = Converging  */
  int atmvolmethod; /* 1 = Lognormal      , 2 = Normal      , 3 = SigmaBeta */
  int calibrationmethod; /* 1 = Swaptions ATM      , 2 = Swaptions Strike eq */
  double dMatStd; /* 0 if no max      , otherwise vol is retrieved for strikes
                     bounded by Fe^{+/-dMaxStd Vol sqrt{t}} */
  long num_subdiscret_obsfreq;
  double shiftvol;
  long nx;
  int nofunding;
  char *vol_id;
  LGMErr (*GetSABRParam)(long, long, char *, SABR_VOL_TYPE, double *);
  /*	LGMErr (*GetATM)(long      , long      , double *);
          LGMErr (*GetAlpha)(long      , long      , double *);
          LGMErr (*GetBeta)(long      , long      , double *);
          LGMErr (*GetRho)(long      , long      , double *); */
  double pv_fixedleg;
  /* Record of fwd vol */
  long fwdVolCpn;
  long fwdVolEx;
  double **fwdVolMat;
  /*	double	*reval_times;	*/
} SrtCallTimeSwap, *SrtCallTimeSwapPtr;

#endif
