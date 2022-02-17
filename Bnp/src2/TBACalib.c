#include "math.h"
#include "stdio.h"

#include "MBSPTCalib.h"
#include "TBACalib.h"
#include "TBAUtil.h"

/*****  ExpDecay
 * ******************************************************************/
/*
 *       Exponential decay function: ExpDecay (a  ,t) = (1 - exp (at)) / (a*t)
 */
double ExpDecay(double a, double t) { return (VolFactor(a, t)); } /* ExpDecay */

/*****  Zero_6m  ************************************************************/
/*
 *       Calculate zeros and their corresponding maturity  , at 6 month
 * interval.
 */
void Zero_6m(TERM_DATA *term_data) /* Structure of term structure data */
{
  double *SwapYield, /* Array of swap yields at 6mo interval  , both benchmark
                        and interpolated ones */
      *YieldS;       /* Array of discount yields at 6mo interval */
  int IndexSwapYield[] = {4,  6,  8, 10, 14,
                          20, 40, 60}, /* Swap benchmark points in 1/2 years */
      i;

  SwapYield = mbs_dvector(0, 2 * MAXMAT);
  YieldS = mbs_dvector(0, 2 * MAXMAT);

  //	term_data->ValueDate = Nxtbusday (      term_data->Today  , /* Value
  //date is today shifted by (Spotdays) business days */ 						term_data->Spotdays  ,
  //term_data->Holidays  , term_data->numHolidays);

  for (i = 0; i < 8;
       i++) { /* We fill the swap yield array with the benchmark yields */
    SwapYield[IndexSwapYield[i]] = term_data->SwapYield[i];
    YieldS[IndexSwapYield[i]] = term_data->YieldS[i];
  }

  yldconv(term_data->Zero,  /* Conversion of the money market rates. */
          term_data->ZeroS, /* Conversion of MM rates for discounts. */
          term_data->Period, term_data->MMYield, term_data->AoS,
          term_data->ValueDate, term_data->MMB);

  coupint(SwapYield); /* Interpolation of Swap yields by 6mo increments.  */
  coupint(YieldS);    /* interpolation of D Yields by 6m increments. */
  spotrate(term_data->Zero +
               3, /* Bootstapping of swap curve to get zeros at 6mo intervals */
           term_data->ZeroS + 3,
           term_data->Period + 3, /* We pass the pointers Zero and Period
                                     starting at 3 so that (Zero + 3)[1] */
           SwapYield, /* = Zero[4] which is the 6m zero  , the first used for
                         bootstrapping.        */
           YieldS, term_data->AoS, term_data->ValueDate, term_data->YearBase);

  mbs_free_dvector(SwapYield, 0, 2 * MAXMAT);
  mbs_free_dvector(YieldS, 0, 2 * MAXMAT);

  return;

} /* Zero_6m */

/*****  yldconv  ************************************************************/
/*
 *       Convert Money Market rates into Bond Equivalent Yields (ACT/365).
 */
void yldconv(
    double *Zero,    /* output: zeros at money market points */
    double *ZeroS,   /* output: zeros at mm points */
    long *Period,    /* output: length of money market periods */
    double *MMYield, /* Money Market rates (O/N  , 1  ,2  ,3  ,6  ,12 m) */
    char AoS,        /* 'A'nnual or 'S'emi-Annual frequency */
    long ValueDate,  /* Value date for the benchmark instruments */
    int MMB)         /* Money Market basis (360 or 365) */
{
  long NewDate, /* Maturity (date) of the current Money Market rate */
      IndexMMYield[] = {0, 1, 2,
                        3, 6, 12}; /* Money Market benchmark points in months */
  int i;

  Zero[0] = pow(1. + MMYield[0] / 100. / (double)MMB, 365.) -
            1.; /* O/N Money Market rate */
  Period[0] = 1;
  ZeroS[0] = Zero[0];

  for (i = 1; i < 6; i++) /* 1  ,2  ,3  ,6  ,12 month Money Market rates */
  {
    NewDate = Nxtmth(ValueDate, IndexMMYield[i], (long)1);
    Period[i] = Daysact(ValueDate, NewDate);
    Zero[i] = pow(1. + MMYield[i] / 100. * Period[i] / (double)MMB,
                  365. / (double)Period[i]) -
              1.;
    ZeroS[i] = Zero[i]; /* for now ZeroS = Zero */
  }                     /* for i */

  if (AoS == 'S') /* Convert Annual rates to Semi-Annual if needed */
    for (i = 0; i < 6; i++) {
      Zero[i] = 2. * (pow(1. + Zero[i], 0.5) - 1.);
      ZeroS[i] = Zero[i];
    }
  return;

} /* yldconv */

/*****  coupint  ************************************************************/
/*
 *       Linear interpolation of the swap benchmark yields by 6-month steps.
 */
void coupint(
    double *SwapYield) /* Array of interpolated and benchmark swap yields */
{
  double yield; /* Current interpolated swap yield */
  int IndexSwapYield[] = {4,  6,  8, 10, 14,
                          20, 40, 60}, /* Swap benchmark points in 1/2 years */
      term, /* Current interpolated swap yield maturity in 1/2 years */
      i;

  term = 5;
  for (i = 0; i < 7; i++) {
    do {
      linterp((double)term, &yield, (double)IndexSwapYield[i],
              (double)IndexSwapYield[i + 1], SwapYield[IndexSwapYield[i]],
              SwapYield[IndexSwapYield[i + 1]]);

      SwapYield[term] = yield;

      term++;

    } while (term < IndexSwapYield[i + 1]);

    term++;

  } /* for i */

  return;

} /* coupint */

/*****  spotrate  ***********************************************************/
/*
 *       Bootstrapping of the swap yield curve. Fill the array Zero and
 *       Period with the zeros and their maturity at 6 month interval.
 */
void spotrate(
    double *Zero,      /* output: Zero rates at 6m interval */
    double *ZeroS,     /* output: ZeroS rates at 6m interval */
    long *Period,      /* output: length of periods at 6m interval */
    double *SwapYield, /* array of interpolated and benchmark swap yields */
    double *YieldS,    /* array of interpolated and benchmark yieldS */
    char AoS,          /* 'A'nnual or 'S'emi-Annual frequency */
    long ValueDate,    /* Value date for the benchmark instruments */
    char *YBase)       /* Year basis for the benchmark swaps ("ACT" or "365") */
{
  double temp, temp2, temp3, temp1;
  long newdate; /* Maturity of the current zero */
  int i, j;

  if (AoS == 'A') /* Swaps have annual coupon */
  {
    for (j = 2; j <= MAXMAT; j++) {
      if (YBase[0] == '3') /* In the "365" convention the actual coupon you  */
      {                    /* receive is the quoted one x number of days in  */
        temp1 = (SwapYield[2 * j] * Period[2] /
                 365.) /* the coupon period / 365. For annual it makes   */
                / pow(1. + Zero[2],
                      (double)(Period[2] /
                               365.)); /* a difference only for leap years. */
        temp3 = (YieldS[2 * j] * Period[2] /
                 365.) /* the coupon period / 365. For annual it makes   */
                / pow(1. + ZeroS[2],
                      (double)(Period[2] /
                               365.)); /* a difference only for leap years. */
        for (i = 2; i < j; i++) { /* Note that we use Period[2] instead of */
          temp1 +=
              (SwapYield[2 * j] * (Period[2 * i] - Period[2 * (i - 1)]) /
               365.) /* Period[2*i] - Period[2*(i-1)] when i=1 because */
              /
              pow(1. + Zero[2 * i],
                  (double)(Period[2 * i] / 365.)); /* we need Period[0]=0 which
                                                      is not the case here */
          temp3 +=
              (YieldS[2 * j] * (Period[2 * i] - Period[2 * (i - 1)]) /
               365.) /* Period[2*i] - Period[2*(i-1)] when i=1 because */
              /
              pow(1. + ZeroS[2 * i],
                  (double)(Period[2 * i] / 365.)); /* we need Period[0]=0 which
                                                      is not the case here */
        }                                          /* for */
      } else {
        temp1 = 0.;
        temp3 = 0.;
        for (i = 1; i < j;
             i++) { /* In the "ACT" convention the coupon is always   */
          temp1 += SwapYield[2 * j] /* equal to the quoted one. */
                   / pow(1. + Zero[2 * i], (double)(Period[2 * i] / 365.));
          temp3 += YieldS[2 * j] /* equal to the quoted one. */
                   / pow(1. + ZeroS[2 * i], (double)(Period[2 * i] / 365.));
        }
      } /* if then else */

      newdate = Nxtmth(ValueDate, /* Maturity date of the current zero */
                       (long)j * 12, (long)1);
      Period[2 * j] =
          Daysact(ValueDate, /* Maturity of the current zero in days */
                  newdate);

      if (YBase[0] == '3') {
        Zero[2 * j] =
            pow((100. +
                 SwapYield[2 * j] *
                     (Period[2 * j] -
                      Period[2 * (j - 1)]) /* (1+zero rate)^(days/365) = */
                     / 365.) /
                    (100. - temp1),
                365. / (double)Period[2 * j]) -
            1.; /* (100 + coupon) / (100 - annuity) (cf 3+i). */
        ZeroS[2 * j] =
            pow((100. +
                 YieldS[2 * j] *
                     (Period[2 * j] -
                      Period[2 * (j - 1)]) /* (1+zero rate)^(days/365) = */
                     / 365.) /
                    (100. - temp3),
                365. / (double)Period[2 * j]) -
            1.; /* (100 + coupon) / (100 - annuity) (cf 3+i). */
      } else {
        Zero[2 * j] = pow((100. + SwapYield[2 * j]) / (100. - temp1),
                          365. / (double)Period[2 * j]) -
                      1.;
        ZeroS[2 * j] = pow((100. + YieldS[2 * j]) / (100. - temp3),
                           365. / (double)Period[2 * j]) -
                       1.;
      }

      newdate = Nxtmth(ValueDate, /* Maturity date of the zero maturing 6m
                                     before the current one */
                       (long)(2 * j - 1) * 6, (long)1);
      Period[2 * j - 1] =
          Daysact(ValueDate, /* Maturity date of the zero maturing 6m before the
                                current one  , in days */
                  newdate);
      linterp((double)Period[2 * j - 1], /* Interpolation of the zero maturing
                                            6m before the current one */
              &temp, (double)Period[2 * j - 2], (double)Period[2 * j],
              Zero[2 * j - 2], Zero[2 * j]);
      Zero[2 * j - 1] = temp;

      linterp((double)Period[2 * j - 1], /* Interpolation of the zero maturing
                                            6m before the current one */
              &temp2, (double)Period[2 * j - 2], (double)Period[2 * j],
              ZeroS[2 * j - 2], ZeroS[2 * j]);
      ZeroS[2 * j - 1] = temp2;
    } /* for j */
  } else {
    temp1 = 1. / pow(1. + Zero[1] / 2., /*  Interpolation of the 18m rate */
                     (double)(Period[1] / 182.5));
    temp = 1. / pow(1. + Zero[2] / 2., (double)(Period[2] / 182.5));
    SwapYield[3] = (200. * (1. - temp) / (temp1 + temp) + SwapYield[4]) / 2.;

    temp3 = 1. / pow(1. + ZeroS[1] / 2., /*  Interpolation of the 18m rate */
                     (double)(Period[1] / 182.5));
    temp2 = 1. / pow(1. + ZeroS[2] / 2., (double)(Period[2] / 182.5));
    YieldS[3] = (200. * (1. - temp2) / (temp3 + temp2) + YieldS[4]) / 2.;

    for (j = 3; j <= 2 * MAXMAT; j++) {
      if (YBase[0] == '3') /* In the semi-annual case the "365" convention */
      {                    /* makes a difference as the 6mo coupon will    */
        temp1 =
            (SwapYield[j] * Period[1] / 365.) /* vary from period to period */
            / pow(1. + Zero[1] / 2., (double)(Period[1] / 182.5));
        temp3 =
            (YieldS[j] * Period[1] / 365.) /* vary from period to period */
            / pow(1. + ZeroS[1] / 2., (double)(Period[1] / 182.5));
        for (i = 2; i < j; i++) {
          temp1 += (SwapYield[j] * (Period[i] - Period[i - 1]) / 365.) /
                   pow(1. + Zero[i] / 2., (double)(Period[i] / 182.5));
          temp3 += (YieldS[j] * (Period[i] - Period[i - 1]) / 365.) /
                   pow(1. + ZeroS[i] / 2., (double)(Period[i] / 182.5));
        }
      } else {
        temp1 = 0.; /* temp1 is the price of an annuity paying a coupon of
                       SwapYield[j] every 6m    */
        temp3 = 0.;
        for (i = 1; i < j; i++) {
          temp1 += (SwapYield[j] / 2.) /
                   pow(1. + Zero[i] / 2., (double)(Period[i] / 182.5));
          temp3 += (YieldS[j] / 2.) /
                   pow(1. + ZeroS[i] / 2., (double)(Period[i] / 182.5));
        }
      } /* if then else */

      newdate = Nxtmth(ValueDate, /* Maturity date of the current zero */
                       (long)j * 6, (long)1);
      Period[j] = Daysact(ValueDate, /* Maturity of the current zero in days */
                          newdate);

      if (YBase[0] == '3') {
        Zero[j] =
            2. *
            (pow((100. + SwapYield[j] * (Period[j] - Period[j - 1]) / 365.) /
                     (100. - temp1),
                 (182.5 / (double)Period[j])) -
             1.);
        ZeroS[j] =
            2. * (pow((100. + YieldS[j] * (Period[j] - Period[j - 1]) / 365.) /
                          (100. - temp3),
                      (182.5 / (double)Period[j])) -
                  1.);
      } else {
        Zero[j] = 2. * (pow((100. + SwapYield[j] / 2.) / (100. - temp1),
                            (182.5 / (double)Period[j])) -
                        1.);
        ZeroS[j] = 2. * (pow((100. + YieldS[j] / 2.) / (100. - temp3),
                             (182.5 / (double)Period[j])) -
                         1.);
      }
    } /* for j */
  }   /* if then else */

  return;

} /* spotrate */

/*****  TimeSche  ***********************************************************/
/*
 *       Set up the time step array so that a node falls on each exercise
 *       dates. Set up the accrued interest and exercise arrays. Allocate
 *       memory for the tree arrays.
 */
void Time(long ValueDate,       /* Value date */
          int Ppy,              /* Number of period per year in the tree */
          long Maturity,        /* Maturity of the deal */
          int NbExer,           /* Number of exercise dates */
          long *Exer,           /* Array of exercise dates */
          int NbFwd,            /* Number of forward prices */
          long *FwdMat,         /* Array of forward maturities */
          double *Strike,       /* Array of the corresponding strikes */
          char Freq,            /* Frequency of the underlying cash-flows */
          char EoA,             /* European or American option */
          TREE_DATA *tree_data) /* Output: Structure of tree data */
{
  int NbCoupon,  /* Number of coupon dates from the value date to maturity */
      NbDate,    /* Number of critical dates on which we want a node to fall */
      DateIndex, /* Index of the current critical date */
      ExerIndex, /* Index of the current exercise date position */
      FwdIndex,  /* Index of the current forward maturity position */
      NbNodes,   /* Total number of nodes in the tree */
      Node1,     /* Node position used in interpolation of strike */
      Node2,     /* Node position used in interpolation of strike */
      i, j, k;
  int *NodePerPeriod, /* Number of Nodes between two critical dates */
      *ExerNode,      /* TRUE if exercise node  , FALSE otherwise */
      *FwdNode;       /* TRUE if forward node  , FALSE otherwise */
  double F,           /* Frequency of underlying cash-flows as a double */
      *CDate,         /* Array of time from value date to each critical date */
      *Coupon,        /* Array of time from value date to each coupon date */
      *ExerDate,      /* Array of time from value date to each exercise date */
      Accrued, x;
  long Date, /* Working date */
      Time;

  F = Conv_Freq(Freq);

  NbCoupon = 0; /* Number of coupon dates from the value date to maturity. It
                   also includes the */
  /* coupon falling before today which is used for accrued interest calculation
   */
  do {
    Date = Nxtmth(Maturity, /* Calculate the current coupon date by going
                               backward from maturity */
                  (long)(-12 * NbCoupon++ / F), (long)1);
    Time = Daysact(ValueDate, Date);
  } while (Time > 0);

  NbDate =
      NbExer + NbFwd + NbCoupon + 1; /* The dates are: the value date  , the
                                        exercise dates  , the coupon dates */
  /* and the forward maturities.                                         */
  CDate = mbs_dvector(0, NbDate);
  Coupon = mbs_dvector(0, NbCoupon);
  ExerDate = mbs_dvector(0, NbExer);
  ExerNode = mbs_ivector(0, NbDate);
  FwdNode = mbs_ivector(0, NbDate);
  NodePerPeriod = mbs_ivector(0, NbDate);

  CDate[0] = 0.; /* Value date is the first date and is not an exercise date */

  DateIndex = 1;

  for (i = 0; i < NbExer; i++) /* Exercise dates */
  {
    CDate[DateIndex] = Daysact(ValueDate, /* Time from value date to current
                                             exercise date (ACT/365) */
                               Exer[i]) /
                       365.;
    ExerDate[i] = CDate[DateIndex];

    ExerNode[DateIndex++] = TRUE;
  } /* for i */

  for (i = 0; i < NbFwd; i++) /* Forward maturities */
  {
    CDate[DateIndex] = Daysact(ValueDate, /* Time from value date to current
                                             forward maturity (ACT/365) */
                               FwdMat[i]) /
                       365.;
    FwdNode[DateIndex++] = TRUE;
  } /* for i */

  for (i = 0; i < NbCoupon; i++) /* Coupon dates. The last CDate is only used
                                    for accrued interest calculation. */
  { /* It correspond to the coupon falling before the value date */
    Date = Nxtmth(Maturity, /* Current coupon date */
                  (long)(-12 * i / F), (long)1);
    CDate[DateIndex] = Daysact(ValueDate, /* Time from value date to current
                                             coupon date (ACT/365) */
                               Date) /
                       365.;
    Coupon[i] = CDate[DateIndex++];
  } /* for i */

  NbDate--; /* We get rid of the coupon date falling before the value date */

  pksort3(NbDate, /* Sort the CDate in ascending order */
          CDate, ExerNode, FwdNode);

  for (i = 0, NbNodes = 0, ExerIndex = 0, FwdIndex = 0; i < NbDate - 1;
       i++) /* Total number of periods and number of periods between critical
               dates */
  {
    x = Ppy *
        (CDate[i + 1] - CDate[i]); /* The number of periods between two critical
                                      dates is the smallest integer */
    NodePerPeriod[i] = (int)ceil(
        x); /* that gives a period length smaller than 1/Ppy (in years) */

    if (ExerNode[i]) /* If the current date is an exercise date we store its
                        position */
      tree_data->ExerPosition[ExerIndex++] = NbNodes;

    if (FwdNode[i]) /* If the current date is a forward maturity we store its
                       position */
      tree_data->FwdPosition[FwdIndex++] = NbNodes;

    NbNodes += NodePerPeriod[i];
  } /* for i */

  tree_data->NbNodes = NbNodes;
  tree_data->ExEnd = tree_data->ExerPosition[NbExer - 1];

  Tree_Alloc(tree_data, /* Allocate memory for the arrays constituting the tree.
                           We could not do it */
             NbNodes);  /* before because we did not know the total number of
                           nodes NbNodes.        */

  for (i = 0, k = 0; i < NbDate - 1;
       i++) /* Calculate the length of each time step (period) which is store in
               the global */
  {         /* structure tree_data->Length         */
    if (!NodePerPeriod[i]) /* If two consecutive dates were equal  ,
                              NodePerPeriod[i] = 0. We ignore it      */
      continue;

    x = (CDate[i + 1] - CDate[i]) /
        NodePerPeriod[i]; /* All periods between CDate[i] and CDate[i+1] have
                             the same length x */

    for (j = 0; j < NodePerPeriod[i]; j++)
      tree_data->Length[k++] = x;
  } /* for i */

  tree_data->Length[k] = x; /* We add one more time step which will be used for
                               probability calculation */

  for (i = 0, k = NbNodes; i < NbCoupon - 1;
       i++) /* Calculation of accrued interest at each node */
  {
    x = Coupon[i] - Coupon[i + 1]; /* Time between current coupon and preceding
                                      one  , going backward in time */
    Accrued = 1.;                  /* Coupon[0] is the final maturity                  */
    while ((Accrued > .0001) && (k > 0)) {
      tree_data->Accrued[k] = Accrued;
      Accrued -= tree_data->Length[k - 1] / x;
      k--;
    }
  }

  tree_data->Accrued[0] =
      max(0., Accrued); /* By convention today's accrued is never 1. */

  for (i = 0; i < NbExer;
       i++) /*  Initialization of the strike prices and flags. */
  {
    tree_data->FwdStrike[tree_data->ExerPosition[i]] = Strike[i];
    tree_data->ExerFlag[tree_data->ExerPosition[i]] = 1;
  } /* for i */

  if ((EoA == 'A') && (NbExer > 1)) /* For American option  , we interpolate the
                                       strike linearly between two dates */
    for (i = 0; i < NbExer - 1; i++) {
      Node1 = tree_data->ExerPosition[i];
      Node2 = tree_data->ExerPosition[i + 1];

      if (Node2 > Node1 + 1)
        for (j = Node1 + 1; j < Node2; j++) {
          x = (Strike[i + 1] - Strike[i]) * tree_data->Length[j - 1] /
              (ExerDate[i + 1] - ExerDate[i]);
          tree_data->FwdStrike[j] = tree_data->FwdStrike[j - 1] + x;
          tree_data->ExerFlag[j] = 1;
        } /* for j */
    }     /* for i */

  mbs_free_ivector(NodePerPeriod, 0, NbDate);
  mbs_free_ivector(FwdNode, 0, NbDate);
  mbs_free_ivector(ExerNode, 0, NbDate);
  mbs_free_dvector(ExerDate, 0, NbExer);
  mbs_free_dvector(Coupon, 0, NbCoupon);
  mbs_free_dvector(CDate, 0, NbDate);

  return;

} /* Time */

/*****  Tree_Alloc  *********************************************************/
/*
 *       Allocation of memory for the arrays constituting the tree.
 */
void Tree_Alloc(TREE_DATA *tree_data, int NbNodes) {
  tree_data->ZeroCoupon =
      mbs_dvector(0, NbNodes + 1); /* As we have to calculate probabilities up
                                      to NbNodes  , we need an extra */
  tree_data->FwdRate = mbs_dvector(
      0,
      NbNodes +
          1); /* forward and  , consequently  , an extra zero coupon price. */
  tree_data->FwdRateS = mbs_dvector(0, NbNodes + 1);
  tree_data->Delta = mbs_dvector(0, NbNodes + 1);
  tree_data->ZeroCouponS = mbs_dvector(0, NbNodes + 1);
  tree_data->FwdVol = mbs_dvector(0, NbNodes);
  tree_data->Drift = mbs_dvector(0, NbNodes);
  tree_data->FwdPrepay = mbs_dvector(0, NbNodes);
  tree_data->PrepayVol = mbs_dvector(0, NbNodes);
  tree_data->Length = mbs_dvector(0, NbNodes);
  tree_data->Accrued = mbs_dvector(0, NbNodes);
  tree_data->FwdStrike = mbs_dvector(0, NbNodes);
  tree_data->ExerFlag = mbs_ivector(0, NbNodes);

  return;

} /* Tree_Alloc */

/*****  Build_Tree  *********************************************************/
/*
 *       Main routine for the construction of the tree.
 */
char *Build_Tree(TERM_DATA *term_data, /* Structure of term structure data */
                 TREE_DATA *tree_data, /* Structure of tree data */
                 PREPAY_DATA *prepay_data) /* Structure of prepayment data */
{
  char *err;

  Forward(
      tree_data->ZeroCoupon, /* Calculation of zero coupon price at each node in
                                the tree */
      tree_data
          ->FwdRate, /* Calculation of forward rate at each node in the tree */
      term_data->Zero, term_data->Period, tree_data->Length, tree_data->NbNodes,
      term_data->AoS);

  //**NOT USED**
  Forward(
      tree_data->ZeroCouponS, /* Calculation of zero coupon price at each node
                                 in the tree */
      tree_data
          ->FwdRateS, /* Calculation of forward rate at each node in the tree */
      term_data
          ->ZeroS, /* These zero coupons are for the porpose of discounting */
      term_data->Period, tree_data->Length, tree_data->NbNodes, term_data->AoS);

  //**NOT USED**
  Forward_Prepay(tree_data->FwdPrepay, /* Calculation of forward prepayment at
                                          each node in the tree */
                 prepay_data->Spot, tree_data->NbNodes);

  if (err = Forward_Vol(
          tree_data->FwdVol, /* Calculation of forward rate volatility at each
                                node in the tree */
          &(tree_data->Jump), tree_data->FwdRate, tree_data->Length,
          term_data->Zero, term_data->Period, term_data->SwaptionVol,
          term_data->SwaptionExp, term_data->Beta, term_data->IndexMat,
          term_data->IndexFreq, term_data->AoS, tree_data->Ppy,
          tree_data->NbNodes))
    return (nrerror(err));

  //**NOT USED**
  Prepay_Vol(tree_data->PrepayVol, /* Calculation of instantaneous volatility of
                                      forward prepayment */
             prepay_data->SpotVol, tree_data->NbNodes);

  //**NOT USED**
  Conv_Adj(
      tree_data
          ->FwdRate, /* Calculate convexity adjustment for discounting rates */
      tree_data->FwdRateS, tree_data->FwdVol, tree_data->Length,
      tree_data->Delta, tree_data->NbNodes, term_data->Beta);

  if (err = Drift(tree_data->Drift, /* Calculation of the 1 period rate drift in
                                       the tree */
                  &(tree_data->NodeMax), tree_data->Jump, tree_data->ZeroCoupon,
                  tree_data->FwdRate, tree_data->FwdVol, tree_data->Length,
                  term_data->Beta, tree_data->NbNodes))
    return (nrerror(err));

  //**NOT USED**
  NodeMax2(&(tree_data->Jump2), /* Calculate the jump in the prepayment tree and
                                   the maximum node ever attained */
           &(tree_data->NodeMax2), prepay_data->SpotVol, tree_data->PrepayVol,
           prepay_data->Beta, tree_data->Length, tree_data->NbNodes);
  return (0);
} /* Build_Tree */

/*****  Forward  ************************************************************/
/*
 *       Interpolate zero and calculate forward at each node point i.e.
 *       every Length[i] interval.
 */
void Forward(
    double *ZeroCoupon, /* Output: zero coupon price at each node in the tree */
    double *FwdRate,    /* Output: forward rate at each node in the tree */
    double *Zero,       /* Zero rates at money market points and 6m interval */
    long *Period,       /* Length of money market periods and 6m intervals */
    double *Length,     /* Length of each time step */
    int NbNodes,        /* Total number of nodes (or time steps) */
    char AoS)           /* Yield curve type (annual or semi-annual) */
{
  double Days,   /* Maturity of the current zero in days from value date */
      zero,      /* Interpolated zero rate value */
      Discount1, /* Discount factor at the beginning of the current time step */
      Discount2, /* Discount factor at the end of the current time step */
      maturity,  /* Maturity of the current zero (ACT/365) */
      Freq;      /* Coupon frequency of yield curve */
  int i, k;

  for (i = 0, k = 1, Days = 0.; i <= NbNodes;
       i++) /* Interpolation of zero rates at each node (time step). Note that
               because we need */
  { /* to calculate probabilities up to NbNodes we need one extra forward and */
    Days +=
        Length[i] * 365.; /* consequently one extra zero at node NbNodes+1 */

    while ((Days > Period[k]) && (k < 2 * MAXMAT + 3))
      k++;

    linterp(Days, &zero, (double)Period[k - 1], (double)Period[k], Zero[k - 1],
            Zero[k]);

    FwdRate[i + 1] = zero; /* For the time being FwdRate[i] is not the forward
                              rate but the zero rate */

  } /* for i */

  FwdRate[0] = FwdRate[1]; /* Artificial value for the initial node */

  Freq = Conv_AoS(AoS);

  for (i = 0, maturity = 0.; i <= NbNodes;
       i++) /* Calculation of forward rates from interpolated zero rates */
  {
    Discount1 = pow(1. + FwdRate[i] / Freq, Freq * maturity);
    ZeroCoupon[i] =
        1. / Discount1; /* We store the prices of zero coupon maturing at each
                           node. We will use it later */

    maturity += Length[i];

    Discount2 = pow(1. + FwdRate[i + 1] / Freq, Freq * maturity);

    FwdRate[i] = Discount2 / Discount1 - 1.;
  } /* for i */

  ZeroCoupon[NbNodes + 1] = 1. / Discount2;

  return;

} /* Forward */

/*****  Forward_Prepay  *****************************************************/
/*
 *       Interpolate forward prepayment at each node point i.e. every
 *       Length[i] interval.
 */
void Forward_Prepay(
    double *FwdPrep, /* Output: prepayment at each node in the tree */
    double Spot,     /* Spot prepayment */
    int NbNodes)     /* Total number of nodes (or time steps) */
{
  int i;

  for (i = 0; i <= NbNodes; i++)
    FwdPrep[i] = Spot; /* The prepayment has no drift for the moment being */

  return;

} /* Forward_Prepay */

/*****  Forward_Vol  ********************************************************/
/*
 *       Calculate forward rate volatility at each node in the tree as well as
 *       the jump size.
 *       We have NBSWAPTION swaption volatility inputs. All of the swaptions
 *       are CMS on a swap with the same maturity as the index. They expire at
 *       different points given by SwaptionExp. From these volatilities we can
 *       get the base volatility at each SwaptionExp
 *       in a one mean reverting factor model. Then we assume that the spot
 *       volatility is constant between SwaptionExp  , and calculate it.
 */
char *Forward_Vol(
    double
        *FwdVol,  /* Output: forward rate volatility at each node in the tree */
    double *Jump, /* Output: Size of the jump (in log space) */
    double *FwdRate,     /* Forward rate at each node */
    double *Length,      /* Length of each time step */
    double *Zero,        /* Zero rates at money market points and 6m interval */
    long *Period,        /* Length of money market periods and 6m intervals */
    double *SwaptionVol, /* Volatility of the benchmark swaptions */
    int *SwaptionExp,    /* Expiration (in years) of the benchmark swaptions */
    double Beta,         /* Mean reversion coefficient */
    int IndexMat,        /* Index maturity in years */
    char IndexFreq,      /* Index frequency */
    char AoS,            /* Yield curve type (annual or semi-annual) */
    int Ppy,             /* Number of periods per year in the tree */
    int NbNodes)         /* Total number of nodes (or time steps) */
{
  double BaseVol[NBSWAPTION + 1], /* Base volatility (squared) at expiration of
                                     each of the benchmark swaptions */
      SpotVol[NBSWAPTION + 1], /* Spot volatility at expiration of each of the
                                  benchmark swaptions */
      FwdYield, /* Yield of the forward swap underlying the current swaption */
      *FwdZeroYield, /* Yield of the zero going from expiration of the swaption
                        to one of the coupons */
      FwdZero, /* Value of the forward zero going from expiration to one of the
                  current coupon */
      Annuity, /* Price of an annuity paying $1 at each coupon of the forward
                  swap */
      MaxSpotVol, /* Maximum of the spot volatility. Used to calculate the jump
                     size */
      Freq,       /* Coupon frequency of yield curve as a double */
      Time,       /* Time position of the current node in the tree */
      T1, T2,     /* Time periods used for base volatility */
      x, y;       /* Doubles used for intermediate calculations */
  int IndexF,     /* Frequency of index payment as a integer */
      NbCoupon,   /* Number of coupon in the CMS swap */
      ZeroIndex,  /* Index of the current spot zero used to calculate forward
                     zero */
      i, j;
  char string[MAXBUFF];

  Freq = Conv_AoS(AoS);
  IndexF = Conv_Freq(IndexFreq);
  NbCoupon = IndexMat * IndexF; /* Each of the CMT swaption is going to have
                                   IndexMat*IndexF coupons */

  FwdZeroYield = mbs_dvector(0, NbCoupon);

  BaseVol[0] = 0.;

  for (i = 0; i < NBSWAPTION; i++) {
    for (j = 1, Annuity = 0.; j <= NbCoupon; j++) {
      ZeroIndex = 2 * SwaptionExp[i] +
                  3; /* Index of the zero maturing at the start of the */
      FwdZero =
          pow(1. + Zero[ZeroIndex] / Freq,
              Freq * Period[ZeroIndex] /
                  365.); /* forward swap. Note the 3 offset due to money   */
                         /* market instruments.                            */
      ZeroIndex = 2 * SwaptionExp[i] + 3 +
                  j; /* Index of the zero maturing at coupon j of the */
      FwdZero /=
          pow(1. + Zero[ZeroIndex] / Freq,
              Freq * Period[ZeroIndex] /
                  365.); /* forward swap underlying the current swaption. */
                         /* FwdZero is the long zero divided by the short */
      Annuity += FwdZero;

      FwdZeroYield[j] = IndexF * (pow(FwdZero, -1. / j) - 1.);

    } /* for j */
      /* Yield of the forward swap. Note that FwdZero has the value of the last
       * zero */
    FwdYield = IndexF * (1. - FwdZero) /
               Annuity; /* calculated above i.e. the one maturing at maturity of
                           the forward swap.     */

    for (j = 1, x = y = 0.; j <= NbCoupon;
         j++) /* Calculation of the intermediate doubles x & y */
    {
      T1 = Period[2 * SwaptionExp[i] + 3] /
           365.; /* Maturity of the current base volatility point (i) */
      T2 = Period[2 * SwaptionExp[i] + 3 + j] /
           365.; /* Maturity of the current coupon (j) */

      x += FwdYield * pow(1. + FwdYield / IndexF, -j - 1.) * j * FwdYield;
      y += FwdYield * pow(1. + FwdZeroYield[j] / IndexF, -j - 1.) * j *
           FwdZeroYield[j] * ExpDecay(Beta, T2 - T1);

    } /* for j */

    x += pow(1. + FwdYield / IndexF, -NbCoupon - 1.) * NbCoupon * FwdYield;
    y += pow(1. + FwdZeroYield[NbCoupon] / IndexF, -NbCoupon - 1.) * NbCoupon *
         FwdZeroYield[NbCoupon] * ExpDecay(Beta, T2 - T1);

    BaseVol[i + 1] = pow(SwaptionVol[i] / 100. * x / y,
                         2.); /* use the relation between the base vol  ,
                                 swaption vol  , x & y */

  } /* for i */

  for (i = 1, T1 = 0, MaxSpotVol = 0.; i <= NBSWAPTION;
       i++) /* The i index is now running from 1 to NBSWAPTION */
  {
    T2 = Period[2 * SwaptionExp[i - 1] + 3] /
         365.; /* Maturity of the current base volatility point */

    SpotVol[i] =
        (T2 * BaseVol[i] - T1 * BaseVol[i - 1] * exp(-2. * Beta * (T2 - T1))) /
        (T2 - T1) / ExpDecay(2. * Beta, T2 - T1);

    if (SpotVol[i] < 0.00001 || SpotVol[i] > 1.0) {
      sprintf(string, "SpotVol[%d]=%f!", i, SpotVol[i]);
      mbs_free_dvector(FwdZeroYield, 0, NbCoupon);
      return (nrerror(string));
    }

    SpotVol[i] = sqrt(SpotVol[i]);

    MaxSpotVol = max(MaxSpotVol, SpotVol[i]);

    T1 = T2;

  } /* for i */

  Time = 0.;
  i = 1;
  j = 0;

  while ((i <= NBSWAPTION) && (j <= NbNodes)) {
    FwdVol[j] = SpotVol[i]; /* The spot volatility is constant between the
                               benchmark swaption points */

    Time += Length[j++];

    if (Time > Period[2 * SwaptionExp[i - 1] + 3] / 365. -
                   MBS_ERROR) /* If Time is greater than the swaption expiration
                                 we switch to the next swaption */
      i++;
  } /* while */

  for (; j <= NbNodes;
       j++) /* If there are nodes after the last benchmark point we fill the */
    FwdVol[j] =
        SpotVol[NBSWAPTION]; /* spot volatility with the last value found. */

  *Jump = MaxSpotVol *
          sqrt(JUMPCOEFF / (double)Ppy); /* The jump size depends on maximum
                                            volatility and time step  , 1/ppy */
  /* so that we avoid negative probabilities.                         */

  if (Beta > MBS_ERROR) /* Because FwdRate is not an instantaneous but a    */
    for (i = 0; i <= NbNodes;
         i++) /* Length[i] rate the volatility has to be adjusted */
      FwdVol[i] *= (1. + FwdRate[i]) * log(1. + FwdRate[i]) / FwdRate[i] *
                   (1. - exp(-Beta * Length[i])) / (Beta * Length[i]);
  else /* If Beta=0 we only adjust for compounding frequency. */
    for (i = 0; i <= NbNodes; i++)
      FwdVol[i] *= (1. + FwdRate[i]) * log(1. + FwdRate[i]) / FwdRate[i];

  mbs_free_dvector(FwdZeroYield, 0, NbCoupon);

  return (0);

} /* Forward_Vol */

/*****  Prepay_Vol  *********************************************************/
/*
 *       Interpolate forward stock volatility at each node point i.e. every
 *       Length[i] interval.
 */
void Prepay_Vol(double *PrepayVol, /* Output: instantaneous volatility of
                                      forward prepayment */
                double SpotVol,    /* Spot volatility of prepayment */
                int NbNodes)       /* Total number of nodes (or time steps) */
{
  int i;

  for (i = 0; i <= NbNodes; i++)
    PrepayVol[i] =
        SpotVol / 100.; /* The spot volatility is constant for the time being */

  return;

} /* Prepay_Vol */

/******Conv_Adj*************************************************************/
/*
 *       Calculate convexity adjustment when discounting rates and
 *     rates on tree nodes differ.
 *
 */
void Conv_Adj(
    double *FwdRate, /* Calculate convexity adjustment for discounting rates */
    double *FwdRateS, double *FwdVol, double *Length, double *Delta,
    int NbNodes, double Beta) {
  int i, j;
  double timelapse = 0.0;

  Delta[0] = 0.0;

  for (i = 1; i <= NbNodes; i++) {

    Delta[i] =
        Delta[i - 1] *
        (1.0 -
         Beta * Length[i]); /* multiplicative correction due to mean rev */

    for (j = i - 1; j >= 0; j--) {

      timelapse += Length[j]; /* time between cuurent (i) and j */

      Delta[i] += Length[i] * FwdVol[j] * FwdVol[j] * /* sum over history */
                  exp(-2. * Beta * timelapse) * (FwdRateS[i] - FwdRate[i]);

    } /* for j */

    timelapse = 0.0;

  } /* for i */
} /*Conv_Adj*/

char *Drift_new(TREE_DATA *old_tree, MBS_BKTree *new_tree) {
  char *err;
  int i;
  if (err = MBS_BKTree_SolveDrift(new_tree))
    return (err);
  for (i = 0; i < new_tree->NbTimeSlices; ++i)
    old_tree->Drift[i] = new_tree->Drift[i];
  old_tree->NodeMax = new_tree->NodeMax;
  return (0);
}
/*****  Drift  **************************************************************/
/*
 *       Calculate the interest rate drift
 */
char *Drift(
    double *Drift, /* Output: drift of the 1 period interest rate in the tree */
    int *NodeMax,  /* Output: maximum node reached in the tree (both top and
                      bottom) */
    double Jump,   /* Size of the jump (in log space) */
    double *ZeroCoupon, /* Zero coupon price at each node in the tree */
    double *FwdRate,    /* Forward rate at each node in the tree */
    double *FwdVol,     /* Forward rate volatility at each node in the tree */
    double *Length,     /* Length of each time step */
    double Beta,        /* Mean reversion coefficient (constant) */
    int NbNodes)        /* Total number of nodes (or time steps) */
{
  double Vol,   /* Volatility of the short term rate over the period */
      FwdRate1, /* 1 period rate in one period */
      pu,       /* Probability in the up state */
      pd,       /* Probability in the down state */
      p0,       /* Probability in the flat state */
      x1, x2, x3, y1, y2, y3, z1, z2, z3, a, b,
      c; /* Doubles used in the intermediate calculations of the drift */
  double *Discount, /* 1 period discount factor at each node (i.e.
                       1/(1+FwdRate[i]) */
      *StatePr,   /* State Prices (i.e. Discounted Expected Value of $1 received
                     at this node) */
      *Discount1, /* 1 period discount factor  , in 1 period at each node */
      *StatePr1;  /* State Prices in one period */
  int Bottom,     /* Bottom node index in the tree */
      Top,        /* Top node index in the tree */
      CutFlag,    /* =TRUE if tree is already cut  , =FALSE otherwise */
      i,          /* Node index (i.e. vertical) */
      k, /* If k=0 the tree is not cut; if k=1 the tree is cut and branches up
            or down */
      NbPer; /* Time step index (i.e. horizontal) */
  char string[MAXBUFF];

  Discount = mbs_dvector(0, 2 * NbNodes + 2);
  StatePr = mbs_dvector(0, 2 * NbNodes + 2);
  Discount1 = mbs_dvector(0, 2 * NbNodes + 2);
  StatePr1 = mbs_dvector(0, 2 * NbNodes + 2);

  if (*NodeMax < 10)
    *NodeMax = 10000; /* 10000 = Big! 10 = Small! */
  CutFlag = 0;
  NbPer = 0;
  StatePr[0] = 1.;
  Discount[0] = ZeroCoupon[1];
  do {
    Vol = FwdVol[NbPer] * sqrt(Length[NbPer]); /* Volatility over the period =
                                                  sigma * sqrt(time step) */

    Bottom =
        (int)max(0., NbPer - *NodeMax); /* If NbPer>NodeMax the tree will be cut
                                           at the top and the bottom */
    Top = (int)min(2 * NbPer, NbPer + *NodeMax);

    for (i = Bottom; i <= Top + 2;
         i++) /* We don't need all the nodes if the tree is cut */
    {
      FwdRate1 =
          FwdRate[NbPer + 1] *
          exp(Jump *
              (i - (NbPer + 1))); /* Value of the short term rate at node i */
      Discount1[i] =
          1. /
          (1. +
           FwdRate1); /* 1 period discount factor in one period at node i */
    }                 /* for i */

    x1 = x2 = x3 = y1 = y2 = y3 = z1 = z2 = z3 =
        0.; /* Intermediate calculations for the drift */
    for (i = Bottom + 1; i <= Top - 1;
         i++) /* We don't include the top and bottom of the tree: it is dealt
                 with afterwards */
    {
      x1 += StatePr[i] * Discount[i] * Discount1[i + 2];
      x2 += StatePr[i] * Discount[i] * Discount1[i + 2] * (i - NbPer);
      x3 += StatePr[i] * Discount[i] * Discount1[i + 2] * (i - NbPer) *
            (i - NbPer);
      y1 += StatePr[i] * Discount[i] * Discount1[i];
      y2 += StatePr[i] * Discount[i] * Discount1[i] * (i - NbPer);
      y3 += StatePr[i] * Discount[i] * Discount1[i] * (i - NbPer) * (i - NbPer);
      z1 += StatePr[i] * Discount[i] * Discount1[i + 1];
      z2 += StatePr[i] * Discount[i] * Discount1[i + 1] * (i - NbPer);
      z3 += StatePr[i] * Discount[i] * Discount1[i + 1] * (i - NbPer) *
            (i - NbPer);
    } /* for i */

    c = pow(Vol / Jump, 2.) * (.5 * x1 + .5 * y1 - z1) +
        .5 * Beta * Length[NbPer] * (y2 - x2) +
        pow(Beta * Length[NbPer], 2.) * (.5 * x3 + .5 * y3 - z3) + z1 -
        ZeroCoupon[NbPer + 2];
    b = .5 * x1 - .5 * y1 + Beta * Length[NbPer] * (2. * z2 - y2 - x2);
    a = .5 * x1 + .5 * y1 - z1;

    for (i = Bottom; i <= Top;
         i += Top - Bottom) /* 2 cases: Top and Bottom of the tree */
    {
      k = (NbPer + 1 > *NodeMax) *
          ((i == Bottom) -
           (i == Top)); /* If the number of nodes in the next period is greater
                           than the maximum   */
      /* number  , the tree is cut i.e. the top node (i=Top) branches down */
      /* (i.e. to Top+1  , Top  , Top-1 instead of Top+2  , Top+1  , Top) and
       * the bottom */
      /* node branches up (i.e. to bottom+1  , bottom+2  , bottom+3). */
      /* The shift in branching is given by k. */

      c += StatePr[i] * Discount[i] *
           (.5 * Discount1[i + 2 + k] *
                (k * (k - 1.) + pow(Vol / Jump, 2.) +
                 (2. * k - 1.) * Beta * Length[NbPer] * (i - NbPer) +
                 pow(Beta * Length[NbPer] * (i - NbPer), 2.)) -
            Discount1[i + 1 + k] *
                (k * k - 1. + pow(Vol / Jump, 2.) +
                 2. * k * Beta * Length[NbPer] * (i - NbPer) +
                 pow(Beta * Length[NbPer] * (i - NbPer), 2.)) +
            .5 * Discount1[i + k] *
                (k * (k + 1.) + pow(Vol / Jump, 2.) +
                 (2. * k + 1.) * Beta * Length[NbPer] * (i - NbPer) +
                 pow(Beta * Length[NbPer] * (i - NbPer), 2.)));

      b += StatePr[i] * Discount[i] *
           (Discount1[i + 2 + k] *
                ((1. - 2. * k) / 2. - Beta * Length[NbPer] * (i - NbPer)) +
            Discount1[i + 1 + k] *
                (2. * k + 2. * Beta * Length[NbPer] * (i - NbPer)) -
            Discount1[i + k] *
                ((1. + 2. * k) / 2. + Beta * Length[NbPer] * (i - NbPer)));

      a += StatePr[i] * Discount[i] *
           (.5 * Discount1[i + 2 + k] - Discount1[i + 1 + k] +
            .5 * Discount1[i + k]);

      if (Top ==
          0) /* For NbPer=0  , Top=Bottom=0 i.e. there is only one node. */
        break;
    } /* for i */

    if (b * b - 4. * a * c < 0. || fabs(a) < 1E-14) {
      sprintf(string, "Quadratic equation for Drift[%d] has no real root",
              NbPer);
      mbs_free_dvector(StatePr1, 0, 2 * NbNodes + 2);
      mbs_free_dvector(Discount1, 0, 2 * NbNodes + 2);
      mbs_free_dvector(StatePr, 0, 2 * NbNodes + 2);
      mbs_free_dvector(Discount, 0, 2 * NbNodes + 2);
      return (nrerror(string));
    }

    Drift[NbPer] = (-b - sqrt(b * b - 4. * a * c)) / 2. / a;

    for (i = Bottom; i <= Top; i++) /* Calculate state prices in one period */
    {
      k = (NbPer + 1 > *NodeMax) * ((i == Bottom) - (i == Top));

      pu = .5 * (k * (k - 1.) + pow(Vol / Jump, 2.) +
                 (1. - 2. * k) *
                     (Drift[NbPer] - Beta * Length[NbPer] * (i - NbPer)) +
                 pow(Drift[NbPer] - Beta * Length[NbPer] * (i - NbPer), 2.));
      pd = .5 * (k * (k + 1.) + pow(Vol / Jump, 2.) -
                 (1. + 2. * k) *
                     (Drift[NbPer] - Beta * Length[NbPer] * (i - NbPer)) +
                 pow(Drift[NbPer] - Beta * Length[NbPer] * (i - NbPer), 2.));
      p0 = 1. - pu - pd;

      StatePr1[i + 2 + k] +=
          pu * StatePr[i] *
          Discount[i]; /* There are 3 branches originating from node i and they
                          contribute to 3 state */
      StatePr1[i + 1 + k] +=
          p0 * StatePr[i] * Discount[i]; /* prices in one period. */
      StatePr1[i + k] += pd * StatePr[i] * Discount[i];
    } /* for i */

    for (i = Bottom; i <= Top + 4;
         i++) /* Permute the values so that they can be used for the next
                 iteration of the loop. */
    {         /* Again we don't use all the nodes if the tree is cut.         */
      Discount[i] = Discount1[i];
      StatePr[i] = StatePr1[i];
      Discount1[i] = StatePr1[i] = 0.;
    } /* for i */

    if (!CutFlag) /* If tree has not been cut yet we check that the
                     probabilities at the top and */
    {             /* bottom of the tree are not too big or too small.             */
      pu = .5 *
           (pow(Vol / Jump, 2.) + Drift[NbPer] + Beta * Length[NbPer] * NbPer +
            pow(Drift[NbPer] + Beta * Length[NbPer] * NbPer, 2.));
      pd = .5 *
           (pow(Vol / Jump, 2.) - Drift[NbPer] - Beta * Length[NbPer] * NbPer +
            pow(Drift[NbPer] + Beta * Length[NbPer] * NbPer, 2.));
      p0 = 1. - pu - pd;

      if ((pu > .999) || (p0 < .001) ||
          (pd < .001)) /* Bottom of the tree (i=0): pu might be close to 1 and
                          pd close to 0 */
      {
        CutFlag = TRUE;
        *NodeMax = NbPer + 1; /* We cut the tree at the next period */
      }                       /* if */

      pu = .5 *
           (pow(Vol / Jump, 2.) + Drift[NbPer] - Beta * Length[NbPer] * NbPer +
            pow(Drift[NbPer] - Beta * Length[NbPer] * NbPer, 2.));
      pd = .5 *
           (pow(Vol / Jump, 2.) - Drift[NbPer] + Beta * Length[NbPer] * NbPer +
            pow(Drift[NbPer] - Beta * Length[NbPer] * NbPer, 2.));
      p0 = 1. - pu - pd;

      if ((pu < .001) || (p0 < .001) ||
          (pd > .999)) /* Top of the tree (i=2*NbPer): pd might be close to 1
                          and pu close to 0 */
      {
        CutFlag = TRUE;
        *NodeMax = NbPer + 1;
      } /* if */
    }   /* if */

    NbPer++;

  } while (NbPer <=
           NbNodes - 1); /* We need Drift[NbNodes-1] => We need FwdRate[NbNodes]
                            , FwdVol[NbNodes-1] and */
                         /* Discount[NbNodes+1].                        */
  Drift[NbPer] = 0.;

  mbs_free_dvector(StatePr1, 0, 2 * NbNodes + 2);
  mbs_free_dvector(Discount1, 0, 2 * NbNodes + 2);
  mbs_free_dvector(StatePr, 0, 2 * NbNodes + 2);
  mbs_free_dvector(Discount, 0, 2 * NbNodes + 2);

  return (0);

} /* Drift */

/*****  NodeMax2  ***********************************************************/
/*
 *       Calculate the jump in the prepayment tree and the maximum node ever
 *       attained.
 */
void NodeMax2(
    double *Jump2,     /* Output: size of the prepayment jump (in log space) */
    int *NodeMax2,     /* Output: maximum node reached in the tree (both top and
                          bottom) */
    double SpotVol,    /* Spot volatility of prepayment */
    double *PrepayVol, /* Instantaneous volatility of forward prepayment */
    double Beta,       /* Mean reversion of prepayment */
    double *Length,    /* Length of each time step */
    int NbNodes)       /* Total number of nodes (or time steps) */
{
  double pu,   /* Probability in the up state */
      pd,      /* Probability in the down state */
      p0,      /* Probability in the flat state */
      Vol,     /* Volatility of the prepayment over the period */
      Drift;   /* Drift of the prepayment */
  int CutFlag, /* =TRUE if tree is already cut  , =FALSE otherwise */
      NbPer;   /* Time step index (i.e. horizontal) */

  *Jump2 =
      SpotVol / 100. *
      sqrt(JUMPCOEFF / JUMPCOEFF2); /* The jump in the second dimension is not
                                       linked to Ppy but to JUMPCOEFF2 */

  *NodeMax2 = 10000; /* 10000 = Big! */

  if (Beta < MBS_ERROR) /* If Beta = 0 the tree will not be cut */
    return;

  CutFlag = 0;
  NbPer = 0;
  /* If tree has not been cut yet we check that the probabilities at the top and
   */
  while ((!CutFlag) &&
         (NbPer <
          NbNodes)) /* bottom of the tree are not too big or too small. */
  {
    Vol =
        PrepayVol[NbPer] * sqrt(Length[NbPer]); /* Volatility over the period =
                                                   sigma * sqrt(time step) */
    Drift = -.5 * Vol * Vol / *Jump2;

    pu = .5 * (pow(Vol / *Jump2, 2.) + Drift + Beta * Length[NbPer] * NbPer +
               pow(Drift + Beta * Length[NbPer] * NbPer, 2.));
    pd = .5 * (pow(Vol / *Jump2, 2.) - Drift - Beta * Length[NbPer] * NbPer +
               pow(Drift + Beta * Length[NbPer] * NbPer, 2.));
    p0 = 1. - pu - pd;

    if ((pu > .999) || (p0 < .001) ||
        (pd < .001)) /* Bottom of the tree (i=0): pu might be close to 1 and pd
                        close to 0 */
    {
      CutFlag = TRUE;
      *NodeMax2 = NbPer + 1; /* We cut the tree at the next period */
    }                        /* if */

    pu = .5 * (pow(Vol / *Jump2, 2.) + Drift - Beta * Length[NbPer] * NbPer +
               pow(Drift - Beta * Length[NbPer] * NbPer, 2.));
    pd = .5 * (pow(Vol / *Jump2, 2.) - Drift + Beta * Length[NbPer] * NbPer +
               pow(Drift - Beta * Length[NbPer] * NbPer, 2.));
    p0 = 1. - pu - pd;

    if ((pu < .001) || (p0 < .001) ||
        (pd > .999)) /* Top of the tree (i=2*NbPer): pd might be close to 1 and
                        pu close to 0 */
    {
      CutFlag = TRUE;
      *NodeMax2 = NbPer + 1;
    } /* if */

    NbPer++;

  } /* while */

  return;

} /* NodeMax2 */

/*****  Print_Term  *********************************************************/
/*
 *       Print term structure in an ascii file.
 */
void Print_Term(TERM_DATA *term_data, TREE_DATA *tree_data) {
  long date;
  int i;
  double Forward, /* Forward rate for the current period */
      ZeroRate, /* Yield of the zero maturing at the current node on an ACT/365
                   basis */
      days,     /* Number of days from today to the current node */
      discount, Freq; /* Coupon frequency of yield curve */
  FILE *stream;

  stream = 0; // fopen ("TERM.prn"  , "w");

  Freq = Conv_AoS(term_data->AoS);

  fprintf(stream, "Maturity    Zero Rate     Discount Factor\n\n");

  for (i = 4; i < 2 * MAXMAT + 4; i++) {
    date = Nxtmth(term_data->ValueDate, /* Maturity date of the current zero */
                  (long)((i - 3) * 6), (long)1);
    discount =
        1. / pow(1. + term_data->Zero[i] /
                          Freq, /* Discount factor up to the current date */
                 (double)(term_data->Period[i] / 365. * Freq));
    fprintf(stream, "%ld    %9.6lf     %6.4lf\n", date,
            term_data->Zero[i] * 100., /* Zero rates at 6m interval */
            discount);
  } /* for i */

  fprintf(stream, "\nNode    days   forward  zero discount   Drift   fwdvol  "
                  "Acc  flag strike\n");

  for (i = 0, days = 0., discount = 1., ZeroRate = 0.; i <= tree_data->NbNodes;
       i++) {
    Forward = 100. *
              (pow(1. + tree_data->FwdRate[i], /* Conversion: 1.+FwdRate[i]=(1.+Forward/Freq)^(Length[i]*Freq/365.)
                                                */
                   1. / tree_data->Length[i] / Freq) -
               1.) *
              Freq;

    fprintf(
        stream,
        "[%3d] %8.2lf %7.4lf %7.4lf %6.4lf %9.6lf %6.3lf %5.3lf  %d   %6.2lf\n",
        i, days, Forward, ZeroRate, tree_data->ZeroCoupon[i],
        tree_data->Drift[i], tree_data->FwdVol[i] * 100., tree_data->Accrued[i],
        tree_data->ExerFlag[i], tree_data->FwdStrike[i]);

    days +=
        tree_data->Length[i] * 365.; /* Length was stored as # of days / 365. */
    ZeroRate =
        100. *
        (pow(1 / tree_data->ZeroCoupon[i + 1], 365. / Freq / days) - 1.) * Freq;
  } /* for i */

  fclose(stream);

} /* Print_Term */

void Lattice_new(double *Discount, /* Output: one period discount factor at each
                                      node at period NbPer */
                 double *DiscountS, /* Output: one period discount factor
                                       including Option Adjusted spread */
                 double *pu, /* Output: array of probabilities in the up state
                                in the interest rate tree */
                 double *pd, /* Output: array of probabilities in the down state
                                in the interest rate tree */
                 double *p0, /* Output: array of probabilities in the flat state
                                in the interest rate tree */
                 int NbPer,  /* Current time period */
                 MBS_BKTree *tree) {
  int sz, i;
  double *ptr;
  double probs[3];
  int Bottom, Top;
  Bottom = tree->Bottom[NbPer];
  Top = tree->Top[NbPer];

  MBS_BKTree_Calc1PrdDF(&sz, &ptr, 0, NbPer, tree);
  for (i = Bottom; i <= Top; ++i)
    Discount[i] = ptr[i];

  MBS_BKTree_Calc1PrdDF(&sz, &ptr, 1, NbPer, tree);
  for (i = Bottom; i <= Top; ++i)
    DiscountS[i] = ptr[i];

  for (i = Bottom; i <= Top; ++i) {
    MBS_BKTree_GetProb(probs, 0, tree, NbPer, i);
    pu[i] = probs[2];
    p0[i] = probs[1];
    pd[i] = probs[0];
  }
}

/*****  Lattice  ************************************************************/
/*
 *       Calculates the arrays of probabilities and the one period discount
 *       factor in the one factor tree at the current time period.
 */
void Lattice(
    double *Discount,  /* Output: one period discount factor at each node at
                          period NbPer */
    double *DiscountS, /* Output: one period discount factor including Option
                          Adjusted spread */
    double *pu,        /* Output: array of probabilities in the up state in the
                          interest rate tree */
    double *pd, /* Output: array of probabilities in the down state in the
                   interest rate tree */
    double *p0, /* Output: array of probabilities in the flat state in the
                   interest rate tree */
    int NbPer,  /* Current time period */
    int N,      /* Number of additional periods at the beginning of the tree */
    TERM_DATA *term_data, /* Structure of term structure data */
    TREE_DATA *tree_data) /* Structure of tree data */
{
  double FwdYield, /* Value of the short term rate at current node */
      Drift,  /* Drift of the 1 period interest rate at the current period */
      Length, /* Length of the current time step */
      Beta,   /* interest rate mean reversion coefficient */
      Vol,    /* Volatility of the short term rate over the period = SigmaR *
                 sqrt(time step) */
      Multiplier, /* multiplier relating tree rates and discounting rates */
      Jump;       /* Size of the jump (in log space) */
  int Bottom,     /* Bottom node index in the interest rate tree */
      Top,        /* Top node index in the interest rate tree */
      i,          /* Node index for the interest rate tree */
      k; /* If k=0 the tree is not cut; if k=1 the tree is cut and branches up
            or down */

  Drift = tree_data->Drift[NbPer - N];
  Length = tree_data->Length[NbPer - N];
  Beta = term_data->Beta;
  Jump = tree_data->Jump;
  Vol = tree_data->FwdVol[NbPer - N] * sqrt(Length);
  Bottom = (int)max(
      0., NbPer - tree_data->NodeMax); /* If NbPer>NodeMax the tree is cut at
                                          the top and the bottom */
  Top = (int)min(2 * NbPer, NbPer + tree_data->NodeMax);
  Multiplier = tree_data->FwdRateS[NbPer - N] / tree_data->FwdRate[NbPer - N] *
               exp(tree_data->Delta[NbPer - N]);

  for (i = Bottom; i <= Top; i++) {
    FwdYield = tree_data->FwdRate[NbPer - N] * exp(Jump * (i - NbPer));
    Discount[i] = 1. / (1. + FwdYield);
    if (term_data->TStree !=
        term_data->TSdiscount)          /*include  convexity adjustments */
      FwdYield = FwdYield * Multiplier; /* if discounting is different from
                                           tree-building  , we redo it */

    //#ifdef USECONSISTENTOAS
    //		DiscountS[i] = 1.0/(1.0 + FwdYield)/(1.0 + term_data->OAS * Length /
    //10000.0); #else
    FwdYield +=
        term_data->OAS * Length / 10000.; /* We add the OAS to the current short
                                             term rate to calculate DiscountS */
    DiscountS[i] = 1. / (1. + FwdYield);
    //#endif

    k = (NbPer + 1 > tree_data->NodeMax) * ((i == Bottom) - (i == Top));

    pu[i] = .5 * (k * (k - 1.) + pow(Vol / Jump, 2.) +
                  (1. - 2. * k) * (Drift - Beta * Length * (i - NbPer)) +
                  pow(Drift - Beta * Length * (i - NbPer), 2.));
    pd[i] = .5 * (k * (k + 1.) + pow(Vol / Jump, 2.) -
                  (1. + 2. * k) * (Drift - Beta * Length * (i - NbPer)) +
                  pow(Drift - Beta * Length * (i - NbPer), 2.));
    p0[i] = 1. - pu[i] - pd[i];

  } /* for i */

  return;

} /* Lattice */

void Zero_Calc_new(
    double ***Zero, /* Set of zeros maturing on each coupon date */
    //			int             *NbZero  ,                                /* Current
    //number of zeros being priced */ 			int             TotNbZero  , /* Total
    //number of zeros needed to calculate the par yield index */ 			int Reset  , /*
    //Reset the zero price to 1 if TRUE */
    int NbPer, /* Current time period */
    //			int             TotPer  ,                                 /* Total
    //number of period in the rate tree including the N additional periods */
    int useOAS, MBS_BKTree *tree) {

  double **tmp;
  int lastFwdIndex;
  MBS_BKTree_ZeroPrices(&lastFwdIndex, &tmp, useOAS, NbPer, tree);
  *Zero = tmp;
  return;

} //* Zero_Calc_new

/*****  Zero_Calc  **********************************************************/
/*
 *       Set the zeros to 1 at their respective maturity and then take the
 *       discounted expected value.
 *       Note that if we need TotNbZero to calculate the par yield we will in
 *       fact end up with NbZero = TotNbZero + 1 because Zero[0] is always
 *       equal to 1 on a coupon payment and cannot be used.
 */
void Zero_Calc(double **Zero, /* Set of zeros maturing on each coupon date */
               int *NbZero,   /* Current number of zeros being priced */
               int TotNbZero, /* Total number of zeros needed to calculate the
                                 par yield index */
               int Reset,     /* Reset the zero price to 1 if TRUE */
               double *Discount, /* One period discount factor at each node at
                                    period NbPer */
               double *pu, /* Array of probabilities in the up state in the
                              interest rate tree */
               double *pd, /* Array of probabilities in the down state in the
                              interest rate tree */
               double *p0, /* Array of probabilities in the flat state in the
                              interest rate tree */
               int NbPer,  /* Current time period */
               int TotPer, /* Total number of period in the rate tree including
                              the N additional periods */
               TREE_DATA *tree_data) /* Structure of tree data */
{
  double *TempPointer; /* Temporary pointer to doubles used in the permutation
                          of the zeros */
  int i;

  if (NbPer == TotPer) /* If we are at the last period we just have to set the
                          final zero to 1. */
  {
    Zero_Price(Zero[0], 1, Discount, pu, pd, p0, NbPer, tree_data);
    return;

  } /* if */

  for (i = 0; i < *NbZero;
       i++) /* Discount all the zeros (if we are not at the last period) */
    Zero_Price(
        Zero[i],
        0, /* Reset = 0 because these zeros have been set to 1 earlier. */
        Discount, pu, pd, p0, NbPer, tree_data);

  if (Reset) /* If Reset = TRUE  , we permute the zeros. */
  {
    TempPointer = Zero[*NbZero - 1]; /* Store the address of the zero with the
                                        longest maturity ("rightmost" zero) */

    for (i = *NbZero - 1; i > 0;
         i--) /* Shift all the zeros to the "right" by one. */
      Zero[i] = Zero[i - 1];

    if (*NbZero <= TotNbZero) /* If we don't have all the zeros we need to
                                 calculate the par yield (i.e.   */
    { /* TotNbZero + 1)  , then we need to keep an extra zero: the one which
         address */
      Zero[0] = Zero[*NbZero]; /* is stored in TempPointer. Then we give a new
                                  address to Zero[0] and       */
      Zero[*NbZero] = TempPointer; /* increment the number of zeros by one. */
      (*NbZero)++;
    } /* If we already had enough zeros to calculate the par yield we just have
         to shift */
    else /* the zeros to the right and drop the last one. We did the shift 10
            lines above.  */
      Zero[0] = TempPointer; /* Here we drop the last one by doing
                                Zero[0]=TempPointer=Zero[NbZero-1].          */

    Zero_Price(Zero[0], /* Zero[0] is the zero coupon maturing now. We reset it
                           to 1. */
               1, Discount, pu, pd, p0, NbPer, tree_data);
  } /* if */

  return;

} /* Zero_Calc */

void Zero_Price_new(
    double *Zero, /* Array of zero prices at the current period */
    int Reset,    /* Reset the zero prices to 1 if TRUE */
    int NbPer,    /* Current time period */
    MBS_BKTree *tree) {
  int i; /* Node index */

  if (Reset)
    for (i = 0; i <= 2 * NbPer; i++)
      Zero[i] = 1.;
  else
    MBS_BKTree_BackwardInduct(Zero, 0, NbPer, tree);

  return;

} // Zero_Price_new

/*****  Zero_Price  *********************************************************/
/*
 *       Calculate the zero coupon price in the lattice using the dev routine.
 */
void Zero_Price(double *Zero, /* Array of zero prices at the current period */
                int Reset,    /* Reset the zero prices to 1 if TRUE */
                double *Discount, /* One period discount factor at each node at
                                     period NbPer */
                double *pu, /* Array of probabilities in the up state in the
                               interest rate tree */
                double *pd, /* Array of probabilities in the down state in the
                               interest rate tree */
                double *p0, /* Array of probabilities in the flat state in the
                               interest rate tree */
                int NbPer,  /* Current time period */
                TREE_DATA *tree_data) /* Structure of tree data */
{
  int i; /* Node index */

  if (Reset)
    for (i = 0; i <= 2 * NbPer; i++)
      Zero[i] = 1.;
  else
    Dev(Zero, /* Discounted expected value function */
        Discount, pu, pd, p0, NbPer, tree_data);

  return;

} /* Zero_Price */

/*****  Dev  ****************************************************************/
/*
 *       Discounted expected value.
 */
void Dev(double *Price,    /* Array of prices that has to be discounted */
         double *Discount, /* One period discount factor at each node at period
                              NbPer */
         double *pu, /* Array of probabilities in the up state in the interest
                        rate tree */
         double *pd, /* Array of probabilities in the down state in the interest
                        rate tree */
         double *p0, /* Array of probabilities in the flat state in the interest
                        rate tree */
         int NbPer,  /* Current time period */
         TREE_DATA *tree_data) /* Structure of tree data */
{
  double Temp; /* Temporary variable */
  int Bottom,  /* Bottom node index in the tree */
      Top,     /* Top node index in the tree */
      i,       /* Node index (i.e. vertical) */
      k; /* If k=0 the tree is not cut; if k=1 the tree is cut and branches up
            or down */

  Bottom = (int)max(
      0., NbPer - tree_data->NodeMax); /* If NbPer>NodeMax the tree is cut at
                                          the top and the bottom */
  Top = (int)min(2 * NbPer, NbPer + tree_data->NodeMax);

  Temp = Price[Top - 1]; /* We store the price at node Top -1. */

  for (i = Bottom; i < Top; i++) {
    k = (NbPer + 1 > tree_data->NodeMax) * ((i == Bottom) - (i == Top));

    Price[i] = pu[i] * Price[i + 2 + k] + p0[i] * Price[i + 1 + k] +
               pd[i] * Price[i + k];
    Price[i] *= Discount[i];

  } /* for i */
  /* i=Top is a special case because if the tree is cut we modified Price[Top-1]
   */
  k = -(NbPer + 1 >
        tree_data->NodeMax); /* in the previous loop and we need the old value
                                (Temp) to calculate Price[Top] */

  Price[Top] = pu[Top] * Price[Top + 2 + k] + p0[Top] * Price[Top + 1 + k] +
               pd[Top] * (Price[Top] * (1 + k) - Temp * k);
  Price[Top] *= Discount[Top];

  return;

} /* Dev */

/*****  Par_Yield  **********************************************************/
/*
 *       Calculate an array of par yields for a specific maturity.
 */
void Par_Yield(
    double
        *ParYield, /* Output: par yield at every node at the current period */
    double **Zero, /* Set of zeros maturing on each coupon date */
    int TotNbZero, /* Total number of zeros needed to calculate the par yield
                      index */
    int F,         /* Frequency of underlying payment as a integer */
    int IndexF,    /* Frequency of index payment as a integer */
    int NodeMax,   /* Maximum node reached in the tree (both top and bottom) */
    int NbPer)     /* Current time period */
{
  double Annuity = 0.; /* Price of an annuity paying $1 at every coupon date */
  int Bottom,          /* Bottom node index in the tree */
      Top,             /* Top node index in the tree */
      i, j;

  Bottom = (int)max(0., NbPer - NodeMax); /* If NbPer>NodeMax the tree is cut at
                                             the top and the bottom */
  Top = (int)min(2 * NbPer, NbPer + NodeMax);

  for (i = Bottom; i <= Top; i++) {
    for (j = F / IndexF; j <= TotNbZero;
         j += F / IndexF) /* Because the SFS and the index might have different
                             maturity we don't use all  */
      Annuity += Zero[j][i]; /* the zeros only zeros at F/IndexF intervals. Note
                                also that the index stops at */
                             /* TotNbZero not TotNbZero-1.                          */
    ParYield[i] = 100. * IndexF * (1. - Zero[TotNbZero][i]) / Annuity;

    Annuity = 0.; /* Reset Annuity to zero for next iteration */

  } /* for i */

  return;
} /* Par_Yield */

char *copyTermStruct(MBS_BKTree *tree, MBSPT_DealStruct *deal_struct,
                     TREE_DATA *tree_data) {
  int i, j, sz;
  MBSPT_CashFlowStruct *cashflow_struct = &(deal_struct->cashflow_struct);
  // translate time division:
  i = 0;
  while (tree->tree_times[i] < cashflow_struct->duePeriodEndTimes[0])
    ++i;
  tree_data->ExerPosition[0] = i;
  while (tree->tree_times[i] <
         cashflow_struct->duePeriodEndTimes[cashflow_struct->numCashFlows - 1])
    ++i;
  tree_data->ExerPosition[1] = i;
  tree_data->NbNodes =
      tree->NbTimeSlices -
      2; // without 0 <-> 0.0 AND tree_data->Length is [0  , NbNodes]
  Tree_Alloc(tree_data, tree_data->NbNodes);
  tree_data->ExEnd = tree_data->ExerPosition[1];

  for (i = 0; i <= tree_data->NbNodes; ++i)
    tree_data->Length[i] = tree->tree_times[i + 1] - tree->tree_times[i];

  tree_data->Accrued[0] =
      cashflow_struct->accrualFactor /
      (cashflow_struct->accrualFactor + tree->tree_times[1]);
  j = 1;
  for (i = 0; i < cashflow_struct->numCashFlows; ++i) {
    while (tree->tree_times[j] < cashflow_struct->duePeriodEndTimes[i])
      ++j;
    tree_data->Accrued[j] = 1.0;
  }
  for (i = j; i <= tree_data->NbNodes; ++i)
    tree_data->Accrued[i] = 1.0; // these are needed because it affects how
                                 // zeros are computed in backward induction
  // translate forward
  sz = tree_data->NbNodes + 1;

  for (i = 0; i < sz; ++i) {
    tree_data->ZeroCoupon[i] = tree->DFS[i];
    tree_data->FwdRate[i] = tree->DFS[i] / tree->DFS[i + 1] - 1.0;
  }
  tree_data->ZeroCoupon[sz] = tree->DFS[sz];

  // fill in fwd_vol
  tree_data->Jump = tree->DX;
  for (i = 0; i <= tree_data->NbNodes; ++i) {
    tree_data->FwdVol[i] = tree->sigmas[i];
  }
  ///
  for (i = 0; i < tree->NbTimeSlices; ++i)
    tree_data->Drift[i] = tree->Drift[i];
  tree_data->NodeMax = tree->NodeMax;
  //
  return (0);
}