#include "LGM2FRefInstrs.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_lgmeval.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"

#define TREE_MESH_SPACING                                                      \
  1.73205080756888 /* 1.73205080756888 theoretically > 4  (1.414213562373) */
#define TREE_MESH_SPACING_SQUARE 3.0
#define ONE_HALF 0.5
#define TREE_LIM_IN_STDEV 5

#define ALLOC_CONV_TREE1                                                       \
  {                                                                            \
    zeta11Deal = dvector(0, NumTimeStep);                                      \
    zeta12Deal = dvector(0, NumTimeStep);                                      \
    zeta22Deal = dvector(0, NumTimeStep);                                      \
  }

#define FREE_CONV_TREE1                                                        \
  {                                                                            \
    free_dvector(zeta11Deal, 0, NumTimeStep);                                  \
    free_dvector(zeta12Deal, 0, NumTimeStep);                                  \
    free_dvector(zeta22Deal, 0, NumTimeStep);                                  \
    srt_free(tDeal);                                                           \
  }

#define ALLOC_CONV_TREE2                                                       \
  {                                                                            \
    max_index0 = lngvector(0, NumTimeStep);                                    \
    max_index1 = lngvector(0, NumTimeStep);                                    \
    dual_basis00 = dvector(0, NumTimeStep);                                    \
    dual_basis01 = dvector(0, NumTimeStep);                                    \
    dual_basis10 = dvector(0, NumTimeStep);                                    \
    dual_basis11 = dvector(0, NumTimeStep);                                    \
    inverse_basis00 = dvector(0, NumTimeStep);                                 \
    inverse_basis01 = dvector(0, NumTimeStep);                                 \
    inverse_basis10 = dvector(0, NumTimeStep);                                 \
    inverse_basis11 = dvector(0, NumTimeStep);                                 \
    spacing0 = dvector(0, NumTimeStep);                                        \
    spacing1 = dvector(0, NumTimeStep);                                        \
    IsVol1Null = lngvector(0, NumTimeStep);                                    \
    IsVol2Null = lngvector(0, NumTimeStep);                                    \
  }

#define FREE_CONV_TREE2                                                        \
  {                                                                            \
    free_lngvector(max_index0, 0, NumTimeStep);                                \
    free_lngvector(max_index1, 0, NumTimeStep);                                \
    free_dvector(dual_basis00, 0, NumTimeStep);                                \
    free_dvector(dual_basis01, 0, NumTimeStep);                                \
    free_dvector(dual_basis10, 0, NumTimeStep);                                \
    free_dvector(dual_basis11, 0, NumTimeStep);                                \
    free_dvector(inverse_basis00, 0, NumTimeStep);                             \
    free_dvector(inverse_basis01, 0, NumTimeStep);                             \
    free_dvector(inverse_basis10, 0, NumTimeStep);                             \
    free_dvector(inverse_basis11, 0, NumTimeStep);                             \
    free_dvector(spacing0, 0, NumTimeStep);                                    \
    free_dvector(spacing1, 0, NumTimeStep);                                    \
    free_lngvector(IsVol1Null, 0, NumTimeStep);                                \
    free_lngvector(IsVol2Null, 0, NumTimeStep);                                \
  }

/* MID-AT PAYOFF AT THE TIME STEP lj0 */

LGMErr
Midat2DPayoff(double *payoff, double *swap, double X, double Y, void *dealPtr,
              Date tNow, long jEx, String ycName, double *Zeta1, double *Zeta2,
              double *Zeta12, double *H1, double *H2, double gamma,
              LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double *),
              LGMErr (*GetBeta)(Date, Date, double *)) {
  SrtSimMidAtPtr MidAt;
  double Zeta1tEx, Zeta2tEx, Zeta12tEx, H1tStart, H2tStart, H1tPay, H2tPay;
  double ZCtStart, ZCtPay, SwapLevel, Strike;
  Date tStart, *MidAttPay, tPay;
  double *MidAtCvgPay, MidAtCvgFirst, ak;
  double DtStart, DtPay;
  double CallPut;
  long MidAtnPay, MidAtiFirst, MidAtnCpn, k, nEx;
  LGMErr err = NULL;

  MidAt = (SrtSimMidAt *)dealPtr;
  CallPut = (MidAt->PayRec == SRT_PAYER) ? -1.0 : 1.0;

  Zeta1tEx = Zeta1[jEx];
  Zeta2tEx = Zeta2[jEx];
  Zeta12tEx = Zeta12[jEx];

  nEx = MidAt->nEx;
  MidAtnPay = MidAt->nPay;
  MidAttPay = MidAt->tPay;
  MidAtCvgPay = MidAt->MidAtCvgPay;
  MidAtCvgFirst = MidAt->MidAtCvgFirst[jEx];

  MidAtiFirst = MidAt->FirstPay[jEx];
  tStart = MidAt->tStart[jEx];
  Strike = MidAt->Strike[jEx];

  H1tStart = H1[jEx];
  H2tStart = H2[jEx];

  DtStart = swp_f_df(tNow, tStart, ycName);
  if (DtStart == SRT_DF_ERROR)
    return (err);

  ZCtStart = DtStart *
             exp(H1tStart * X + H2tStart * Y - H1tStart * H2tStart * Zeta12tEx -
                 0.5 * (H1tStart * H1tStart * Zeta1tEx +
                        H2tStart * H2tStart * Zeta2tEx));

  (*payoff) = -Strike * ZCtStart;

  MidAtnCpn = (MidAtnPay - MidAtiFirst - 1);

  SwapLevel = 0.0;
  for (k = 0; k <= MidAtnCpn; k++) {
    tPay = MidAttPay[MidAtiFirst + k];
    DtPay = swp_f_df(tNow, tPay, ycName);
    if (DtPay == SRT_DF_ERROR)
      return (err);

    H1tPay = H1_Func(tPay, MidAt->HDates, tNow, nEx, H1);
    H2tPay = H2_Func(tPay, MidAt->HDates, tNow, nEx, H1, gamma);

    ZCtPay =
        DtPay *
        exp(H1tPay * X + H2tPay * Y - H1tPay * H2tPay * Zeta12tEx -
            0.5 * (H1tPay * H1tPay * Zeta1tEx + H2tPay * H2tPay * Zeta2tEx));

    ak = MidAt->Payment[MidAtiFirst + k] - MidAt->RedFirstPay[jEx];

    (*payoff) += ak * ZCtPay;

    /* Computation of the SwapLevel */
    if ((k == 0) && (MidAtnCpn > 0))
      SwapLevel += ZCtPay * MidAtCvgFirst;
    else if ((k > 0) && (k < MidAtnCpn))
      SwapLevel += ZCtPay * MidAtCvgPay[MidAtiFirst + k - 1];
    else if (k == MidAtnCpn) {
      if (MidAtnCpn == 0)
        SwapLevel += ZCtPay * MidAtCvgFirst;
      else
        SwapLevel += ZCtPay * MidAtCvgPay[MidAtiFirst + k - 1];
    }

  } /* END OF LOOP ON k */

  (*swap) = (ZCtStart - ZCtPay) / SwapLevel;
  (*payoff) = max(CallPut * (*payoff), 0.0);

  return (err);
}

static void
LinearInterpTS(Date tNow, long nEx, Date *tEx,
               double *zeta11, /* [0        , 1        , ...        , nEx-1] ,
                            values of zeta11 at the exercise dates */
               double *zeta12, /* [0        , 1        , ...        , nEx-1] ,
                            values of zeta12 at the exercise dates */
               double *zeta22, /* [0        , 1        , ...        , nEx-1] ,
                            values of zeta22 at the exercise dates */
               long j, Date *tDeal,
               double **zeta11Deal, /* [0        , 1        , ...        , j] ,
                                 values of zeta_kl the deal dates */
               double **zeta12Deal, /*  zeta_kl(0)=0.0 */
               double **zeta22Deal) {
  double Weight;
  long k, index;
  LGMErr error = NULL;

  k = 1;
  while ((k < j) && (tDeal[k] < tEx[0])) {
    Weight = ((double)(tDeal[k] - tNow)) / (tEx[0] - tNow);
    (*zeta11Deal)[k] = zeta11[0] * Weight;
    (*zeta22Deal)[k] = zeta22[0] * Weight;
    (*zeta12Deal)[k] = zeta12[0] * Weight;

    k++;
  }

  (*zeta11Deal)[k] = zeta11[0];
  (*zeta22Deal)[k] = zeta22[0];
  (*zeta12Deal)[k] = zeta12[0];

  index = 1;
  while ((k < j)) {
    if (tDeal[k] == tEx[index]) {
      (*zeta11Deal)[k] = zeta11[index];
      (*zeta12Deal)[k] = zeta12[index];
      (*zeta22Deal)[k] = zeta22[index];
    } else {
      while (tDeal[k] > tEx[index])
        index++; /* tDeal[k] <= tEx[index] */

      Weight =
          ((double)(tDeal[k] - tEx[index - 1])) / (tEx[index] - tEx[index - 1]);
      (*zeta11Deal)[k] =
          zeta11[index - 1] + Weight * (zeta11[index] - zeta11[index - 1]);
      (*zeta12Deal)[k] =
          zeta12[index - 1] + Weight * (zeta12[index] - zeta12[index - 1]);
      (*zeta22Deal)[k] =
          zeta22[index - 1] + Weight * (zeta22[index] - zeta22[index - 1]);
    }

    k++;
  }

  (*zeta11Deal)[j] = zeta11[nEx - 1]; /* zeta11(tDeal[j]=tEx[nEx]) */
  (*zeta22Deal)[j] = zeta22[nEx - 1];
  (*zeta12Deal)[j] = zeta12[nEx - 1];
}

LGMErr LGMAutoCal2DTree(
    /* info about convolutions */
    long nEx,                  /* number of exercises */
    Date *tEx, double *zeta11, /* [0        , 1        , ...        , nEx-1] ,
                      values of zeta11 at the exercise dates */
    double *zeta12, /* [0        , 1        , ...        , nEx-1]        ,
                 values of zeta12 at the exercise dates */
    double *zeta22, /* [0        , 1        , ...        , nEx-1]        ,
                 values of zeta22 at the exercise dates */
    double *G1, double *G2, double gamma,
    /* info about today's discount curve and swaption/cap vols */
    Date EvalDate, /* tNow for deal evaluation */
    String ycName, /* yield curve name */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *),              /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* swaption/cap exponents (beta) */
                                             /* information about the deal */
    void *dealPtr, LGMErr (*payofffunc)(),
    /* output */
    double *answer,   /* value of the deal */
    double **xExBdry) /* array of exercise points (in x) */
{
  /* declarations */
  LGMErr error = NULL;

  double *dual_basis00, *dual_basis01, *dual_basis10, *dual_basis11,
      *inverse_basis00, *inverse_basis11, *inverse_basis01, *inverse_basis10,
      *spacing0, *spacing1, payoff, **TheMatCur, **TheMatPrev, **TheMatTemp;
  /* TS at tDeal [0        ,1        ,...        ,NumTimeStep-1] */
  double *zeta11Deal, *zeta12Deal, *zeta22Deal;

  double p[5];
  long *IsVol1Null, *IsVol2Null, *max_index0, *max_index1;
  Date *tDeal;

  double local_var0, local_var1, local_cov, local_correl, V1, V2, C, temp1,
      temp2, temp3, temp4, temp5, temp6, temp7, det, global_var_dir0,
      global_var_dir1, temp_global_var_dir0, temp_global_var_dir1, std, lim,
      StateVar0, StateVar1, forward0, forward1, vX0, vX1, dx, dxx, dy, dyy, xx,
      StatevarMid0, StatevarMid1, StatevarUp0, StatevarUp1, StatevarDown0,
      StatevarDown1;

  double TimeStep = 1.0 / 24, dt;
  long MinNode = 350, MaxNode = 1500, n;

  double Start, End;
  long NumInterval;
  long max_max_index0, max_max_index1,
      IndexEx = 0, IndexTemp = 0, IndexCurrent = 0, i, j, k, l, max_trunc_index,
      IndexMid0, IndexMid1, IndexUp0, IndexUp1, IndexDown0, IndexDown1,
      NumTimeStep;
  double swap, l_stop;
  SRT_Boolean exercise;

  SrtSimMidAtPtr Midat;
  SrtReceiverType PayRec;
  long min_l_index, max_l_index;
  double *sum_prob_vX0, prob_vX0, X0;

  Midat = (SrtSimMidAt *)dealPtr;
  PayRec = Midat->PayRec;

  /* ---------------------- STEP 1: INITIALISATION -------------------------- */
  /* Eliminate trivial case */
  if (nEx < 1) {
    *answer = 0.0;
    return (NULL);
  }
  if (xExBdry) {
    *xExBdry = (double *)srt_calloc(nEx, sizeof(double));
    if (*xExBdry == NULL)
      return ("allocation failed in treeautocal2d");
  }

  /** FIRST STEP INITIALISATION OF THE DATES **/
  TimeStep = max(1.0 / DAYS_IN_YEAR,
                 min((double)(tEx[nEx - 1] - EvalDate) / DAYS_IN_YEAR / MinNode,
                     TimeStep));
  tDeal = (Date *)malloc((MaxNode) * sizeof(Date));

  sum_prob_vX0 = dvector(0, nEx - 1);

  n = 0;
  tDeal[0] = EvalDate;
  Start = (double)EvalDate;
  for (i = 0; i < nEx; i++) {
    End = (double)tEx[i];
    NumInterval = ((long)((End - Start) / DAYS_IN_YEAR / TimeStep)) + 1;
    NumInterval = min(NumInterval, (long)(End - Start));
    dt = (End - Start) / DAYS_IN_YEAR / NumInterval;
    for (j = 1; j < NumInterval; j++) {
      n++;
      tDeal[n] = (Date)((double)tDeal[n - 1] + dt * DAYS_IN_YEAR);
    }
    n++;
    tDeal[n] = tEx[i];
    Start = End;
  }
  NumTimeStep = n;

  /* Initialise the TS at the tDeal[0        ,1        ,...        ,NumTimeStep]
   and Interpolate Linearly zetaDeal[0        ,1        ,... ,NumTimeStep-1]*/
  ALLOC_CONV_TREE1

  if (zeta11Deal == NULL || zeta12Deal == NULL || zeta22Deal == NULL) {
    smessage("Warning: Allocation failed 2D Conv-Tree");
    FREE_CONV_TREE1
    return "Allocation Failed";
  }

  LinearInterpTS(EvalDate, nEx, tEx, zeta11, zeta12, zeta22, NumTimeStep, tDeal,
                 &zeta11Deal, &zeta12Deal, &zeta22Deal);

  /* SECOND STEP THE NEW BASIS */

  /* allocation of all following pointers */
  ALLOC_CONV_TREE2

  if (max_index0 == NULL || max_index1 == NULL || dual_basis00 == NULL ||
      dual_basis01 == NULL || dual_basis10 == NULL || dual_basis11 == NULL ||
      inverse_basis00 == NULL || inverse_basis01 == NULL ||
      inverse_basis10 == NULL || inverse_basis11 == NULL || spacing0 == NULL ||
      spacing1 == NULL || IsVol1Null == NULL || IsVol2Null == NULL) {
    smessage("Warning: Allocation failed 2D Conv-Tree");
    FREE_CONV_TREE1
    FREE_CONV_TREE2
    return "Allocation Failed";
  }

  /* No change of basis on first step: Transfer matrices are equal to 1 */
  dual_basis00[0] = dual_basis11[0] = inverse_basis00[0] = inverse_basis11[0] =
      1;
  dual_basis10[0] = dual_basis01[0] = inverse_basis01[0] = inverse_basis10[0] =
      0;
  ;

  /* Init maximum indexation */
  max_max_index0 = max_max_index1 = 0;

  for (j = 1; j <= NumTimeStep; j++) {
    /* Local */
    local_var0 = zeta11Deal[j] - zeta11Deal[j - 1];
    local_var1 = zeta22Deal[j] - zeta22Deal[j - 1];
    local_cov = zeta12Deal[j] - zeta12Deal[j - 1];

    /* Check for numerical errors */
    local_correl = local_cov / sqrt(local_var0 * local_var1);
    if (local_correl > 1.0 - EPS) {
      smessage("Warning: local correl > 1 in 2D Conv-Tree");
      local_correl = 1.0 - EPS;
      local_cov = local_correl * sqrt(local_var0 * local_var1);
    } else if (local_correl < -1.0 + EPS) {
      smessage("Warning: local correl < -1 in 2D Conv-Tree");
      local_correl = -1.0 + EPS;
      local_cov = local_correl * sqrt(local_var0 * local_var1);
    }

    /* Diagonalise Local var matrix        , compute sqrt(cov) */
    V1 = local_var0;
    V2 = local_var1;
    C = local_cov;
    /* Covariance is significative: need to diagonalise and multiply by
     * eigenvalues */
    if (fabs(C) > EPS) {
      temp1 = sqrt(4 * C * C + (V2 - V1) * (V2 - V1));
      temp2 = sqrt(1.0 + pow(V2 - V1 + temp1, 2) / (4 * C * C));
      temp3 = 2 * sqrt(2) * C * temp2;
      temp4 = sqrt(2) * temp2;
      temp5 = sqrt(1.0 + pow(V2 - V1 - temp1, 2) / (4 * C * C));
      temp6 = 2 * sqrt(2) * C * temp5;
      temp7 = sqrt(2) * temp5;
      dual_basis00[j] = sqrt(V1 + V2 - temp1) * (V1 - V2 - temp1) / temp3;
      dual_basis10[j] = sqrt(V1 + V2 - temp1) / temp4;
      dual_basis01[j] = sqrt(V1 + V2 + temp1) * (V1 - V2 + temp1) / temp6;
      dual_basis11[j] = sqrt(V1 + V2 + temp1) / temp7;
    } else
    /*	covariance is already diagonal */
    {
      dual_basis00[j] = sqrt(V1);
      dual_basis01[j] = dual_basis10[j] = 0.0;
      dual_basis11[j] = sqrt(V2);
    }

    /* DUAL BASIS : COORDINATES OF EIGENVECTORS IN TERMS OF X1        , X2
       MULTIPLIED BY EIGENVALUES DUAL BASIS  = SQRT(D)*P IN THE DECOMPOSITION
       PD(tP) */
    /* compute Inverse basis
       = coordinates of x1 and x2 in terms of eigenvectors */
    det = dual_basis00[j] * dual_basis11[j] - dual_basis01[j] * dual_basis10[j];
    inverse_basis00[j] = dual_basis11[j] / det;
    inverse_basis10[j] = -dual_basis10[j] / det;
    inverse_basis01[j] = -dual_basis01[j] / det;
    inverse_basis11[j] = dual_basis00[j] / det;

    /* Trimming */
    /* Compute global variance in the direction of local eigenvectors */
    temp_global_var_dir0 =
        zeta11Deal[j] * inverse_basis00[j] * inverse_basis00[j] +
        zeta22Deal[j] * inverse_basis01[j] * inverse_basis01[j] +
        2.0 * zeta12Deal[j] * inverse_basis00[j] * inverse_basis01[j];

    /* Compute global variance in the direction of local eigenvectors */
    temp_global_var_dir1 =
        zeta11Deal[j] * inverse_basis10[j] * inverse_basis10[j] +
        zeta22Deal[j] * inverse_basis11[j] * inverse_basis11[j] +
        2.0 * zeta12Deal[j] * inverse_basis10[j] * inverse_basis11[j];

    /* Test if the volatility in one direction is Null */
    if (V1 <= 1e-8 || V2 <= 1e-8) {
      if (fabs(temp_global_var_dir0 - global_var_dir0) < 1e-8)
        IsVol1Null[j] = 1;
      if (fabs(temp_global_var_dir1 - global_var_dir1) < 1e-8)
        IsVol2Null[j] = 1;
    }

    global_var_dir0 = temp_global_var_dir0;
    global_var_dir1 = temp_global_var_dir1;
    /* Compute Spacing */
    spacing0[j] = spacing1[j] = TREE_MESH_SPACING;

    /* Compute standard deviations and limits (TREE_LIM_IN_STDEV * std) */
    std = sqrt(global_var_dir0);
    lim = TREE_LIM_IN_STDEV * std;
    /* Compute max index */
    max_trunc_index = (long)DTOL(lim / spacing0[j]) + 1;

    /* Store the information at the time step tree info */
    max_index0[j] = max_trunc_index;

    /* Update maximum indexation ever met */
    max_max_index0 = IMAX(max_max_index0, max_trunc_index);

    /* Compute standard deviations and limits (TREE_LIM_IN_STDEV * std) */
    std = sqrt(global_var_dir1);
    lim = TREE_LIM_IN_STDEV * std;
    /* Compute max index */
    max_trunc_index = (int)DTOL(lim / spacing1[j]) + 1;

    /* Store the information at the time step tree info */
    max_index1[j] = max_trunc_index;

    /* Update maximum indexation ever met */
    max_max_index1 = IMAX(max_max_index1, max_trunc_index);

    /* Now we can allocate the matrix at time j */
  }

  /*---------------- THIRD STEP : BACKWARD INDUCTION -------------------*/

  TheMatCur =
      dmatrix(-max_max_index0, max_max_index0, -max_max_index1, max_max_index1);
  TheMatPrev =
      dmatrix(-max_max_index0, max_max_index0, -max_max_index1, max_max_index1);

  /** At the maturity date: initialisation of the Grid... **/
  /** (vX0        ,vX1)' IS A TWO DIM. GUAUSSIAN VECT.- P*SQRT(D)*(vX0 ,vX1)'
   * HAS A VAR EQUAL TO  */
  for (k = -max_index0[NumTimeStep]; k <= max_index0[NumTimeStep]; k++)
    for (l = -max_index1[NumTimeStep]; l <= max_index1[NumTimeStep]; l++) {
      vX0 = spacing0[NumTimeStep] * k;
      vX1 = spacing1[NumTimeStep] * l;
      StateVar0 =
          dual_basis00[NumTimeStep] * vX0 + dual_basis01[NumTimeStep] * vX1;
      StateVar1 =
          dual_basis10[NumTimeStep] * vX0 + dual_basis11[NumTimeStep] * vX1;
      /* Compute payoff for this date j */
      error = (*payofffunc)(&payoff, &swap, StateVar0, StateVar1, dealPtr,
                            EvalDate, nEx - 1, ycName, zeta11, zeta22, zeta12,
                            G1, G2, gamma, GetVol, GetBeta);
      if (error)
        return error;
      TheMatPrev[k][l] = payoff;
    }
  /* last exercise boundary */

  if (xExBdry)
    (*xExBdry)[nEx - 1] = ((SrtSimMidAt *)dealPtr)->MidAtStrike[nEx - 1];

  /** The main LOOP   **/
  IndexEx = nEx - 2; /* exercice date before the last exercice */
  for (j = NumTimeStep - 1; j >= 0; j--) {
    swap = 0.0;
    for (k = -max_index0[j]; k <= max_index0[j]; k++) /* X[k] */
    {
      exercise = SRT_TRUE;
      if (PayRec == SRT_RECEIVER) /* REVERSE THE ORDER */
      {
        min_l_index = max_index1[j];
        max_l_index = -max_index1[j];
      } else {
        min_l_index = -max_index1[j];
        max_l_index = max_index1[j];
      }

      l = min_l_index;
      l_stop = 0;
      while (l_stop == 0) {
        /* Computes neighbours points in dual basis */
        vX0 = spacing0[j] * k;
        vX1 = spacing1[j] * l;
        StateVar0 = dual_basis00[j] * vX0 + dual_basis01[j] * vX1;
        StateVar1 = dual_basis10[j] * vX0 + dual_basis11[j] * vX1;
        forward0 = inverse_basis00[j + 1] * StateVar0 +
                   inverse_basis01[j + 1] * StateVar1;
        forward1 = inverse_basis10[j + 1] * StateVar0 +
                   inverse_basis11[j + 1] * StateVar1;

        xx = (forward0 < 0.0) ? -0.5 : 0.5;
        IndexMid0 = (int)(forward0 / spacing0[j + 1] + xx);
        StatevarMid0 = IndexMid0 * spacing0[j + 1];
        IndexMid0 = DMAX(min(IndexMid0, max_index0[j + 1]), -max_index0[j + 1]);

        IndexUp0 =
            (IndexMid0 < max_index0[j + 1]) ? IndexMid0 + 1 : max_index0[j + 1];
        StatevarUp0 = StatevarMid0 + spacing0[j + 1];

        IndexDown0 = (IndexMid0 > -max_index0[j + 1]) ? IndexMid0 - 1
                                                      : -max_index0[j + 1];
        StatevarDown0 = StatevarMid0 - spacing0[j + 1];

        xx = (forward1 < 0.0) ? -0.5 : 0.5;
        IndexMid1 = (int)(forward1 / spacing1[j + 1] + xx);
        StatevarMid1 = IndexMid1 * spacing1[j + 1];
        IndexMid1 = DMAX(min(IndexMid1, max_index1[j + 1]), -max_index1[j + 1]);

        IndexUp1 =
            (IndexMid1 < max_index1[j + 1]) ? IndexMid1 + 1 : max_index1[j + 1];
        StatevarUp1 = StatevarMid1 + spacing1[j + 1];

        IndexDown1 = (IndexMid1 > -max_index1[j + 1]) ? IndexMid1 - 1
                                                      : -max_index1[j + 1];
        StatevarDown1 = StatevarMid1 - spacing1[j + 1];

        /* Computes the 5 connection probabilities */
        /* Points are like this : 	0: center ; 1 : right (x1 Up/x2 Mid) ; 2
         * : up ; 3 : left ; 4 : down */
        dx = forward0 - StatevarMid0;
        dy = forward1 - StatevarMid1;
        dxx = 1.0 + dx * dx;
        dyy = 1.0 + dy * dy;
        p[1] = ONE_HALF * (dxx + dx * TREE_MESH_SPACING) /
               TREE_MESH_SPACING_SQUARE;
        p[3] = ONE_HALF * (dxx - dx * TREE_MESH_SPACING) /
               TREE_MESH_SPACING_SQUARE;
        p[2] = ONE_HALF * (dyy + dy * TREE_MESH_SPACING) /
               TREE_MESH_SPACING_SQUARE;
        p[4] = ONE_HALF * (dyy - dy * TREE_MESH_SPACING) /
               TREE_MESH_SPACING_SQUARE;

        if (IsVol2Null[j])
          p[2] = p[4] = 0.0;
        if (IsVol1Null[j])
          p[1] = p[3] = 0.0;

        p[0] = 1.0 - p[1] - p[2] - p[3] - p[4];

        /* conditional expectation */
        TheMatCur[k][l] = p[0] * TheMatPrev[IndexMid0][IndexMid1] +
                          p[1] * TheMatPrev[IndexUp0][IndexMid1] +
                          p[2] * TheMatPrev[IndexMid0][IndexUp1] +
                          p[3] * TheMatPrev[IndexDown0][IndexMid1] +
                          p[4] * TheMatPrev[IndexMid0][IndexDown1];

        /* At exercise dates */
        if ((tDeal[j] == tEx[IndexEx]) && (IndexEx >= 0)) {
          error = (*payofffunc)(&payoff, &swap, StateVar0, StateVar1, dealPtr,
                                EvalDate, IndexEx, ycName, zeta11, zeta22,
                                zeta12, G1, G2, gamma, GetVol, GetBeta);

          if (error)
            return error;

          /* StateVar0 =	sqrt(Zeta11).U         , U : N(0        ,1) */

          if ((payoff > TheMatCur[k][l])) {
            TheMatCur[k][l] = payoff;

            if ((exercise == SRT_TRUE) &&
                (xExBdry)) /* i.e 1( Y == f(X) ) == 1 */
            {
              X0 = inverse_basis00[j] * StateVar0 +
                   inverse_basis01[j] * StateVar1;

              temp_global_var_dir0 =
                  zeta11Deal[j] * inverse_basis00[j] * inverse_basis00[j] +
                  zeta22Deal[j] * inverse_basis01[j] * inverse_basis01[j] +
                  2.0 * zeta12Deal[j] * inverse_basis00[j] * inverse_basis01[j];

              prob_vX0 =
                  gauss(X0 / sqrt(temp_global_var_dir0)) / temp_global_var_dir0;
              sum_prob_vX0[IndexEx] += prob_vX0;

              (*xExBdry)[IndexEx] += prob_vX0 * swap;
              exercise = SRT_FALSE;
            }
          }
          /*	COMPUTE E(Z) WHERE Z = E( S(X        ,Y).1( Y = f(X) ) | X ) */
          /* Y = f(X) DEFINES THE EXER. FRONTIER */
        }

        if (PayRec == SRT_RECEIVER) {
          if (l == max_l_index)
            l_stop = 1.0;
          if (l_stop == 0)
            l--;
        } else {
          if (l == max_l_index)
            l_stop = 1.0;
          if (l_stop == 0)
            l++;
        }

      } /* END OF LOOP ON l */

    } /* END OF LOOP ON k */

    if ((tDeal[j] == tEx[IndexEx]) && (IndexEx >= 0)) {
      if ((sum_prob_vX0[IndexEx] > 0) && (xExBdry))
        (*xExBdry)[IndexEx] /= sum_prob_vX0[IndexEx];
      IndexEx--;
    }

    TheMatTemp = TheMatPrev;
    TheMatPrev = TheMatCur;
    TheMatCur = TheMatTemp;
  } /* End loop on time steps */

  (*answer) = TheMatPrev[0][0];

  /* Main Grid UnAllocation */
  free_dmatrix(TheMatCur, -max_max_index0, max_max_index0, -max_max_index1,
               max_max_index1);
  free_dmatrix(TheMatPrev, -max_max_index0, max_max_index0, -max_max_index1,
               max_max_index1);
  if (sum_prob_vX0)
    free_dvector(sum_prob_vX0, 0, nEx - 1);
  sum_prob_vX0 = NULL;

  FREE_CONV_TREE2
  FREE_CONV_TREE1

  return error;
}

/*---------------------------------------------------------------------*/

LGMErr LGMAutoCal2DVega(
    int nScenarii,              /* number of scenarii different than spot */
    long nEx,                   /* number of exercises */
    Date *tEx, double **zeta11, /* for each scenario        , values of */
    double **zeta12,            /* zeta11 at the exercise dates */
    double **zeta22, double **G1, double **G2, double gamma,
    /* info about today's discount curve and swaption/cap vols */
    Date EvalDate, /* tNow for deal evaluation */
    String ycName, /* yield curve name */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean,
                     double *),              /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double *), /* swaption/cap exponents (beta) */
                                             /* information about the deal */
    void *dealPtr, LGMErr (*payofffunc)(),
    /* output */
    double *answer,   /* value of the deal */
    double **xExBdry) /* array of exercise points (in x) */
{
  LGMErr error = NULL;
  int i;
  char buf[200];

  for (i = 1; i <= nScenarii; i++) {
    error =
        LGMAutoCal2DTree(nEx, tEx, zeta11[i], zeta12[i], zeta22[i], G1[i],
                         G2[i], gamma, EvalDate, ycName, GetVol, GetBeta,
                         (void *)dealPtr, Midat2DPayoff, &(answer[i]), xExBdry);
    if (error)
      return error;
    sprintf(buf, "%lf\n", answer[i]);
    smessage(buf);
  }

  return error;
}

#undef ALLOC_CONV_TREE1
#undef FREE_CONV_TREE1
#undef ALLOC_CONV_TREE2
#undef FREE_CONV_TREE2
#undef TREE_MESH_SPACING
#undef TREE_MESH_SPACING_SQUARE
#undef ONE_HALF
#undef TREE_LIM_IN_STDEV

/* ========= END OF FILE =================================================== */
