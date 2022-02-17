/* ==========================================================================================

   FILENAME:       srt_f_rrf_accrual.c

   PURPOSE:        computation of the price of an Accrual Resettable Range
   Floater in a Black-Scholes like framework (with input of Forward Volatility )
   ==========================================================================================
 */
#include "math.h"
#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <num_h_allhdr.h"
#include <srt_h_rrf_accrual.h"

#define Tolerance 1e-6
#define TINY 1e-20

/* A useful macro to free all the used memory */

#define FREE_ALL_ACCRUAL_RRF_MEMORY                                            \
  {                                                                            \
    if (HermiteWeights)                                                        \
      free_dvector(HermiteWeights, 1, NumHermitePoints);                       \
    if (HermiteAbscissas)                                                      \
      free_dvector(HermiteAbscissas, 1, NumHermitePoints);                     \
    if (SlopeStdT0T1)                                                          \
      free_dvector(SlopeStdT0T1, 2, (s_N + 1));                                \
    if (OptStrike)                                                             \
      free_dmatrix(OptStrike, 1, NumHermitePoints, 1, NumHermitePoints);       \
    if (s_Theta)                                                               \
      free_dvector(s_Theta, 2, (s_N + 1));                                     \
    if (s_SlopeVarT0T1)                                                        \
      free_dvector(s_SlopeVarT0T1, 2, (s_N + 1));                              \
    if (s_FwdVol)                                                              \
      free_dvector(s_FwdVol, 2, (s_N + 1));                                    \
    if (s_FwdStDev)                                                            \
      free_dvector(s_FwdStDev, 2, (s_N + 1));                                  \
    if (s_SlopeT0)                                                             \
      free_dvector(s_SlopeT0, 2, (s_N + 1));                                   \
  }

/* ------------------------
 * ----------------------------------------------------------- */

/* Static variables required for the numerical optimisation of the bandwidth */

static double *s_Theta = NULL;
static double *s_FwdVol = NULL;
static double *s_FwdStDev = NULL;
static double *s_SlopeT0 = NULL;
static double *s_SlopeVarT0T1 = NULL;
static long s_Today;
static long s_FirstFixingDate;
static long s_LastFixingDate;
static double s_DRST1T1;
static double s_XT1Tn, s_DRST1Tn;
static double s_BandWidth;
static double s_Margin;
static long s_N;

/* ---------------------------------------------------------------------------------------------
 */

/* Function to minimise for the LOGNORMAL case (== k for band position centered
 * around k*Libor)  */

static double logdiff_num_f_golden(double K) {
  double S;
  double Alpha, Beta;
  double YT1Tn;
  double d, d1, d2;
  int i;

  /* Upper and lower values of the Band (as a function of K ) */
  Alpha = K - s_BandWidth / (2.0 * s_DRST1T1);
  Beta = K + s_BandWidth / (2.0 * s_DRST1T1);

  YT1Tn = 1.0 / s_XT1Tn - 1;
  S = 0.0;

  /* Loop on all the digitals that form the band */
  for (i = 2; i <= (s_N + 1); i++) {

    /* Linearly interpolate the Forwards between T1 and Tn */
    d = -log(1.0 + s_Theta[i] * YT1Tn) / s_FwdStDev[i] + 0.5 * s_FwdStDev[i];

    d1 = d + log(Beta) / s_FwdStDev[i];
    d2 = d + log(Alpha) / s_FwdStDev[i];
    S += (norm(d2) - norm(d1));
  }

  return S;
}

/* ---------------------------------------------------------------------------------------------
 */

/* Function to minimise for the NORMAL case (== k for band position centered
 * around k*Libor)  */

static double normdiff_num_f_golden(double K) {
  double S;
  double d, d1, d2;
  int i;
  double *SlopeT1;

  SlopeT1 = dvector(2, (s_N + 1));
  S = 0.0;

  for (i = 2; i <= (s_N + 1); i++) {

    /* Linearly interpolate the Forwards between T1 and Tn */
    SlopeT1[i] = s_Theta[i] * s_DRST1Tn - (K - 1 + s_Theta[i]) * s_DRST1T1;
    d = -SlopeT1[i] / s_FwdStDev[i];

    d1 = d + s_BandWidth / 2.0 / s_FwdStDev[i];
    d2 = d - s_BandWidth / 2.0 / s_FwdStDev[i];
    S += (norm(d2) - norm(d1));
  }

  return S;
}

/* ---------------------------------------------------------------------------------------------
 */

/* Integrand for the integration with respect to the joint probability
 * distribution of (DRST1T1  ,DRST1Tn)*/

static double srt_f_accrual_gauss_hermite_integrand(
    double u, double v, double Correl, double Var_1, double Var_2,
    double Mean_1, double Mean_2, double OptStrike, double Bandwidth,
    SrtResOptType ValResetType, SrtDiffusionType ValDiffusion)

{
  double Alpha, Beta;
  double LnXT1Tn, LnDRST1T1;
  double SqCorrel, CorrResid;
  double MeanLnDRST1T1, MeanLnXT1n, SqMeanLnXT1n, SqMeanLnDRST1T1;
  double LnDRST1Var, LnXT1nVar, LnXT1nStd, DRST1Std, DRST1T1Var, DRST1TnVar;
  double d1, d2;
  double SDeltaN;
  double Integrand;
  double d;
  int i;

  SqCorrel = Correl * Correl;
  CorrResid = 1 - SqCorrel;
  s_BandWidth = Bandwidth;

  if (ValDiffusion == SRT_LOGNORMAL) {

    /* LOGNORMAL Diffusion Case */

    LnXT1nVar = Var_1;
    LnDRST1Var = Var_2;
    MeanLnXT1n = Mean_1;
    MeanLnDRST1T1 = Mean_2;

    SqMeanLnXT1n = MeanLnXT1n * MeanLnXT1n;
    SqMeanLnDRST1T1 = MeanLnDRST1T1 * MeanLnDRST1T1;
    LnXT1nStd = sqrt(LnXT1nVar);
    DRST1Std = sqrt(LnDRST1Var);

    LnXT1Tn = sqrt(2 * CorrResid) * LnXT1nVar * u + MeanLnXT1n;
    s_XT1Tn = exp(LnXT1Tn);
    LnDRST1T1 = sqrt(2 * CorrResid) * LnDRST1Var * v + MeanLnDRST1T1;
    s_DRST1T1 = exp(LnDRST1T1);

    if (ValResetType == SRT_OPTIMISED) {

      Beta = OptStrike + s_BandWidth / (2.0 * s_DRST1T1);
      Alpha = OptStrike - s_BandWidth / (2.0 * s_DRST1T1);
    }

    else if (ValResetType == SRT_AUTORESET) {

      Beta = 1.0 + s_BandWidth / (2.0 * s_DRST1T1);
      Alpha = 1.0 - s_BandWidth / (2.0 * s_DRST1T1);
    }

    SDeltaN = 0.0;

    for (i = 2; i <= (s_N + 1); i++) {

      d = -log(1 + s_Theta[i] * ((1.0 / s_XT1Tn) - 1.0)) / s_FwdStDev[i] +
          0.5 * s_FwdStDev[i];
      d1 = d + log(Beta) / s_FwdStDev[i];
      d2 = d + log(Alpha) / s_FwdStDev[i];
      SDeltaN += norm(d1) - norm(d2);
    }

    /*
    Integrand = SDeltaN*sqrt(CorrResid)*exp(2*Correl*u*v)/(SRT_PI);
    */
    Integrand = SDeltaN * sqrt(CorrResid) * (s_DRST1T1 + s_Margin) *
                exp(2 * Correl * u * v) / (SRT_PI);

  }

  else if (ValDiffusion == SRT_NORMAL) {

    /* NORMAL Diffusion Case */

    DRST1T1Var = Var_1;
    DRST1TnVar = Var_2;
    s_DRST1T1 = Mean_1 + u * sqrt(2 * CorrResid * DRST1T1Var);
    s_DRST1Tn = Mean_2 + v * sqrt(2 * CorrResid * DRST1TnVar);

    SDeltaN = 0.0;

    for (i = 2; i <= (s_N + 1); i++) {

      d = (-s_Theta[i] * s_DRST1Tn + (OptStrike - 1 + s_Theta[i]) * s_DRST1T1) /
          s_FwdStDev[i];
      d1 = d + s_BandWidth / 2.0 / s_FwdStDev[i];
      d2 = d - s_BandWidth / 2.0 / s_FwdStDev[i];
      SDeltaN += (norm(d1) - norm(d2));
    }

    Integrand = SDeltaN * sqrt(CorrResid) * (s_DRST1T1 + s_Margin) *
                exp(2 * Correl * u * v) / (SRT_PI);
  }

  return Integrand;

} /* END static double srt_f_accrual_gauss_hermite_integrand(...) */

/* -------------------------------------------------------------------------------
 */

/* The MAIN function to compute the price of an Accrual Resettable Range Floater
 */

Err srt_f_rrf_accrual(String ResetType, long Today, long FirstFixingDate,
                      long LastFixingDate, double Margin, double BandWidth,
                      double DRS0T1, double DRS0Tn, double Correlation,
                      String VolType, double CapBsVolT1, double CapBsVolTn,
                      double FwdVolTn, double DfPaymentDate,
                      int NumHermitePoints, double *Rrf_Accrual_Premium)

{

  Err err;
  double A, B, C, BB, MTerm;
  double DRSnVolT1, SqDRSnVolT1, LnXT1Tn, LnDRST1T1;
  double LnXT1nVar, LnDRST1Var, DRST1Var;
  double LnXT1nStd, LnDRST1Std, DRST1Std;
  double ChangeMeasureTerm;
  double Correl, SqCorrel, CorrResid;
  double ax, bx, cx;
  double SqCapBsVolT1, SqCapBsVolTn;
  double P1, P2, Term1, Term2;
  double MeanLnXT1n, NewMeanLnXT1n;
  double MeanLnDRST1T1, NewMeanLnDRST1T1;
  double MatT0Tn, MatT0T1, MatT1Tn, MatT1Tk;
  int k, i, j;
  double Premium;
  SrtResOptType ValResetType;
  SrtDiffusionType ValDiffusion;
  double *SlopeStdT0T1 = NULL;
  double **OptStrike = NULL;
  double *HermiteWeights = NULL;
  double *HermiteAbscissas = NULL;

  /* Transforms strings into tyep: vol type  , reset type */
  err = interp_reset_optimised(ResetType, &ValResetType);
  if (err)
    return err;
  err = interp_diffusion_type(VolType, &ValDiffusion);
  if (err)
    return err;

  /* Memory allocation */

  HermiteWeights = dvector(1, NumHermitePoints);
  HermiteAbscissas = dvector(1, NumHermitePoints);

  s_N = LastFixingDate - FirstFixingDate;
  SlopeStdT0T1 = dvector(2, (s_N + 1));
  OptStrike = dmatrix(1, NumHermitePoints, 1, NumHermitePoints);

  /* Some Static Initialisations */

  s_FirstFixingDate = FirstFixingDate;
  s_LastFixingDate = LastFixingDate;
  s_BandWidth = BandWidth;
  s_Margin = Margin;
  s_Today = Today;

  /* Memory allocation for the statics vectors */

  s_Theta = dvector(2, (s_N + 1));
  s_SlopeVarT0T1 = dvector(2, (s_N + 1));
  s_FwdVol = dvector(2, (s_N + 1));
  s_FwdStDev = dvector(2, (s_N + 1));
  s_SlopeT0 = dvector(2, (s_N + 1));

  /* Sets the maturities of the different periods */
  MatT0T1 = (s_FirstFixingDate - s_Today) * YEARS_IN_DAY;
  MatT0Tn = (s_LastFixingDate - s_Today) * YEARS_IN_DAY;
  MatT1Tn = (s_LastFixingDate - s_FirstFixingDate) * YEARS_IN_DAY;

  SqCapBsVolT1 = CapBsVolT1 * CapBsVolT1;
  SqCapBsVolTn = CapBsVolTn * CapBsVolTn;

  /* Sets the time ratio (theta)  , the fwdvol  , the cap vol and the square */
  for (k = 2; k <= (s_N + 1); k++) {
    s_Theta[k] = (double)(k - 1) / s_N;
    s_FwdVol[k] = FwdVolTn;
    s_FwdStDev[k] = FwdVolTn * sqrt((k - 1) * YEARS_IN_DAY);
  }

  /* Compute Volatility of the DRStTn until T1 (caplet vol - forward vol ) */

  SqDRSnVolT1 = (SqCapBsVolTn * MatT0Tn -
                 MatT1Tn * s_FwdVol[s_N + 1] * s_FwdVol[s_N + 1]) /
                MatT0T1;
  if (SqDRSnVolT1 < 0) {
    FREE_ALL_ACCRUAL_RRF_MEMORY;
    return serror("Fatal Error in DRSnVolT1 Computation");
  }
  DRSnVolT1 = sqrt(SqDRSnVolT1);

  /* Compute the Hermite points for the numerical integration */
  err = gauss_hermite(HermiteAbscissas, HermiteWeights, NumHermitePoints);
  if (err) {
    FREE_ALL_ACCRUAL_RRF_MEMORY;
    return err;
  }

  /* Computation of price */
  if (ValDiffusion == SRT_LOGNORMAL) {
    /* ----------------- LOGNORMAL CASE ------------------------ */

    DRSnVolT1 = sqrt(SqDRSnVolT1);
    LnXT1nVar = (SqCapBsVolT1 + SqDRSnVolT1 -
                 2 * Correlation * CapBsVolT1 * DRSnVolT1) *
                MatT0T1;
    LnXT1nStd = sqrt(LnXT1nVar);
    LnDRST1Var = SqCapBsVolT1 * MatT0T1;
    LnDRST1Std = sqrt(LnDRST1Var);

    ChangeMeasureTerm =
        (SqCapBsVolT1 - Correlation * CapBsVolT1 * DRSnVolT1) * MatT0T1;
    MeanLnXT1n =
        log(DRS0T1 / DRS0Tn) + 0.5 * (SqCapBsVolT1 - SqDRSnVolT1) * MatT0T1;
    NewMeanLnXT1n = MeanLnXT1n + ChangeMeasureTerm;
    MeanLnDRST1T1 = log(DRS0T1) - 0.5 * LnDRST1Var;
    NewMeanLnDRST1T1 = log(DRS0T1) + 0.5 * LnDRST1Var;

    Correl = (LnDRST1Var - DRSnVolT1 * CapBsVolT1 * MatT0T1 * Correlation) /
             (LnDRST1Std * LnXT1nStd);
    SqCorrel = Correl * Correl;
    CorrResid = 1 - SqCorrel;

    /* LOGNORMAL / RESET & OPTIMISED */

    P1 = 0.0;
    P2 = 0.0;

    /* Double integration with respect to: ( DRST1Tn   , XT1Tn = DRST1T1 /
     * DRST1Tn) */
    for (i = 1; i <= NumHermitePoints; i++) {
      Term1 = 0.0;
      Term2 = 0.0;

      for (j = 1; j <= NumHermitePoints; j++) {

        if (ValResetType == SRT_OPTIMISED) {
          /* OPTIMISED case */

          LnXT1Tn = sqrt(2 * CorrResid) * LnXT1nVar * HermiteAbscissas[i] +
                    MeanLnXT1n;
          s_XT1Tn = exp(LnXT1Tn);
          LnDRST1T1 = sqrt(2 * CorrResid) * LnDRST1Var * HermiteAbscissas[j] +
                      MeanLnDRST1T1;
          s_DRST1T1 = exp(LnDRST1T1);

          s_DRST1Tn = s_DRST1T1 / s_XT1Tn;

          /* Sets proper upper and lower values for the band position */
          if (s_DRST1Tn >= s_DRST1T1) {
            ax = 1 - 0.5 * s_BandWidth / s_DRST1T1;
            bx = (s_DRST1Tn + 0.5 * s_BandWidth) / s_DRST1T1;
          } else if (s_DRST1T1 >= s_DRST1Tn)

          {

            bx = 1 + 0.5 * s_BandWidth / s_DRST1T1;
            ax = (s_DRST1Tn - 0.5 * s_BandWidth) / s_DRST1T1;
          }
          if (ax * s_DRST1T1 - 0.5 * s_BandWidth < 0.0)
            ax = (TINY + 0.5 * s_BandWidth) / s_DRST1T1;

          /* Finds the optimal position of the band */
          cx = (ax + bx) * 0.5;
          golden_section(ax, cx, bx, logdiff_num_f_golden, Tolerance,
                         &OptStrike[i][j]);

        } /* END Optimised case */

        else if (ValResetType == SRT_AUTORESET) {
          /* In the auto reset  , the band is centerd on the initial Libor ( k =
           * 1 ) */
          OptStrike[i][j] = 1;
        }

        /* Increment the integral by the current value (Hermite integration)
         * (first term for Margin) */
        Term1 += HermiteWeights[j] *
                 srt_f_accrual_gauss_hermite_integrand(
                     HermiteAbscissas[i], HermiteAbscissas[j], Correl,
                     LnXT1nVar, LnDRST1Var, MeanLnXT1n, MeanLnDRST1T1,
                     OptStrike[i][j], BandWidth, ValResetType, ValDiffusion);
      }
      /* Multiplies the current value by the appropriate Hermite Weight */
      Term1 *= HermiteWeights[i];

      /* Adds the current value to the integral value */
      P1 += Term1;

    } /* END of first loop for integration */

    /* The Final premium value  */
    Premium = DfPaymentDate * P1;

  } /* END of Lognormal case */

  else if (ValDiffusion == SRT_NORMAL) {
    /* ----------------- NORMAL CASE ------------------------ */

    DRST1Var = MatT0T1 * SqCapBsVolT1;
    DRST1Std = sqrt(DRST1Var);
    s_SlopeVarT0T1[s_N + 1] =
        MatT0T1 *
        (SqDRSnVolT1 - 2 * Correlation * CapBsVolT1 * DRSnVolT1 + SqCapBsVolT1);
    SlopeStdT0T1[s_N + 1] = sqrt(s_SlopeVarT0T1[s_N + 1]);
    Correl = (DRST1Var - DRSnVolT1 * CapBsVolT1 * MatT0T1 * Correlation) /
             (DRST1Std * SlopeStdT0T1[s_N + 1]);
    SqCorrel = Correl * Correl;
    CorrResid = 1 - SqCorrel;

    P1 = 0.0;
    P2 = 0.0;

    if (ValResetType == SRT_AUTORESET) {

      /* AUTO_RESET case */

      for (k = 2; k <= (s_N + 1); k++) {

        s_SlopeT0[k] = s_Theta[k] * (DRS0Tn - DRS0T1);

        if (k != (s_N + 1)) {
          s_SlopeVarT0T1[k] =
              s_Theta[k] * s_Theta[k] * MatT0T1 *
              (SqDRSnVolT1 - 2 * Correlation * CapBsVolT1 * DRSnVolT1 +
               SqCapBsVolT1);
          SlopeStdT0T1[k] = sqrt(s_SlopeVarT0T1[k]);
        }

        MatT1Tk = (double)((k - 1) * YEARS_IN_DAY);
        A = (s_BandWidth / 2.0 - s_SlopeT0[k]) / (Correl * SlopeStdT0T1[k]);
        B = sqrt(MatT1Tk * s_FwdVol[k] * s_FwdVol[k] +
                 (1 - SqCorrel) * s_SlopeVarT0T1[k]) /
            (Correl * SlopeStdT0T1[k]);
        C = -(s_BandWidth / 2.0 + s_SlopeT0[k]) / (Correl * SlopeStdT0T1[k]);

        /* The nice part of the equations ... (see document) */
        MTerm = srt_f_intnrm_m2(A, B, 0, 1) - srt_f_intnrm_m1(A, B, 0, 1) -
                srt_f_intnrm_m2(C, B, 0, 1) + srt_f_intnrm_m1(C, B, 0, 1);

        BB = sqrt(MatT1Tk * s_FwdVol[k] * s_FwdVol[k] + s_SlopeVarT0T1[k]) /
             (Correl * SlopeStdT0T1[k]);

        P1 += (DRS0T1 + Margin) * (norm(A / BB) - norm(C / BB));
        P2 += sqrt(MatT0T1) * CapBsVolT1 * MTerm;
      }

      /* The Premium is the sum of both terms */
      Premium = DfPaymentDate * (P1 + P2);

    } /* END of AutoReset (Normal) case */

    else if (ValResetType == SRT_OPTIMISED) {
      /* Normal OPTIMISED case */

      P1 = 0.0;

      for (i = 1; i <= NumHermitePoints; i++) {
        Term1 = 0;

        for (j = 1; j <= NumHermitePoints; j++) {

          /* Sets the Values of the DRS's on start date (T1) and end date (Tn)
           */
          s_DRST1T1 = DRS0T1 + HermiteAbscissas[i] *
                                   sqrt(2 * CorrResid * MatT0T1 * SqCapBsVolT1);
          s_DRST1Tn = DRS0Tn + HermiteAbscissas[j] *
                                   sqrt(2 * CorrResid * MatT0Tn * SqCapBsVolTn);

          /* Sets proper upper and lower values for the band position */
          if (s_DRST1Tn >= s_DRST1T1) {
            ax = 1 - 0.5 * s_BandWidth / s_DRST1T1;
            bx = (s_DRST1Tn + 0.5 * s_BandWidth) / s_DRST1T1;
          } else if (s_DRST1T1 >= s_DRST1Tn)

          {
            bx = 1 + 0.5 * s_BandWidth / s_DRST1T1;
            ax = (s_DRST1Tn - 0.5 * s_BandWidth) / s_DRST1T1;
          }

          /* Finds the optimal band position (== OptStrike) using Golden section
           */
          cx = 0.5 * (ax + bx);
          golden_section(ax, cx, bx, normdiff_num_f_golden, Tolerance,
                         &OptStrike[i][j]);

          /* Computes the value with the optimal strike */
          Term1 += HermiteWeights[j] *
                   srt_f_accrual_gauss_hermite_integrand(
                       HermiteAbscissas[i], HermiteAbscissas[j], Correl,
                       MatT0T1 * SqCapBsVolT1, MatT0Tn * SqCapBsVolTn, DRS0T1,
                       DRS0Tn, OptStrike[i][j], s_BandWidth, ValResetType,
                       ValDiffusion);
        }

        /* Multiplies the current integrand value by the first Hermite weight */
        Term1 *= HermiteWeights[i];

        /* Increment the integral by the current value */
        P1 += Term1;

      } /* END of first loop of numerical integration (DRS T1) */

      /* The final premium */
      Premium = DfPaymentDate * P1;
    }
  }

  *Rrf_Accrual_Premium = Premium;

  /* Free all used memory */
  FREE_ALL_ACCRUAL_RRF_MEMORY;

  /* Returns a success message */
  return NULL;

} /* END Err srt_f_accrual_rrf(...) */

#undef Tolerance
#undef SHIFT3
#undef GOLD
#undef GLIMIT
#undef TINY
