
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "FxSabrAdi.h"
#include "LGMSVPDE.h"
#include "cpdcalib.h"
#include "lgmsvclosedform.h"
#include "math.h"
#include "opfnctns.h"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define PI2 6.28318530717958647692528676655900576839433879875021164194988918462
#define LGM_VEGA_SHIFT 0.025
#define MAX_CPN 600
#define ONE_MONTH 0.083333333

void LGMSV_prod_comp(double re1, double im1, double re2, double im2,
                     double *re3, double *im3) {
  *re3 = re1 * re2 - im1 * im2;
  *im3 = re1 * im2 + im1 * re2;
}

void LGMSV_div_comp(double re1, double im1, double re2, double im2, double *re3,
                    double *im3) {
  static double div;

  div = re2 * re2 + im2 * im2;
  LGMSV_prod_comp(re1, im1, re2, -im2, re3, im3);
  *re3 /= div;
  *im3 /= div;
}

void LGMSV_sqr_comp(double re1, double im1, double *re3, double *im3) {
  static double norm;

  norm = sqrt(re1 * re1 + im1 * im1);
  *re3 = sqrt((re1 + norm) / 2.0);
  *im3 = sqrt((-re1 + norm) / 2.0);
  if (im1 < 0.0) {
    *im3 *= -1.0;
  }
}

void LGMSVMomentInitOnVol(
    /* Inputs */
    double dLambdaX, double dAlpha, double dLambdaEps, double dRho,

    LGMSVSolFunc *FuncPsi, LGMSVSolFunc *FuncPsi2, LGMSVSolFunc *FuncV2,
    LGMSVSolFunc *FuncPsiV) {
  FuncPsi->dLambda = 0.0;
  FuncPsi->bIsft1 = 1;
  FuncPsi->bIsgt1 = 1;
  FuncPsi->b = 0;
  FuncPsi->bIsht1 = 1;
  FuncPsi->c = 0;

  FuncPsi2->dLambda = 0.0;
  FuncPsi2->bIsft1 = 0;
  FuncPsi2->pft = FuncPsiV;
  FuncPsi2->bIsgt1 = 1;
  FuncPsi2->b = 0.0;
  FuncPsi2->bIsht1 = 1;
  FuncPsi2->c = 0.0;

  FuncV2->dLambda = 2.0 * dLambdaEps;
  FuncV2->bIsft1 = 1;
  FuncV2->a = 2.0 * dLambdaEps + dAlpha * dAlpha;
  FuncV2->bIsgt1 = 1;
  FuncV2->b = 0.0;
  FuncV2->bIsht1 = 1;
  FuncV2->c = 0.0;

  FuncPsiV->dLambda = 2.0 * dLambdaX + dLambdaEps;
  FuncPsiV->bIsft1 = 0;
  FuncPsiV->pft = FuncPsi;
  FuncPsiV->a = dLambdaEps;
  FuncPsiV->bIsgt1 = 0;
  FuncPsiV->pgt = FuncV2;
  FuncPsiV->bIsht1 = 1;
  FuncPsiV->c = 0.0;
}

void LGMSVMomentCalculationOnVol(

    double Sig2, double dt,

    double InitPsi, double InitPsi2, double InitV2, double InitPsiV,

    LGMSVSolFunc *FuncPsi, LGMSVSolFunc *FuncPsi2, LGMSVSolFunc *FuncV2,
    LGMSVSolFunc *FuncPsiV, double *lambda,

    double *ResPsi, double *ResPsi2, double *ResV2, double *ResPsiV) {
  /* Initialisation */
  FuncPsi->a = Sig2;
  FuncPsi->dXt1 = InitPsi;
  FuncPsi2->a = 2.0 * Sig2;
  FuncPsi2->dXt1 = InitPsi2;
  FuncV2->dXt1 = InitV2;
  FuncPsiV->b = Sig2;
  FuncPsiV->dXt1 = InitPsiV;

  /* Calculation */
  LGMSVFuncValue2(FuncPsi, dt, lambda, 0, ResPsi);
  LGMSVFuncValue2(FuncPsi2, dt, lambda, 0, ResPsi2);
  LGMSVFuncValue2(FuncV2, dt, lambda, 0, ResV2);
  LGMSVFuncValue2(FuncPsiV, dt, lambda, 0, ResPsiV);
}

/* -------------------------------------------------------------------------------------------------------------
        LGMSVFillPayoff

  --------------------------------------------------------------------------------------------------------------
*/
void LGMSVFillPayoffOnVol(/* Input */
                          int iNbPsi, int iNbDrift, double dLambdaX,
                          double dRho, double dPsitMean, int iIndexPsiMean,
                          double dPsiStep, int iIndexDriftMean,
                          double dDriftStep,
                          double dExTime, /* Exercice of the swaption in years
                                             from today */
                          double dTStar,  /* Tstar */
                          int iNbCoupon,  /* Description of the cashflows */
                          double *CouponTime, double *Coupon,

                          /* Outputs */
                          double **Payoff) {
  int iNumDff, iNumPsi, iNumDrift;
  int LimPsi, LimDrift;
  double CoefDriftLim, CoefPsiLim;
  double *beta;
  double Psi, Drift;
  double ystar, ylim;
  double rho2, sqrho2, sqpsi, stdf;
  int dir;

  /* Initialisation */
  beta = dvector(0, iNbCoupon - 1);
  for (iNumDff = 0; iNumDff < iNbCoupon; iNumDff++) {
    beta[iNumDff] =
        (1.0 - exp(-dLambdaX * (CouponTime[iNumDff] - dTStar))) / dLambdaX;
  }

  dDriftStep *= dRho;
  LimPsi = iNbPsi - iIndexPsiMean + 1;
  LimDrift = iNbDrift - iIndexDriftMean + 1;
  CoefDriftLim = -iNbDrift * dDriftStep;
  CoefPsiLim = -iNbPsi * dPsiStep;
  rho2 = dRho * dRho;
  sqrho2 = sqrt(1.0 - rho2);

  Psi = dPsitMean;

  for (iNumPsi = 1; iNumPsi <= iNbPsi; iNumPsi++) {
    ystar = static_lgmystar_stochvol(iNbCoupon, Coupon, beta, Psi, &dir);
    sqpsi = sqrt(Psi);
    stdf = sqpsi * sqrho2;

    Drift = 0.0;

    for (iNumDrift = 1; iNumDrift <= iNbDrift; iNumDrift++) {
      ylim = (ystar - Drift) / stdf;

      /*
      Payoff[iNumPsi][iNumDrift] = Drift;
      */

      for (iNumDff = 0; iNumDff < iNbCoupon; iNumDff++) {
        Payoff[iNumPsi][iNumDrift] +=
            Coupon[iNumDff] *
            exp(-beta[iNumDff] * (Drift + 0.5 * beta[iNumDff] * Psi * rho2)) *
            static_lgmsafenorm(dir * (ylim + stdf * beta[iNumDff]));
      }

      Drift += dDriftStep;
      if (iNumDrift == LimDrift) {
        Drift += CoefDriftLim;
      }
    }

    Psi += dPsiStep;
    if (iNumPsi == LimPsi) {
      Psi += CoefPsiLim;
    }
  }
}

void LGMSVOptionPricePrecalcOnVol(
    int iNbPhi, int iNbft,

    double dLambdaX, double dRho, double dTStar, int endi,

    double *dt, double Coef1Re, double Coef2Re, double *Coef2ImT,
    double *Coef3ReT, double *Coef3ImT,

    double dPhitMean, int SaveFile, double ***Density, double **DensitySpeq,
    double **PayOff,

    double dExTime, int iNbCoupon, double *CouponTime, double *Coupon,

    double dPhiFreqStep, double dPhiStep, int iIndexPhiMean, double dftFreqStep,
    double dftStep, int iIndexft0,

    double *Price) {
  int iNumPhi, iNumft;
  double dIntegral;
  clock_t time1, time2, time3;

  time1 = clock();

  LGMSVCalculateDensityPrecalcNew(
      iNbPhi, iNbft, dPhiFreqStep, dftFreqStep, endi, dt, Coef1Re, Coef2Re,
      Coef2ImT, Coef3ReT, Coef3ImT, dPhitMean, SaveFile, Density, DensitySpeq);

  time2 = clock();

  LGMSVFillPayoffOnVol(iNbPhi, iNbft, dLambdaX, dRho, dPhitMean, iIndexPhiMean,
                       dPhiStep, iIndexft0, dftStep, dExTime, dTStar, iNbCoupon,
                       CouponTime, Coupon, PayOff);

  /* Calculation of the Integral of the payoff times the density */
  dIntegral = 0;
  for (iNumPhi = 1; iNumPhi <= iNbPhi; iNumPhi++)
    for (iNumft = 1; iNumft <= iNbft; iNumft++)
      /* Calculation of the payoff */
      dIntegral += Density[1][iNumPhi][iNumft] * PayOff[iNumPhi][iNumft];
  dIntegral *= dPhiStep * dftStep;

  time3 = clock();

  /* Return the result */
  *Price = dIntegral;
}

/* -------------------------------------------------------------------------------------------------------------
        LGMSVClosedForm

  --------------------------------------------------------------------------------------------------------------
*/
void LGMSVClosedFormOnVol(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    double dLambdaX, int iNbSigTime,             /* Term Structure of g(t) */
    double *SigTime, double *Sig, double dTStar, /* Tstar in years from today */
    double dAlpha, double dLambdaEps, double dRho,

    /* Product description */
    long lExDate,   /* Exercice date of the swaption  */
    double dExTime, /* Exercice of the swaption in years from today */
    int iNbCoupon,  /* Description of the cashflows */
    double *CouponTime, long *CouponDate, double *Coupon,
    char *cYieldCurve, /* Yield Curve */

    /* Parameter of grids */
    int iNbPsi,   /* Number of Psi : Should be a power of two */
    int iNbDrift, /* Number of ft : Should be a power of two */
    double iNbSigmaPsiGridLeft, double iNbSigmaPsiGridRight,
    double iNbSigmaftLeft, double iNbSigmaftRight,

    double iRatioPsi, double iRatioFt, int iPriorityFreqPsi,
    int iPriorityFreqFt, double iMaxTime, int SaveFile,

    /* Outputs */
    double *Price) {
  /* Declaration of locals variables */
  Err err = NULL;
  double lambdaArray[10];

  /* For moments calculation */
  LGMSVSolFunc FuncPsi_, *FuncPsi = &FuncPsi_;
  LGMSVSolFunc FuncPsi2_, *FuncPsi2 = &FuncPsi2_;
  LGMSVSolFunc FuncV2_, *FuncV2 = &FuncV2_;
  LGMSVSolFunc FuncPsiV_, *FuncPsiV = &FuncPsiV_;
  double ExpectPsi, ExpectPsi2, ExpevtV2, ExpectPsiV;

  double ***Density = NULL, **Densityspeq = NULL, **PayOff = NULL;

  double dPsitMean, dPsitStd, dPsiStep;
  double dDriftStep;
  double dPsiFreqStep, dftFreqStep;
  int iIndexPsiMean, iIndexDriftMean;

  double LimitPsi, Limitft;
  int i, endi;
  double Coef1Re, Coef2Re, Coef2ImInit, Coef3ReInit;
  double t1, t2, sig1, AlphaEq;

  double *dt = NULL, *Coef2ImT = NULL, *Coef3ReT = NULL, *Coef3ImT = NULL;

  double smooth_factor = 0.99999;

  clock_t time1, time2;

  time1 = clock();

  /* Initialisation and memory allocation */
  PayOff = dmatrix(1, iNbPsi, 1, iNbDrift);

  Density = f3tensor(1, 1, 1, iNbPsi, 1, iNbDrift);

  Densityspeq = dmatrix(1, 1, 1, 2 * iNbPsi);

  endi = Get_Index(dExTime, SigTime, iNbSigTime);

  dt = dvector(0, endi);
  Coef2ImT = dvector(0, endi);
  Coef3ReT = dvector(0, endi);
  Coef3ImT = dvector(0, endi);

  /* Gestion of allocation errors */
  if (!Density || !Densityspeq || !dt || !Coef2ImT || !Coef3ReT || !Coef3ImT) {
    err = "Memory allocation error (1) in LGMSVClosedForm2New";
    goto FREE_RETURN;
  }

  memset(lambdaArray, 0, 10 * sizeof(double));

  /* Precalculations */

  /* Constant independent on time */
  Coef1Re = -0.5 * dAlpha * dAlpha;
  Coef2Re = dLambdaEps;
  Coef2ImInit = -dAlpha * smooth_factor;
  Coef3ReInit = 0.5;

  Limitft = log(iRatioFt);
  LimitPsi = log(iRatioPsi);

  LGMSVMomentInitOnVol(0.0, dAlpha, dLambdaEps, dRho, FuncPsi, FuncPsi2, FuncV2,
                       FuncPsiV);
  ExpectPsi = 0.0;
  ExpectPsi2 = 0.0;
  ExpevtV2 = 1.0;
  ExpectPsiV = 0.0;

  t1 = 0.0;
  for (i = 0; i <= endi; i++) {
    /* Precalculation on the option i*/

    /* First the time discretisation */
    if (i < endi) {
      t2 = SigTime[i];
    } else {
      t2 = dExTime;
    }

    /* Calculate constant values */
    sig1 = Sig[i];
    dt[i] = (t2 - t1);
    Coef2ImT[i] = Coef2ImInit * sig1;
    Coef3ReT[i] = Coef3ReInit * sig1 * sig1;
    Coef3ImT[i] = -sig1 * sig1;

    LGMSVMomentCalculationOnVol(sig1 * sig1, t2 - t1, ExpectPsi, ExpectPsi2,
                                ExpevtV2, ExpectPsiV, FuncPsi, FuncPsi2, FuncV2,
                                FuncPsiV, lambdaArray, &ExpectPsi, &ExpectPsi2,
                                &ExpevtV2, &ExpectPsiV);

    t1 = t2;
  }

  AlphaEq = dAlpha / 2.0;
  if (dLambdaEps > 1.0E-16) {
    AlphaEq *= sqrt((1.0 - exp(-dLambdaEps * t2)) / (dLambdaEps * t2));
  }

  dPsitMean = ExpectPsi;
  dPsitStd = sqrt(ExpectPsi2 - ExpectPsi * ExpectPsi);

  /* Find Frequence */
  LGMSVFindFreqNew(iNbPsi, iNbDrift, AlphaEq, 0.0, endi, dt, Coef1Re, Coef2Re,
                   Coef2ImT, Coef3ReT, Coef3ImT, dPsitMean, dPsitStd,
                   sqrt(dPsitMean), iNbSigmaPsiGridLeft, iNbSigmaPsiGridRight,
                   iNbSigmaftLeft, iNbSigmaftRight, LimitPsi, Limitft,
                   iPriorityFreqPsi, iPriorityFreqFt, &dPsiFreqStep, &dPsiStep,
                   &iIndexPsiMean, &dftFreqStep, &dDriftStep, &iIndexDriftMean);

  /* First pricing */
  LGMSVOptionPricePrecalcOnVol(
      iNbPsi, iNbDrift, dLambdaX, dRho, dTStar, endi, dt, Coef1Re, Coef2Re,
      Coef2ImT, Coef3ReT, Coef3ImT, dPsitMean, SaveFile, Density, Densityspeq,
      PayOff, dExTime, iNbCoupon, CouponTime, Coupon, dPsiFreqStep, dPsiStep,
      iIndexPsiMean, dftFreqStep, dDriftStep, iIndexDriftMean, Price);

  /* Save the Payoff and Density */
  if (SaveFile) {
    LGMSVSaveDensity(iNbPsi, iNbDrift, iIndexPsiMean, iIndexDriftMean, dPsiStep,
                     dDriftStep, dPsitMean, dPsitStd, Density);

    LGMSVSavePayOff(iNbPsi, iNbDrift, iIndexPsiMean, iIndexDriftMean, dPsiStep,
                    dDriftStep, dPsitMean, PayOff);
  }

  time2 = clock();

FREE_RETURN:

  /* free memory */
  if (Density)
    free_f3tensor(Density, 1, 1, 1, iNbPsi, 1, iNbDrift);

  if (Densityspeq)
    free_dmatrix(Densityspeq, 1, 1, 1, iNbPsi);

  if (PayOff)
    free_dmatrix(PayOff, 1, iNbPsi, 1, iNbDrift);

  if (dt)
    free_dvector(dt, 0, endi);
  if (Coef2ImT)
    free_dvector(Coef2ImT, 0, endi);
  if (Coef3ReT)
    free_dvector(Coef3ReT, 0, endi);
  if (Coef3ImT)
    free_dvector(Coef3ImT, 0, endi);
}