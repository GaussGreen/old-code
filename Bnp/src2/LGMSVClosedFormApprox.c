
#include "lgmsvclosedformapprox.h"
#include "CTSProdStruct.h"
#include "DiagCalibDLM.h"
#include "Fx3FUtils.h"
#include "FxSabrAdi.h"
#include "LGMSVCalibApprox.h"
#include "LGMSVPDE.h"
#include "lgmsvclosedform.h"
#include "math.h"
#include "opfnctns.h"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define PI2 6.28318530717958647692528676655900576839433879875021164194988918462
#define LGM_VEGA_SHIFT 0.025
#define MAX_CPN 600
#define ONE_MONTH 0.083333333
#define MAX_LGMSHIFT 100

#define LGMSV_MINVOL 0.0001

/* Set the default values for all the CTS autocal parameters */
void lgmsv_app_set_default_params(int *iNbX, double *iNbSigmaXLeft,
                                  double *iNbSigmaXRight, double *dIntegParam,
                                  int *iIntegMethod, double *dVolLimit,
                                  int *iCalibLGM, double *dMinStd,
                                  double *dMaxStd) {
  *iNbX = 350;
  *iNbSigmaXLeft = 1000;
  *iNbSigmaXRight = 100000;
  *dIntegParam = 0.5;
  *iIntegMethod = 3;
  *dVolLimit = 0.0050;
  *iCalibLGM = 1;
  *dMinStd = 7.0;
  *dMaxStd = 7.0;
}

void lgmsv_init_numer_params(LGMSV_NUMERPARAMS NumerParams, int iNbX,
                             double dParam1, double dParam2, double dIntegParam,
                             double dVolLimit, int iCalibLGM, int iIntegMethod,
                             double dMinStd, double dMaxStd) {
  NumerParams->iNbX = iNbX;
  NumerParams->dParam1 = dParam1;
  NumerParams->dParam2 = dParam2;
  NumerParams->dIntegParam = dIntegParam;
  NumerParams->dVolLimit = dVolLimit;
  NumerParams->iCalibLGM = iCalibLGM;
  NumerParams->iIntegMethod = iIntegMethod;
  NumerParams->dMinStd = dMinStd;
  NumerParams->dMaxStd = dMaxStd;
}

void lgmsv_app_set_default_params_struct(LGMSV_NUMERPARAMS NumerParams) {
  lgmsv_app_set_default_params(
      &(NumerParams->iNbX), &(NumerParams->dParam1), &(NumerParams->dParam2),
      &(NumerParams->dIntegParam), &(NumerParams->iIntegMethod),
      &(NumerParams->dVolLimit), &(NumerParams->iCalibLGM),
      &(NumerParams->dMinStd), &(NumerParams->dMaxStd));
}

void LGMSVMomentInitApprox(

    double InitX, double InitV,

    LGMSVSolFunc *FuncX, LGMSVSolFunc *FuncX2, LGMSVSolFunc *FuncV,
    LGMSVSolFunc *FuncV2, LGMSVSolFunc *FuncXV) {
  /*	The equations: */
  /*	dX = -0.5 * (k1 * f(t))^2 * V * dt + k1 * f(t) * Sqr(V) * dW1 */
  /*	dV = (lam2 - (lam2 + rho * alpha * k2 * f(t)) * V) * dt + alpha * Sqr(V)
   * * dW2 */

  /*	dX = -0.5 * (k1 * f(t))^2 * V * dt +... */

  FuncX->dXt1 = InitX;
  FuncX->dLambda = 0.0;
  FuncX->bIsft1 = 0;
  FuncX->pft = FuncV;
  FuncX->bIsgt1 = 1;
  FuncX->b = 0;
  FuncX->bIsht1 = 1;
  FuncX->c = 0;

  /* dX2 = -(k1 * f(t))^2 * XV * dt + (k1 * f(t))^2 * V * dt + ... */
  FuncX2->dXt1 = InitX * InitX;
  FuncX2->dLambda = 0.0;
  FuncX2->bIsft1 = 0;
  FuncX2->pft = FuncXV;
  FuncX2->bIsgt1 = 0;
  FuncX2->pgt = FuncV;
  FuncX2->bIsht1 = 1;
  FuncX2->c = 0;

  /*	dV = (lam2 - (lam2 + rho * alpha * k2 * f(t)) * V) * dt + ... */
  FuncV->dXt1 = InitV;
  FuncV->bIsft1 = 1;
  FuncV->bIsgt1 = 1;
  FuncV->b = 0.0;
  FuncV->bIsht1 = 1;
  FuncV->c = 0.0;

  /*	dV2 = ((2 * lam2 + alpha^2) * V - 2 * (lam2 + rho * alpha * k2 * f(t)) *
   * V^2) * dt + ... */
  FuncV2->dXt1 = InitV * InitV;
  FuncV2->bIsft1 = 0;
  FuncV2->pft = FuncV;
  FuncV2->bIsgt1 = 1;
  FuncV2->b = 0.0;
  FuncV2->bIsht1 = 1;
  FuncV2->c = 0.0;

  /*	dXV = -0.5 * (k1 * f(t))^2 * V^2 * dt + lam2 * X * dt + rho * alpha * k1
  * f(t) * V * dt
  - (lam2 + rho * alpha * k2 * f(t) * XV * dt + ... */
  FuncXV->dXt1 = InitX * InitV;
  FuncXV->bIsft1 = 0;
  FuncXV->pft = FuncV2;
  FuncXV->bIsgt1 = 0;
  FuncXV->pgt = FuncX;
  FuncXV->bIsht1 = 0;
  FuncXV->pht = FuncV;
}

void LGMSVMomentCalculationApprox(

    /* Inputs */
    double dLambdaX, double dAlpha, double dLambdaEps, double dRho,
    double dSigma, double k1, double k2, double k3, double dt,

    LGMSVSolFunc *FuncX, LGMSVSolFunc *FuncX2, LGMSVSolFunc *FuncV,
    LGMSVSolFunc *FuncV2, LGMSVSolFunc *FuncXV,

    double *lambda) {
  double temp_cms, temp_var;
  double temp_lam;

  temp_var = k1 * dSigma;
  temp_var *= temp_var;
  temp_cms = dSigma * dSigma * k3;
  temp_lam = dLambdaEps + dRho * dAlpha * k2 * dSigma;

  /* Initialisation */
  FuncX->a = temp_cms - 0.5 * temp_var;

  FuncX2->a = 2.0 * temp_cms - temp_var;
  FuncX2->b = temp_var;

  FuncV->dLambda = temp_lam;
  FuncV->a = dLambdaEps;

  FuncV2->dLambda = 2.0 * temp_lam;
  FuncV2->a = 2.0 * dLambdaEps + dAlpha * dAlpha;

  FuncXV->dLambda = temp_lam;
  FuncXV->a = temp_cms - 0.5 * temp_var;
  FuncXV->b = dLambdaEps;
  FuncXV->c = dRho * dAlpha * k1 * dSigma;

  /* Calculation */
  LGMSVFuncValue2(FuncX2, dt, lambda, 0, &(FuncX2->dXt1));
  LGMSVFuncValue2(FuncXV, dt, lambda, 0, &(FuncXV->dXt1));
  LGMSVFuncValue2(FuncX, dt, lambda, 0, &(FuncX->dXt1));
  LGMSVFuncValue2(FuncV2, dt, lambda, 0, &(FuncV2->dXt1));
  LGMSVFuncValue2(FuncV, dt, lambda, 0, &(FuncV->dXt1));
}

void LGMSVSolveODE(double FreqRe, double FreqIm, double Fwd, int endi,
                   double *dt, double *CoefInt, double *Coef1ReT,
                   double *Coef1ImT, double *Coef2ReT, double *Coef2ImT,
                   double *Coef3ReT, double *Coef3ImT,

                   /* Outputs */
                   double *LogTFRe, double *LogTFIm) {
  double ARe, AIm, DRe, DIm;
  double Coef2Re, Coef2Im, Coef3Re, Coef3Im;
  double coef1, coef2;
  int i;

  /* Initialisation */
  ARe = 0.0;
  AIm = 0.0;
  DRe = 0.0;
  DIm = 0.0;

  /*
  coef1 = FreqRe * FreqRe - FreqIm * (FreqIm + 1);
  coef2 = FreqRe * (1 + 2.0 * FreqIm);
  */

  coef1 = FreqRe * FreqRe - FreqIm * FreqIm;
  coef2 = 2.0 * FreqRe * FreqIm;

  for (i = endi; i >= 0; i--) {
    /* Solve the ODE between t1 and t2 */
    Coef2Re = Coef2ReT[i] - Coef2ImT[i] * FreqIm;
    Coef2Im = Coef2ImT[i] * FreqRe;
    Coef3Re = Coef3ReT[i] * coef1 - Coef3ImT[i] * FreqIm;
    Coef3Im = Coef3ReT[i] * coef2 + Coef3ImT[i] * FreqRe;

    LGMSVUpdateADClosedForm(dt[i], CoefInt[i], Coef1ReT[i], Coef2Re, Coef2Im,
                            Coef3Re, Coef3Im, &DRe, &DIm, &ARe, &AIm);
  }

  *LogTFRe = ARe + DRe - Fwd * FreqIm;
  *LogTFIm = AIm + DIm + Fwd * FreqRe;
}

void LGMSVSolveODEGride(int iNbX, double *X, double Fwd, int endi, double *dt,
                        double *CoefIntT, double *Coef1ReT, double *Coef1ImT,
                        double *Coef2ReT, double *Coef2ImT, double *Coef3ReT,
                        double *Coef3ImT,

                        double dIntegParam,

                        /* Outputs */
                        double *IntRe, double *IntIm) {
  double FreqRe, FreqIm, ResRe, ResIm;
  double coef1, coef2, temp1, temp2, temp3;
  double dExpReal, dImag;
  int i;

  FreqIm = -(dIntegParam + 1);
  coef1 = dIntegParam * (dIntegParam + 1);
  coef2 = 2.0 * dIntegParam + 1;

  for (i = 0; i < iNbX; i++) {
    FreqRe = X[i];

    LGMSVSolveODE(FreqRe, FreqIm, Fwd, endi, dt, CoefIntT, Coef1ReT, Coef1ImT,
                  Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, &ResRe, &ResIm);

    dExpReal = exp(ResRe);
    dImag = fmod(ResIm, PI2);
    ResRe = dExpReal * cos(dImag);
    ResIm = dExpReal * sin(dImag);

    temp1 = coef1 - FreqRe * FreqRe;
    temp2 = FreqRe * coef2;
    temp3 = temp1 * temp1 + temp2 * temp2;

    IntRe[i] = (ResRe * temp1 + ResIm * temp2) / temp3;
    IntIm[i] = (ResIm * temp1 - ResRe * temp2) / temp3;
  }
}

void LGMSVSolveODEGrideTV(int iNbX, double *X, double Fwd, int endi, double *dt,
                          double *CoefIntT, double *Coef1ReT, double *Coef1ImT,
                          double *Coef2ReT, double *Coef2ImT, double *Coef3ReT,
                          double *Coef3ImT,

                          double dIntegParam,

                          /* Outputs */
                          double *IntRe, double *IntIm) {
  double FreqRe, FreqIm1, FreqIm2, ResRe, ResIm;
  double coef11, coef21, coef12, coef22, temp1, temp2, temp3;
  double dExpReal, dImag;
  int i;

  FreqIm1 = -(1 + dIntegParam);
  FreqIm2 = -(1 - dIntegParam);
  coef11 = dIntegParam * (1 + dIntegParam);
  coef12 = -dIntegParam * (1 - dIntegParam);

  coef21 = 1.0 + 2.0 * dIntegParam;
  coef22 = 1.0 - 2.0 * dIntegParam;

  for (i = 0; i < iNbX; i++) {
    FreqRe = X[i];

    LGMSVSolveODE(FreqRe, FreqIm1, Fwd, endi, dt, CoefIntT, Coef1ReT, Coef1ImT,
                  Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, &ResRe, &ResIm);

    dExpReal = exp(ResRe);
    dImag = fmod(ResIm, PI2);
    ResRe = dExpReal * cos(dImag) - 1.0;
    ResIm = dExpReal * sin(dImag);

    temp1 = coef11 - FreqRe * FreqRe;
    temp2 = FreqRe * coef21;
    temp3 = temp1 * temp1 + temp2 * temp2;

    IntRe[i] = (ResRe * temp1 + ResIm * temp2) / temp3;
    IntIm[i] = (ResIm * temp1 - ResRe * temp2) / temp3;

    LGMSVSolveODE(FreqRe, FreqIm2, Fwd, endi, dt, CoefIntT, Coef1ReT, Coef1ImT,
                  Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, &ResRe, &ResIm);

    dExpReal = exp(ResRe);
    dImag = fmod(ResIm, PI2);
    ResRe = dExpReal * cos(dImag) - 1.0;
    ResIm = dExpReal * sin(dImag);

    temp1 = coef12 - FreqRe * FreqRe;
    temp2 = FreqRe * coef22;
    temp3 = temp1 * temp1 + temp2 * temp2;

    IntRe[i] -= (ResRe * temp1 + ResIm * temp2) / temp3;
    IntIm[i] -= (ResIm * temp1 - ResRe * temp2) / temp3;

    IntRe[i] *= 0.5;
    IntIm[i] *= 0.5;
  }
}

void LGMSVCalculateGride(int iNbX, double XMax, double *X) {
  double dX;
  int i;

  dX = XMax / (iNbX - 1);

  X[0] = 0;
  for (i = 1; i < iNbX; i++) {
    X[i] = X[i - 1] + dX;
  }
}

void LGMSVComputeIntegral(int iNbX, double *X, double *IntRe, double *IntIm,
                          int iNbStrike, double *dLogStrike, double dIntegParam,
                          double *result) {
  int i, j;
  double theta, temp, temp1, temp2;

  for (j = 0; j < iNbStrike; j++) {
    theta = X[0] * dLogStrike[j];
    temp1 = (cos(theta) * IntRe[0] + sin(theta) * IntIm[0]);
    temp = 0.0;

    for (i = 1; i < iNbX; i++) {
      theta = X[i] * dLogStrike[j];
      temp2 = (cos(theta) * IntRe[i] + sin(theta) * IntIm[i]);

      temp += 0.5 * (temp1 + temp2) * (X[i] - X[i - 1]);

      temp1 = temp2;
    }

    result[j] = exp(-dIntegParam * dLogStrike[j]) * temp / PI;
  }
}

void LGMSVComputeIntegralLaguerre(int iNbX, double *X, double *W, double Alpha,
                                  double *IntRe, double *IntIm, int iNbStrike,
                                  double *dLogStrike, double dIntegParam,
                                  double *result) {
  int i, j;
  double theta, temp, temp1;

  if (Alpha == 0) {
    for (j = 0; j < iNbStrike; j++) {
      temp = 0.0;

      for (i = 0; i < iNbX; i++) {
        theta = X[i] * dLogStrike[j];
        temp1 = (cos(theta) * IntRe[i] + sin(theta) * IntIm[i]);
        temp += W[i] * temp1 * exp(X[i]);
      }

      result[j] = exp(-dIntegParam * dLogStrike[j]) * temp / PI;
    }
  } else {
    for (j = 0; j < iNbStrike; j++) {
      temp = 0.0;

      for (i = 0; i < iNbX; i++) {
        theta = X[i] * dLogStrike[j];
        temp1 = (cos(theta) * IntRe[i] + sin(theta) * IntIm[i]);
        if (X[i] > 0) {
          temp += W[i] * temp1 * exp(-Alpha * log(X[i]) + X[i]);
        } else {
          temp += 0.0;
        }
      }

      result[j] = exp(-dIntegParam * dLogStrike[j]) * temp / PI;
    }
  }
}

void LGMSVComputeIntegralLegendre(int iNbX, double *X, double *W, double *IntRe,
                                  double *IntIm, int iNbStrike,
                                  double *dLogStrike, double dIntegParam,
                                  double *result) {
  int i, j;
  double theta, temp, temp1;

  for (j = 0; j < iNbStrike; j++) {
    temp = 0.0;

    for (i = 0; i < iNbX; i++) {
      theta = X[i] * dLogStrike[j];
      temp1 = (cos(theta) * IntRe[i] + sin(theta) * IntIm[i]);
      temp += W[i] * temp1;
    }

    result[j] = exp(-dIntegParam * dLogStrike[j]) * temp / PI;
  }
}

void LGMSVComputeIntegralHermite(int iNbX, double *X, double *W, double *IntRe,
                                 double *IntIm, int iNbStrike,
                                 double *dLogStrike, double dIntegParam,
                                 double *result) {
  int i, j;
  double theta, temp, temp1;

  for (j = 0; j < iNbStrike; j++) {
    temp = 0.0;

    for (i = 0; i < iNbX; i++) {
      theta = X[i] * dLogStrike[j];
      temp1 = (cos(theta) * IntRe[i] + sin(theta) * IntIm[i]);
      temp += W[i] * temp1 * exp(X[i] * X[i]);
    }

    result[j] =
        temp / PI /
        (exp(dIntegParam * dLogStrike[j]) - exp(-dIntegParam * dLogStrike[j]));
  }
}

void LGMSVConstructLegendreGrid(int iNbInt, double *dBreakPoints,
                                int *iNbPoints, double *dX, double *dW) {
  int i;
  int iStartIndex;
  double dStart, dEnd;

  dStart = 0;
  iStartIndex = 0;

  for (i = 0; i < iNbInt; i++) {
    dEnd = dBreakPoints[i];
    gauleg(dStart, dEnd, &(dX[iStartIndex - 1]), &(dW[iStartIndex - 1]),
           iNbPoints[i]);
    dStart = dEnd;
    iStartIndex += iNbPoints[i];
  }
}

void LGMSVClosedFormApprox(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    double dLambdaX,

    int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigma, double *dAlpha, double *dLambdaEps,
    double *dRho,

    double dTStar, /* Tstar in years from today */

    /* Product description */
    double dFwd, int iNbStrike, double *dStrike, double dExTime,
    double dCoefVol,

    double dSwitchTime, double dCoefMeanRev1, double dCoefCMS1,
    double dCoefMeanRev2, double dCoefCMS2,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Outputs */
    double *Price) {
  Err err = NULL;
  double lambdaArray[10];

  double InitX, t1, t2, temp_vol;
  double dXMean, dXStd, dVMean, dVStd, XMax;
  int i, endi, iSwitchIndex, do_switch;

  LGMSVSolFunc FuncX_, *FuncX = &FuncX_, FuncX2_, *FuncX2 = &FuncX2_, FuncV_,
                       *FuncV = &FuncV_, FuncV2_, *FuncV2 = &FuncV2_, FuncXV_,
                       *FuncXV = &FuncXV_;

  double *dt = NULL, *CoefIntT = NULL, *Coef1ReT = NULL, *Coef2ReT = NULL,
         *Coef2ImT = NULL, *Coef3ReT = NULL, *Coef3ImT = NULL;

  double *X = NULL, *W = NULL, *IntRe = NULL, *IntIm = NULL;

  double *LogStrike = NULL;

  double dBreakPoints[4];
  int iBreakNbX[4];

  endi = Get_Index(dExTime, dPWTime, iNbPWTime);
  iSwitchIndex = Get_Index(dSwitchTime, dPWTime, iNbPWTime);

  if (fabs(dSwitchTime - dPWTime[iSwitchIndex]) > 1.0E-08 &&
      dSwitchTime < dExTime) {
    do_switch = 1;
  } else {
    do_switch = 0;

    if (dSwitchTime > dExTime) {
      dSwitchTime = dExTime;
      iSwitchIndex = endi;
    }
  }

  endi += do_switch;

  /* Memory allocation */

  dt = dvector(0, endi);
  CoefIntT = dvector(0, endi);
  Coef1ReT = dvector(0, endi);
  Coef2ReT = dvector(0, endi);
  Coef2ImT = dvector(0, endi);
  Coef3ReT = dvector(0, endi);
  Coef3ImT = dvector(0, endi);

  X = dvector(0, NumerParams->iNbX - 1);
  IntRe = dvector(0, NumerParams->iNbX - 1);
  IntIm = dvector(0, NumerParams->iNbX - 1);

  LogStrike = dvector(0, iNbStrike - 1);

  if (!dt || !CoefIntT || !Coef1ReT || !Coef2ReT || !Coef2ImT || !Coef3ReT ||
      !Coef3ImT || !X || !IntRe || !IntIm || !LogStrike) {
    err = "Memory Allocation failure (1) in LGMSVClosedFormNew";
    goto FREE_RETURN;
  }

  memset(lambdaArray, 0, 10 * sizeof(double));

  /* Precalculations */

  if (dCoefVol < 0.0) {
    dCoefVol *= -1.0;
    NumerParams->dIntegParam = -1.0 - NumerParams->dIntegParam;
  }

  InitX = log(fabs(dFwd));
  LGMSVMomentInitApprox(InitX, 1.0, FuncX, FuncX2, FuncV, FuncV2, FuncXV);

  t1 = 0.0;
  for (i = 0; i <= iSwitchIndex; i++) {
    /* Precalculation on the option i*/

    /* First the time discretisation */
    if (i < iSwitchIndex) {
      t2 = dPWTime[i];
    } else {
      t2 = dSwitchTime;
    }

    /* Calculate constant values */
    temp_vol = dCoefVol * dSigma[i];
    dt[i] = (t2 - t1);
    CoefIntT[i] = dLambdaEps[i];
    Coef1ReT[i] = -0.5 * dAlpha[i] * dAlpha[i];
    Coef2ReT[i] =
        dLambdaEps[i] + dCoefMeanRev1 * dAlpha[i] * dRho[i] * dSigma[i];
    Coef2ImT[i] = -dAlpha[i] * dRho[i] * temp_vol;
    temp_vol *= temp_vol;
    Coef3ReT[i] = 0.5 * temp_vol;
    Coef3ImT[i] = 0.5 * temp_vol - dCoefCMS1 * dSigma[i] * dSigma[i];

    LGMSVMomentCalculationApprox(dLambdaX, dAlpha[i], dLambdaEps[i], dRho[i],
                                 dSigma[i], dCoefVol, dCoefMeanRev1, dCoefCMS1,
                                 dt[i], FuncX, FuncX2, FuncV, FuncV2, FuncXV,
                                 lambdaArray);
    t1 = t2;
  }

  for (i = i; i <= endi; i++) {
    /* Precalculation on the option i*/

    /* First the time discretisation */
    if (i < endi) {
      t2 = dPWTime[i - do_switch];
    } else {
      t2 = dExTime;
    }

    /* Calculate constant values */
    temp_vol = dCoefVol * dSigma[i - do_switch];
    dt[i] = (t2 - t1);
    CoefIntT[i] = dLambdaEps[i - do_switch];
    Coef1ReT[i] = -0.5 * dAlpha[i - do_switch] * dAlpha[i - do_switch];
    Coef2ReT[i] = dLambdaEps[i - do_switch] +
                  dCoefMeanRev2 * dAlpha[i - do_switch] * dRho[i - do_switch] *
                      dSigma[i - do_switch];
    Coef2ImT[i] = -dAlpha[i - do_switch] * dRho[i - do_switch] * temp_vol;
    temp_vol *= temp_vol;
    Coef3ReT[i] = 0.5 * temp_vol;
    Coef3ImT[i] = 0.5 * temp_vol -
                  dCoefCMS2 * dSigma[i - do_switch] * dSigma[i - do_switch];

    LGMSVMomentCalculationApprox(dLambdaX, dAlpha[i - do_switch],
                                 dLambdaEps[i - do_switch], dRho[i - do_switch],
                                 dSigma[i - do_switch], dCoefVol, dCoefMeanRev2,
                                 dCoefCMS2, dt[i - do_switch], FuncX, FuncX2,
                                 FuncV, FuncV2, FuncXV, lambdaArray);
    t1 = t2;
  }

  dXMean = FuncX->dXt1;
  dXStd = sqrt(FuncX2->dXt1 - dXMean * dXMean);
  dVMean = FuncV->dXt1;
  dVStd = sqrt(FuncV2->dXt1 - dVMean * dVMean);

  if (NumerParams->iIntegMethod == 4) {
    /* Select the method */
    if (dXStd < NumerParams->dVolLimit) {
      /* special method for low volatilities */
      NumerParams->iIntegMethod = 3;
    } else {
      /* general method */
      NumerParams->iIntegMethod = 1;
    }
  }

  /* construct the grid */

  if (NumerParams->iIntegMethod == 0) {
    /* simple unifrom grid */

    if (NumerParams->dParam1 == 0) {
      NumerParams->dParam1 = 5;
    }

    if (NumerParams->dParam2 == 0) {
      NumerParams->dParam2 = 5;
    }

    XMax = (NumerParams->dParam1 + NumerParams->dParam2) * dXStd /
           (NumerParams->iNbX - 1);
    XMax = 1.0 / XMax;

    LGMSVCalculateGride(NumerParams->iNbX, XMax, X);
  } else if (NumerParams->iIntegMethod == 1) {
    /* Gauss Laguerre */

    W = dvector(0, NumerParams->iNbX - 1);

    if (!W) {
      err = "Memory Allocation failure (2) in LGMSVClosedFormNew";
      goto FREE_RETURN;
    }

    gaulag(X - 1, W - 1, NumerParams->iNbX, 0.0);
  } else if (NumerParams->iIntegMethod == 2) {
    /* Gauss Laguerre */
    W = dvector(0, NumerParams->iNbX - 1);

    if (!W) {
      err = "Memory Allocation failure (2) in LGMSVClosedFormNew";
      goto FREE_RETURN;
    }

    gauss_hermite(X - 1, W - 1, NumerParams->iNbX);
  } else if (NumerParams->iIntegMethod == 3) {
    /* Gauss Legendre */
    W = dvector(0, NumerParams->iNbX - 1);

    if (!W) {
      err = "Memory Allocation failure (2) in LGMSVClosedFormNew";
      goto FREE_RETURN;
    }

    dBreakPoints[0] = LGMSV_FIRSTBREAK;
    dBreakPoints[1] = NumerParams->dParam1;
    dBreakPoints[2] = NumerParams->dParam2;
    dBreakPoints[3] = LGMSV_LASTBREAK;

    iBreakNbX[0] = (int)(NumerParams->iNbX / 4);
    iBreakNbX[1] = iBreakNbX[0];
    iBreakNbX[2] = iBreakNbX[0];
    iBreakNbX[3] = NumerParams->iNbX - 3 * iBreakNbX[0];

    LGMSVConstructLegendreGrid(4, dBreakPoints, iBreakNbX, X, W);
  }

  /* calculate the integrand */
  if (NumerParams->iIntegMethod <= 1 || NumerParams->iIntegMethod == 3) {
    LGMSVSolveODEGride(NumerParams->iNbX, X, InitX, endi, dt, CoefIntT,
                       Coef1ReT, Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT,
                       Coef3ImT, NumerParams->dIntegParam, IntRe, IntIm);
  } else {
    LGMSVSolveODEGrideTV(NumerParams->iNbX, X, InitX, endi, dt, CoefIntT,
                         Coef1ReT, Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT,
                         Coef3ImT, NumerParams->dIntegParam, IntRe, IntIm);
  }

  /* do the integral */

  for (i = 0; i < iNbStrike; i++) {
    LogStrike[i] = log(fabs(dStrike[i]));
  }

  if (NumerParams->iIntegMethod == 0) {
    LGMSVComputeIntegral(NumerParams->iNbX, X, IntRe, IntIm, iNbStrike,
                         LogStrike, NumerParams->dIntegParam, Price);
  } else if (NumerParams->iIntegMethod == 1) {
    LGMSVComputeIntegralLaguerre(NumerParams->iNbX, X, W, 0.0, IntRe, IntIm,
                                 iNbStrike, LogStrike, NumerParams->dIntegParam,
                                 Price);
  } else if (NumerParams->iIntegMethod == 2) {
    LGMSVComputeIntegralHermite(NumerParams->iNbX, X, W, IntRe, IntIm,
                                iNbStrike, LogStrike, NumerParams->dIntegParam,
                                Price);

    for (i = 0; i < iNbStrike; i++) {
      if (LogStrike[i] < 0) {
        Price[i] += dFwd - dStrike[i];
      } else {
        Price[i] += dStrike[i] - dFwd;
      }
    }
  } else if (NumerParams->iIntegMethod == 3) {
    LGMSVComputeIntegralLegendre(NumerParams->iNbX, X, W, IntRe, IntIm,
                                 iNbStrike, LogStrike, NumerParams->dIntegParam,
                                 Price);
  }

FREE_RETURN:

  if (dt)
    free_dvector(dt, 0, endi);
  if (CoefIntT)
    free_dvector(CoefIntT, 0, endi);
  if (Coef1ReT)
    free_dvector(Coef1ReT, 0, endi);
  if (Coef2ReT)
    free_dvector(Coef2ReT, 0, endi);
  if (Coef2ImT)
    free_dvector(Coef2ImT, 0, endi);
  if (Coef3ReT)
    free_dvector(Coef3ReT, 0, endi);
  if (Coef3ImT)
    free_dvector(Coef3ImT, 0, endi);

  if (X)
    free_dvector(X, 0, NumerParams->iNbX - 1);
  if (W)
    free_dvector(W, 0, NumerParams->iNbX - 1);
  if (IntRe)
    free_dvector(IntRe, 0, NumerParams->iNbX - 1);
  if (IntIm)
    free_dvector(IntIm, 0, NumerParams->iNbX - 1);

  if (LogStrike)
    free_dvector(LogStrike, 0, iNbStrike - 1);
}

void LGMSVFwdClosedFormApprox(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    double dLambdaX,

    int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigma, double *dAlpha, double *dLambdaEps,
    double *dRho,

    double dTStar, /* Tstar in years from today */

    /* Product description */
    double dFwd, double dFixTime, double dCoefVol,

    double dSwitchTime, double dCoefMeanRev1, double dCoefCMS1,
    double dCoefMeanRev2, double dCoefCMS2,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Outputs */
    double *dFwdPrice) {
  Err err = NULL;
  double lambdaArray[10];

  double InitX, t1, t2, temp_vol;
  double dXMean, dXStd, dVMean, dVStd;
  int i, endi, iSwitchIndex, do_switch;

  LGMSVSolFunc FuncX_, *FuncX = &FuncX_, FuncX2_, *FuncX2 = &FuncX2_, FuncV_,
                       *FuncV = &FuncV_, FuncV2_, *FuncV2 = &FuncV2_, FuncXV_,
                       *FuncXV = &FuncXV_;

  double *dt = NULL, *CoefIntT = NULL, *Coef1ReT = NULL, *Coef2ReT = NULL,
         *Coef2ImT = NULL, *Coef3ReT = NULL, *Coef3ImT = NULL;

  double IntRe, IntIm;

  endi = Get_Index(dFixTime, dPWTime, iNbPWTime);
  iSwitchIndex = Get_Index(dSwitchTime, dPWTime, iNbPWTime);

  if (fabs(dSwitchTime - dPWTime[iSwitchIndex]) > 1.0E-08 &&
      dSwitchTime < dFixTime) {
    do_switch = 1;
  } else {
    do_switch = 0;

    if (dSwitchTime > dFixTime) {
      dSwitchTime = dFixTime;
      iSwitchIndex = endi;
    }
  }

  endi += do_switch;

  /* Memory allocation */

  dt = dvector(0, endi);
  CoefIntT = dvector(0, endi);
  Coef1ReT = dvector(0, endi);
  Coef2ReT = dvector(0, endi);
  Coef2ImT = dvector(0, endi);
  Coef3ReT = dvector(0, endi);
  Coef3ImT = dvector(0, endi);

  if (!dt || !CoefIntT || !Coef1ReT || !Coef2ReT || !Coef2ImT || !Coef3ReT ||
      !Coef3ImT) {
    err = "Memory Allocation failure (1) in LGMSVClosedFormNew";
    goto FREE_RETURN;
  }

  /* Precalculations */
  InitX = log(fabs(dFwd));
  LGMSVMomentInitApprox(InitX, 1.0, FuncX, FuncX2, FuncV, FuncV2, FuncXV);

  t1 = 0.0;
  for (i = 0; i <= iSwitchIndex; i++) {
    /* Precalculation on the option i*/

    /* First the time discretisation */
    if (i < iSwitchIndex) {
      t2 = dPWTime[i];
    } else {
      t2 = dSwitchTime;
    }

    /* Calculate constant values */
    temp_vol = dCoefVol * dSigma[i];
    dt[i] = (t2 - t1);
    CoefIntT[i] = dLambdaEps[i];
    Coef1ReT[i] = -0.5 * dAlpha[i] * dAlpha[i];
    Coef2ReT[i] =
        dLambdaEps[i] + dCoefMeanRev1 * dAlpha[i] * dRho[i] * dSigma[i];
    Coef2ImT[i] = -dAlpha[i] * dRho[i] * temp_vol;
    temp_vol *= temp_vol;
    Coef3ReT[i] = 0.5 * temp_vol;
    Coef3ImT[i] = 0.5 * temp_vol - dCoefCMS1 * dSigma[i] * dSigma[i];

    LGMSVMomentCalculationApprox(dLambdaX, dAlpha[i], dLambdaEps[i], dRho[i],
                                 dSigma[i], dCoefVol, dCoefMeanRev1, dCoefCMS1,
                                 dt[i], FuncX, FuncX2, FuncV, FuncV2, FuncXV,
                                 lambdaArray);
    t1 = t2;
  }

  for (i = i; i <= endi; i++) {
    /* Precalculation on the option i*/

    /* First the time discretisation */
    if (i < endi) {
      t2 = dPWTime[i - do_switch];
    } else {
      t2 = dFixTime;
    }

    /* Calculate constant values */
    temp_vol = dCoefVol * dSigma[i - do_switch];
    dt[i] = (t2 - t1);
    CoefIntT[i] = dLambdaEps[i - do_switch];
    Coef1ReT[i] = -0.5 * dAlpha[i - do_switch] * dAlpha[i - do_switch];
    Coef2ReT[i] = dLambdaEps[i - do_switch] +
                  dCoefMeanRev2 * dAlpha[i - do_switch] * dRho[i - do_switch] *
                      dSigma[i - do_switch];
    Coef2ImT[i] = -dAlpha[i - do_switch] * dRho[i - do_switch] * temp_vol;
    temp_vol *= temp_vol;
    Coef3ReT[i] = 0.5 * temp_vol;
    Coef3ImT[i] = 0.5 * temp_vol -
                  dCoefCMS2 * dSigma[i - do_switch] * dSigma[i - do_switch];

    LGMSVMomentCalculationApprox(dLambdaX, dAlpha[i - do_switch],
                                 dLambdaEps[i - do_switch], dRho[i - do_switch],
                                 dSigma[i - do_switch], dCoefVol, dCoefMeanRev2,
                                 dCoefCMS2, dt[i - do_switch], FuncX, FuncX2,
                                 FuncV, FuncV2, FuncXV, lambdaArray);
    t1 = t2;
  }

  dXMean = FuncX->dXt1;
  dXStd = sqrt(FuncX2->dXt1 - dXMean * dXMean);
  dVMean = FuncV->dXt1;
  dVStd = sqrt(FuncV2->dXt1 - dVMean * dVMean);

  /* calculate the integrand */

  LGMSVSolveODE(0.0, -1.0, InitX, endi, dt, CoefIntT, Coef1ReT, Coef1ReT,
                Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, &IntRe, &IntIm);

  *dFwdPrice = exp(IntRe);

  if (dCoefVol < 0.0) {
    (*dFwdPrice) *= -1.0;
  }

FREE_RETURN:

  if (dt)
    free_dvector(dt, 0, endi);
  if (CoefIntT)
    free_dvector(CoefIntT, 0, endi);
  if (Coef1ReT)
    free_dvector(Coef1ReT, 0, endi);
  if (Coef2ReT)
    free_dvector(Coef2ReT, 0, endi);
  if (Coef2ImT)
    free_dvector(Coef2ImT, 0, endi);
  if (Coef3ReT)
    free_dvector(Coef3ReT, 0, endi);
  if (Coef3ImT)
    free_dvector(Coef3ImT, 0, endi);
}

Err Construct_Swap_Schedule(long lTodayDate, char *cYcName, long lStartDate,
                            long lEndDate, SrtCompounding iFreq,
                            SrtBasisCode iBasis, int *iNbCoupon,
                            long **lCouponDate, double **dCouponTime,
                            double **dCouponCvg, double **dCouponDf) {
  long theo_date, act_date, temp_date;
  int i, ncpn;
  Err err = NULL;

  /*	1) First calculate the number of coupon */
  theo_date = lEndDate;
  act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
  ncpn = 1;

  while (act_date >= lStartDate - 5) {
    theo_date =
        add_unit(theo_date, -12 / iFreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
    act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    ncpn++;
  }
  ncpn--;

  if (ncpn < 2) {
    err = "Not enough coupons to construct a schedule";
    goto FREE_RETURN;
  }

  /*	2) Allocate memory */

  *lCouponDate = lvector(0, ncpn - 1);
  *dCouponTime = dvector(0, ncpn - 1);
  *dCouponCvg = dvector(0, ncpn - 1);
  *dCouponDf = dvector(0, ncpn - 1);

  if (!*lCouponDate || !*dCouponTime || !*dCouponCvg || !*dCouponDf) {
    err = "Memory allocation faillure in Construct_Swap_Schedule";
    goto FREE_RETURN;
  }

  /*	3) Fill the schedule informations */

  theo_date = lEndDate;
  act_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
  i = ncpn - 1;

  while (i >= 0) {
    (*dCouponTime)[i] = (act_date - lTodayDate) * YEARS_IN_DAY;
    (*lCouponDate)[i] = act_date;
    (*dCouponDf)[i] = swp_f_df(lTodayDate, act_date, cYcName);

    theo_date =
        add_unit(theo_date, -12 / iFreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

    temp_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    (*dCouponCvg)[i] = coverage(temp_date, act_date, iBasis);
    act_date = temp_date;

    i--;

    if (i == 0) {
      /* force the first coupon to be on start date */
      act_date = lStartDate;
    }
  }

  *iNbCoupon = ncpn;

FREE_RETURN:

  return err;
}

Err Calculate_HestonEquivalent_LGMSV(
    int iNbCoupon, long *lCouponDate, double *dCouponTime, double *dCouponCvg,
    double *dCouponDf, double dTPay, int iIsCMS, long lToday, char *cYcName,
    double dTau, double dTStar, double *dFwdSwapCash, double *dLevel,
    double *dShift, double *dCoefVol, double *dCoefMeanRev, double *dCoefCMS) {
  int i;
  double sumdf, level, swp_cash;
  double wi, xi, Li, theta, shift, temp;
  double expTi, expTi1;

  Err err = NULL;

  level = 0.0;
  sumdf = 0.0;

  for (i = 1; i < iNbCoupon; i++) {
    level += dCouponCvg[i] * dCouponDf[i];
    sumdf += dCouponDf[i];
  }

  swp_cash = (dCouponDf[0] - dCouponDf[iNbCoupon - 1]) / level;

  theta = level / sumdf;
  shift = 1.0 / theta;

  expTi = exp(-(dCouponTime[0] - dTStar) / dTau);

  *dCoefVol = 0.0;
  *dCoefMeanRev = 0.0;

  for (i = 1; i < iNbCoupon; i++) {
    wi = dCouponDf[i] / sumdf;
    Li = (dCouponDf[i - 1] - dCouponDf[i]) / (dCouponDf[i] * dCouponCvg[i]);
    xi = wi * (Li + shift) / (swp_cash + shift);
    expTi1 = exp(-(dCouponTime[i] - dTStar) / dTau);

    *dCoefVol += xi * (expTi - expTi1);
    *dCoefMeanRev += wi * (1.0 - expTi1);

    expTi = expTi1;
  }

  (*dCoefVol) *= dTau;
  (*dCoefMeanRev) *= dTau;

  if (iIsCMS) {
    temp = (1.0 - exp(-(dTPay - dTStar) / dTau)) * dTau;
    *dCoefCMS = (*dCoefMeanRev - temp) * (*dCoefVol);
    *dCoefMeanRev = temp;
  } else {
    *dCoefCMS = 0.0;
  }

  *dFwdSwapCash = swp_cash;
  *dLevel = level;
  *dShift = shift;

  return err;
}

Err Calculate_CoefShiftAndVol_LGMSV(/* product information */
                                    double dExeTime, int iNbCoupon,
                                    double *dCouponTime, double *dCouponCvg,
                                    double *dCouponDf,

                                    /* model information */
                                    double dTau,
                                    int iNbPWTime, /* Piece Wise Term Structures
                                                    */
                                    double *dPWTime, double *dSigmaTS,
                                    double dTStar,

                                    /* output */
                                    double *dLevel, double *dSwpCash,
                                    double *dShift, double *dCoefVol) {
  int i, index;
  double *cpn_G = NULL;

  double strike[3], price[3], vol[3];
  double sumx, sumy, sumxy, sumx2, coef_lin, const_lin;
  double t2, t1;
  double sumdf, swp_cash, level, std, cum_var, ex_zeta, ex_G;
  Err err = NULL;

  cpn_G = dvector(0, iNbCoupon - 1);

  if (!cpn_G) {
    err = "Memory allocation failure in Calculate_ShiftAndVol_LGMSV_FromLGM";
    goto FREE_RETURN;
  }

  cpn_G[0] = (1.0 - exp(-dCouponTime[0] / dTau)) * dTau;

  level = 0.0;
  sumdf = 0.0;

  for (i = 1; i < iNbCoupon; i++) {
    cpn_G[i] = (1.0 - exp(-dCouponTime[i] / dTau)) * dTau;
    level += dCouponDf[i] * dCouponCvg[i];
    sumdf += dCouponDf[i];
  }

  ex_G = (1.0 - exp(-dExeTime / dTau)) * dTau;

  swp_cash = (dCouponDf[0] - dCouponDf[iNbCoupon - 1]) / level;

  index = Get_Index(dExeTime, dPWTime, iNbPWTime);

  /* Compute Zeta */
  cum_var = 0.0;

  for (i = 0; i < index + 1; i++) {
    if (i > 0) {
      t1 = dPWTime[i - 1];
    } else {
      t1 = 0.0;
    }

    if (i == index || 0 == index) {
      /* Last part */
      t2 = dExeTime;
    } else {
      t2 = dPWTime[i];
    }

    cum_var += dSigmaTS[i] * dSigmaTS[i] * (t2 - t1);
  }

  ex_zeta = cum_var * exp(2.0 / dTau * dTStar);

  strike[1] = swp_cash;

  price[1] = lgmsopval1F(iNbCoupon, dCouponDf, dCouponCvg, cpn_G, ex_zeta, ex_G,
                         strike[1]);
  err = srt_f_optimpvol(price[1], swp_cash, strike[1], dExeTime, level, SRT_PUT,
                        SRT_NORMAL, &(vol[1]));

  if (err) {
    goto FREE_RETURN;
  }

  std = vol[1] * sqrt(dExeTime);
  strike[0] = strike[1] - 0.5 * std;
  strike[2] = strike[1] + 0.5 * std;

  price[0] = lgmsopval1F(iNbCoupon, dCouponDf, dCouponCvg, cpn_G, ex_zeta, ex_G,
                         strike[0]);
  err = srt_f_optimpvol(price[0], swp_cash, strike[0], dExeTime, level, SRT_PUT,
                        SRT_NORMAL, &(vol[0]));

  if (err) {
    goto FREE_RETURN;
  }

  price[2] = lgmsopval1F(iNbCoupon, dCouponDf, dCouponCvg, cpn_G, ex_zeta, ex_G,
                         strike[2]);
  err = srt_f_optimpvol(price[2], swp_cash, strike[2], dExeTime, level, SRT_PUT,
                        SRT_NORMAL, &(vol[2]));

  if (err) {
    goto FREE_RETURN;
  }

  /* Linear regression */
  sumx = 0.0;
  sumx2 = 0.0;
  sumy = 0.0;
  sumxy = 0.0;

  for (i = 0; i < 3; i++) {
    sumx += strike[i];
    sumx2 += strike[i] * strike[i];
    sumy += vol[i];
    sumxy += strike[i] * vol[i];
  }

  sumx /= 3.0;
  sumx2 /= 3.0;
  sumy /= 3.0;
  sumxy /= 3.0;

  coef_lin = (sumxy - sumx * sumy) / (sumx2 - sumx * sumx);
  const_lin = sumy - coef_lin * sumx;

  (*dShift) = (const_lin - coef_lin * swp_cash) / (2.0 * coef_lin);

  /* Min  , max on the shift to avoid numerical problems...*/
  (*dShift) = min(MAX_LGMSHIFT, max(-MAX_LGMSHIFT, (*dShift)));

  err = srt_f_optimpvol(price[1], fabs(swp_cash + (*dShift)),
                        fabs(strike[1] + (*dShift)), dExeTime, level, SRT_PUT,
                        SRT_LOGNORMAL, dCoefVol);

  if (err) {
    goto FREE_RETURN;
  }

  (*dCoefVol) = (*dCoefVol) / sqrt(cum_var / dExeTime);

  if (swp_cash + (*dShift) < 0) {
    (*dCoefVol) *= -1.0;
  }

  if (dSwpCash) {
    (*dSwpCash) = swp_cash;
  }
  if (dLevel) {
    (*dLevel) = level;
  }

FREE_RETURN:

  if (cpn_G)
    free_dvector(cpn_G, 0, iNbCoupon - 1);

  return err;
}

void Calculate_CoefCMSandMeanRev_LGMSV(/* product information */
                                       int iNbCoupon, double *dCouponTime,
                                       double *dCouponCvg, double *dCouponDf,
                                       double dTPay, int iIsCMS,
                                       double dTRStart, double dTREnd,
                                       int iIsRatioNum,

                                       /* model information */
                                       double dTau, double dTStar,
                                       double dCoefVol,

                                       /* output */
                                       double *dCoefMeanRev, double *dCoefCMS) {
  double temp, sumdf, expTi1, coefMean, coefCMS;
  int i;

  coefMean = 0.0;
  sumdf = 0.0;

  for (i = 1; i < iNbCoupon; i++) {
    sumdf += dCouponDf[i];
    expTi1 = exp(-(dCouponTime[i] - dTStar) / dTau);
    coefMean += dCouponDf[i] * (1.0 - expTi1);
  }

  coefMean /= sumdf;
  coefMean *= dTau;

  if (iIsCMS || iIsRatioNum) {
    temp = (1.0 - exp(-(dTPay - dTStar) / dTau)) * dTau;
    coefCMS = (coefMean - temp);
    coefMean = temp;

    if (iIsRatioNum) {
      temp =
          (exp(-(dTREnd - dTStar) / dTau) - exp(-(dTRStart - dTStar) / dTau)) *
          dTau;
      coefCMS -= temp;
      coefMean += temp;
    }
  } else {
    coefCMS = 0.0;
  }

  *dCoefMeanRev = coefMean;
  *dCoefCMS = coefCMS * dCoefVol;
}

Err Calculate_HestonEquivalent_LGMSV_FromLGM(/*	Market information */
                                             long lToday, char *cYcName,

                                             /* product information */
                                             long lExeDate, double dExeTime,
                                             int iNbCoupon, long *lCouponDate,
                                             double *dCouponTime,
                                             double *dCouponCvg,
                                             double *dCouponDf, double dTPay,
                                             int iIsCMS, double dTRStart,
                                             double dTREnd, int iIsRatioNum,

                                             /* model information */
                                             double dTau,
                                             int iNbPWTime, /* Piece Wise Term
                                                               Structures  */
                                             double *dPWTime, double *dSigmaTS,
                                             double dTStar,

                                             double *dFwdSwapCash,
                                             double *dLevel, double *dShift,
                                             double *dCoefVol,
                                             double *dCoefMeanRev,
                                             double *dCoefCMS) {
  int i, index;
  double *cpn_G = NULL;

  double strike[3], price[3], vol[3];
  double sumx, sumy, sumxy, sumx2, coef_lin, const_lin;
  double t2, t1;
  double swp_cash, level, std, cum_var, ex_zeta, ex_G;
  double sumdf, expTi1, wi, temp;
  Err err = NULL;

  cpn_G = dvector(0, iNbCoupon - 1);

  if (!cpn_G) {
    err = "Memory allocation failure in Calculate_ShiftedLog_LGM";
    goto FREE_RETURN;
  }

  cpn_G[0] = (1.0 - exp(-dCouponTime[0] / dTau)) * dTau;

  level = 0.0;
  sumdf = 0.0;

  for (i = 1; i < iNbCoupon; i++) {
    cpn_G[i] = (1.0 - exp(-dCouponTime[i] / dTau)) * dTau;
    level += dCouponDf[i] * dCouponCvg[i];
    sumdf += dCouponDf[i];
  }

  ex_G = (1.0 - exp(-dExeTime / dTau)) * dTau;

  swp_cash = (dCouponDf[0] - dCouponDf[iNbCoupon - 1]) / level;

  index = Get_Index(dExeTime, dPWTime, iNbPWTime);

  /* Compute Zeta */
  cum_var = 0.0;

  for (i = 0; i < index + 1; i++) {
    if (i > 0) {
      t1 = dPWTime[i - 1];
    } else {
      t1 = 0.0;
    }

    if (i == index || 0 == index) {
      /* Last part */
      t2 = dExeTime;
    } else {
      t2 = dPWTime[i];
    }

    cum_var += dSigmaTS[i] * dSigmaTS[i] * (t2 - t1);
  }

  ex_zeta = cum_var * exp(2.0 / dTau * dTStar);

  strike[1] = swp_cash;

  price[1] = lgmsopval1F(iNbCoupon, dCouponDf, dCouponCvg, cpn_G, ex_zeta, ex_G,
                         strike[1]);

  err = srt_f_optimpvol(price[1], swp_cash, strike[1], dExeTime, level, SRT_PUT,
                        SRT_NORMAL, &(vol[1]));

  if (err) {
    goto FREE_RETURN;
  }

  std = vol[1] * sqrt(dExeTime);
  strike[0] = strike[1] - 0.5 * std;
  strike[2] = strike[1] + 0.5 * std;

  price[0] = lgmsopval1F(iNbCoupon, dCouponDf, dCouponCvg, cpn_G, ex_zeta, ex_G,
                         strike[0]);

  err = srt_f_optimpvol(price[0], swp_cash, strike[0], dExeTime, level, SRT_PUT,
                        SRT_NORMAL, &(vol[0]));

  if (err) {
    goto FREE_RETURN;
  }

  price[2] = lgmsopval1F(iNbCoupon, dCouponDf, dCouponCvg, cpn_G, ex_zeta, ex_G,
                         strike[2]);

  err = srt_f_optimpvol(price[2], swp_cash, strike[2], dExeTime, level, SRT_PUT,
                        SRT_NORMAL, &(vol[2]));

  if (err) {
    goto FREE_RETURN;
  }

  /* Linear regression */
  sumx = 0.0;
  sumx2 = 0.0;
  sumy = 0.0;
  sumxy = 0.0;

  for (i = 0; i < 3; i++) {
    sumx += strike[i];
    sumx2 += strike[i] * strike[i];
    sumy += vol[i];
    sumxy += strike[i] * vol[i];
  }

  sumx /= 3.0;
  sumx2 /= 3.0;
  sumy /= 3.0;
  sumxy /= 3.0;

  coef_lin = (sumxy - sumx * sumy) / (sumx2 - sumx * sumx);
  const_lin = sumy - coef_lin * sumx;

  (*dShift) = (const_lin - coef_lin * swp_cash) / (2.0 * coef_lin);

  /* Min  , max on the shift to avoid numerical problems...*/
  (*dShift) = min(MAX_LGMSHIFT, max(-MAX_LGMSHIFT, (*dShift)));

  err = srt_f_optimpvol(price[1], swp_cash + (*dShift), strike[1] + (*dShift),
                        dExeTime, level, SRT_PUT, SRT_LOGNORMAL, dCoefVol);

  if (err) {
    goto FREE_RETURN;
  }

  (*dCoefVol) = (*dCoefVol) / sqrt(cum_var / dExeTime);

  *dCoefMeanRev = 0.0;

  for (i = 1; i < iNbCoupon; i++) {
    wi = dCouponDf[i] / sumdf;
    expTi1 = exp(-(dCouponTime[i] - dTStar) / dTau);
    *dCoefMeanRev += wi * (1.0 - expTi1);
  }

  (*dCoefMeanRev) *= dTau;

  if (iIsCMS || iIsRatioNum) {
    temp = (1.0 - exp(-(dTPay - dTStar) / dTau)) * dTau;
    *dCoefCMS = (*dCoefMeanRev - temp);
    *dCoefMeanRev = temp;

    if (iIsRatioNum) {
      temp =
          (exp(-(dTREnd - dTStar) / dTau) - exp(-(dTRStart - dTStar) / dTau)) *
          dTau;
      *dCoefCMS -= temp;
      *dCoefMeanRev += temp;
    }

    *dCoefCMS *= (*dCoefVol);
  } else {
    *dCoefCMS = 0.0;
  }

  *dFwdSwapCash = swp_cash;
  *dLevel = level;

FREE_RETURN:

  if (cpn_G)
    free_dvector(cpn_G, 0, iNbCoupon - 1);

  return err;
}

Err LGMSVCmsRateApprox_TS(

    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lFixingDate, long lStartDate, long lEndDate, long lPayDate,

    double dLambdaX, int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigmaTS, double *dAlphaTS, double *dLambdaEpsTS,
    double *dRhoTS, double dTStar,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *dFwd) {
  int ncpn;
  SrtCompounding ifreq;
  SrtBasisCode ibasis;
  long *cpn_date = NULL;
  double *cpn_time = NULL, *cpn_cvg = NULL, *cpn_df = NULL;

  long today;

  double dFixTime, dPayTime;
  double swp_cash, swp_rate, spread, level;

  SrtCurvePtr yc_ptr;
  Err err = NULL;

  double shift, coefVol, coefMeanRev, coefCMS;
  double dTau, dShiftFwd;
  double *dShiftStrike = NULL;

  /*	1) get the schedule informations */
  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }

  today = get_today_from_curve(yc_ptr);

  err = interp_compounding(swaption_freq, &ifreq);
  if (err) {
    goto FREE_RETURN;
  }

  err = interp_basis(swaption_basis, &ibasis);
  if (err) {
    goto FREE_RETURN;
  }

  /*	2) Construct the schedule */

  err = Construct_Swap_Schedule(today, yc_name, lStartDate, lEndDate, ifreq,
                                ibasis, &ncpn, &cpn_date, &cpn_time, &cpn_cvg,
                                &cpn_df);

  if (err) {
    goto FREE_RETURN;
  }

  /*	3) Calculate the constant for equivalent Heston */

  dTau = 1.0 / dLambdaX;
  dFixTime = (lFixingDate - today) * YEARS_IN_DAY;
  dPayTime = (lPayDate - today) * YEARS_IN_DAY;

  if (NumerParams->iCalibLGM) {
    err = Calculate_CoefShiftAndVol_LGMSV(
        dFixTime, ncpn, cpn_time, cpn_cvg, cpn_df, dTau, iNbPWTime, dPWTime,
        dSigmaTS, dTStar, &level, &swp_cash, &shift, &coefVol);

    Calculate_CoefCMSandMeanRev_LGMSV(ncpn, cpn_time, cpn_cvg, cpn_df, dPayTime,
                                      1, 0, 0, 0, dTau, dTStar, coefVol,
                                      &coefMeanRev, &coefCMS);
  } else {
    err = Calculate_HestonEquivalent_LGMSV(
        ncpn, cpn_date, cpn_time, cpn_cvg, cpn_df, dPayTime, 1, today, yc_name,
        dTau, dTStar, &swp_cash, &level, &shift, &coefVol, &coefMeanRev,
        &coefCMS);
  }

  if (err) {
    goto FREE_RETURN;
  }

  /*	4) shift the parameters */

  err = swp_f_ForwardRate(lStartDate, lEndDate, swaption_freq, swaption_basis,
                          yc_name, ref_rate_name, &swp_rate);

  if (err) {
    goto FREE_RETURN;
  }

  spread = swp_rate - swp_cash;
  shift -= spread;

  dShiftFwd = swp_rate + shift;

  LGMSVFwdClosedFormApprox(1.0 / dTau, iNbPWTime, /* Term Structure of g(t) */
                           dPWTime, dSigmaTS, dAlphaTS, /* Alpha of V = Eps^2 */
                           dLambdaEpsTS,   /* LambdaEps of V = Eps^2 */
                           dRhoTS, dTStar, /* Tstar in years from today */
                           dShiftFwd, dFixTime, coefVol, 1000, coefMeanRev,
                           coefCMS, 0, 0, NumerParams, dFwd);

  *dFwd -= shift;

FREE_RETURN:

  if (cpn_date)
    free_lvector(cpn_date, 0, ncpn - 1);
  if (cpn_time)
    free_dvector(cpn_time, 0, ncpn - 1);
  if (cpn_cvg)
    free_dvector(cpn_cvg, 0, ncpn - 1);

  return err;
}

Err LGMSVOptionApprox_TS(

    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lEndDate, long lPayDate, int iIsCMS,
    int iNbStrike, double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    double dLambdaX, int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigmaTS, double *dAlphaTS, double *dLambdaEpsTS,
    double *dRhoTS, double dTStar,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice) {
  int i, ncpn;
  SrtCompounding ifreq;
  SrtBasisCode ibasis;
  long *cpn_date = NULL;
  double *cpn_time = NULL, *cpn_cvg = NULL, *cpn_df = NULL;

  long today;

  double dExTime, dPayTime;
  double swp_cash, swp_rate, spread, level, numeraire;

  SrtCurvePtr yc_ptr;
  Err err = NULL;

  double shift, coefVol, coefMeanRev, coefCMS;
  double dTau, dShiftFwd, RescaleFactor;
  double *dShiftStrike = NULL;

  /*	1) get the schedule informations */
  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }

  today = get_today_from_curve(yc_ptr);

  err = interp_compounding(swaption_freq, &ifreq);
  if (err) {
    goto FREE_RETURN;
  }

  err = interp_basis(swaption_basis, &ibasis);
  if (err) {
    goto FREE_RETURN;
  }

  /*	2) Construct the schedule */

  err = Construct_Swap_Schedule(today, yc_name, lStartDate, lEndDate, ifreq,
                                ibasis, &ncpn, &cpn_date, &cpn_time, &cpn_cvg,
                                &cpn_df);

  if (err) {
    goto FREE_RETURN;
  }

  /*	3) Calculate the constant for equivalent Heston */

  dTau = 1.0 / dLambdaX;
  dPayTime = (lPayDate - today) * YEARS_IN_DAY;
  dExTime = (lExDate - today) * YEARS_IN_DAY;

  if (NumerParams->iCalibLGM) {
    err = Calculate_CoefShiftAndVol_LGMSV(
        dExTime, ncpn, cpn_time, cpn_cvg, cpn_df, dTau, iNbPWTime, dPWTime,
        dSigmaTS, dTStar, &level, &swp_cash, &shift, &coefVol);

    Calculate_CoefCMSandMeanRev_LGMSV(ncpn, cpn_time, cpn_cvg, cpn_df, dPayTime,
                                      iIsCMS, 0, 0, 0, dTau, dTStar, coefVol,
                                      &coefMeanRev, &coefCMS);
  } else {
    err = Calculate_HestonEquivalent_LGMSV(
        ncpn, cpn_date, cpn_time, cpn_cvg, cpn_df, dPayTime, iIsCMS, today,
        yc_name, dTau, dTStar, &swp_cash, &level, &shift, &coefVol,
        &coefMeanRev, &coefCMS);
  }

  if (err) {
    goto FREE_RETURN;
  }

  /*	4) shift the parameters */

  err = swp_f_ForwardRate(lStartDate, lEndDate, swaption_freq, swaption_basis,
                          yc_name, ref_rate_name, &swp_rate);

  if (err) {
    goto FREE_RETURN;
  }

  spread = swp_rate - swp_cash;
  shift -= spread;

  dShiftStrike = dvector(0, iNbStrike - 1);

  if (!dShiftStrike) {
    err = "Memory allocation faillure in LGMSVOptionApprox";
    goto FREE_RETURN;
  }

  dShiftFwd = swp_rate + shift;

  for (i = 0; i < iNbStrike; i++) {
    dShiftStrike[i] = dStrike[i] + shift;
  }

  RescaleFactor = 1.0;

  /* rescale */

  for (i = 0; i < iNbStrike; i++) {
    dShiftStrike[i] /= RescaleFactor;
  }

  dShiftFwd /= RescaleFactor;

  /*	5) Call the solver */

  if (NumerParams->iIntegMethod != 2) {
    if (pay_rec == -1) {
      NumerParams->dIntegParam = -1 - NumerParams->dIntegParam;
    }
  }

  LGMSVClosedFormApprox(1.0 / dTau, iNbPWTime, /* Term Structure of g(t) */
                        dPWTime, dSigmaTS, dAlphaTS, /* Alpha of V = Eps^2 */
                        dLambdaEpsTS,   /* LambdaEps of V = Eps^2 */
                        dRhoTS, dTStar, /* Tstar in years from today */
                        dShiftFwd, iNbStrike, dShiftStrike, dExTime, coefVol,
                        1000, coefMeanRev, coefCMS, 0, 0, NumerParams,
                        pSwaptionPrice);

  if (iIsCMS) {
    /* change the numeraire */
    numeraire = RescaleFactor * swp_f_df(today, lPayDate, yc_name);
  } else {
    numeraire = RescaleFactor * level;
  }

  if (NumerParams->iIntegMethod == 2) {
    dShiftFwd = swp_rate + shift;

    for (i = 0; i < iNbStrike; i++) {
      if (1.0 > dShiftStrike[i]) {
        pSwaptionPrice[i] *= dShiftFwd;
      } else {
        pSwaptionPrice[i] =
            (pSwaptionPrice[i] + dShiftStrike[i] - 1.0) * dShiftFwd;
      }
    }

    for (i = 0; i < iNbStrike; i++) {
      if (pay_rec == 1) {
        pSwaptionPrice[i] *= numeraire;
      } else {
        /* call put parity */
        pSwaptionPrice[i] =
            (pSwaptionPrice[i] + dStrike[i] - swp_rate) * numeraire;
      }
    }
  } else {
    for (i = 0; i < iNbStrike; i++) {
      pSwaptionPrice[i] *= numeraire;
    }
  }

FREE_RETURN:

  if (cpn_date)
    free_lvector(cpn_date, 0, ncpn - 1);
  if (cpn_time)
    free_dvector(cpn_time, 0, ncpn - 1);
  if (cpn_cvg)
    free_dvector(cpn_cvg, 0, ncpn - 1);

  if (dShiftStrike)
    free_dvector(dShiftStrike, 0, iNbStrike - 1);

  return err;
}

Err LGMSVLiborOptionApprox_TS(

    char *yc_name,                /*	Name of the yield curve */
    char *swaption_ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis, char *libor_ref_rate_name, char *libor_freq,
    char *libor_basis,

    long lExDate, long lStartDate, long lEndDate, long lPayDate, int iIsCMS,
    long lLiborFixDate, long lLiborStartDate, long lLiborEndDate,

    int iNbStrike, double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    double dLambdaX, int iNbPWTime, /* Piece Wise Term Structures  */
    double *dPWTime, double *dSigmaTS, double *dAlphaTS, double *dLambdaEpsTS,
    double *dRhoTS, double dTStar,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice) {
  int i, ncpn;
  SrtCompounding ifreq, ifreq2;
  SrtBasisCode ibasis, ibasis2;
  long *cpn_date = NULL;
  double *cpn_time = NULL, *cpn_cvg = NULL, *cpn_df = NULL, *dPriceAdj = NULL,
         *dPriceNonAdj = NULL;

  long today;

  double dExTime, dPayTime, dFixTime, dTRStart, dTREnd;
  double dLiborFwd, theta;
  double swp_cash, swp_rate, spread, level, numeraire;

  SrtCurvePtr yc_ptr;
  Err err = NULL;

  double shift, coefVol, coefMeanRev1, coefCMS1, coefMeanRev2, coefCMS2;
  double dTau, dShiftFwd;
  double *dShiftStrike = NULL;

  /*	1) get the schedule informations */
  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found";
    goto FREE_RETURN;
  }

  today = get_today_from_curve(yc_ptr);

  err = interp_compounding(swaption_freq, &ifreq);
  if (err) {
    goto FREE_RETURN;
  }

  err = interp_basis(swaption_basis, &ibasis);
  if (err) {
    goto FREE_RETURN;
  }

  err = interp_compounding(libor_freq, &ifreq2);
  if (err) {
    goto FREE_RETURN;
  }

  err = interp_basis(libor_basis, &ibasis2);
  if (err) {
    goto FREE_RETURN;
  }

  theta = coverage(lLiborStartDate, lLiborEndDate, ibasis2);

  /*	2) Construct the schedule */

  err = Construct_Swap_Schedule(today, yc_name, lStartDate, lEndDate, ifreq,
                                ibasis, &ncpn, &cpn_date, &cpn_time, &cpn_cvg,
                                &cpn_df);

  if (err) {
    goto FREE_RETURN;
  }

  /*	3) Calculate the constant for equivalent Heston */

  dTau = 1.0 / dLambdaX;
  dPayTime = (lPayDate - today) * YEARS_IN_DAY;
  dExTime = (lExDate - today) * YEARS_IN_DAY;
  dFixTime = (lLiborFixDate - today) * YEARS_IN_DAY;
  dTRStart = (lLiborStartDate - today) * YEARS_IN_DAY;
  dTREnd = (lLiborEndDate - today) * YEARS_IN_DAY;

  err = Calculate_CoefShiftAndVol_LGMSV(
      dExTime, ncpn, cpn_time, cpn_cvg, cpn_df, dTau, iNbPWTime, dPWTime,
      dSigmaTS, dTStar, &level, &swp_cash, &shift, &coefVol);

  if (err) {
    goto FREE_RETURN;
  }

  Calculate_CoefCMSandMeanRev_LGMSV(ncpn, cpn_time, cpn_cvg, cpn_df, dPayTime,
                                    iIsCMS, 0, 0, 0, dTau, dTStar, coefVol,
                                    &coefMeanRev2, &coefCMS2);

  Calculate_CoefCMSandMeanRev_LGMSV(ncpn, cpn_time, cpn_cvg, cpn_df, dPayTime,
                                    iIsCMS, dTRStart, dTREnd, 1, dTau, dTStar,
                                    coefVol, &coefMeanRev1, &coefCMS1);

  /*	4) shift the parameters */

  err = swp_f_ForwardRate(lStartDate, lEndDate, swaption_freq, swaption_basis,
                          yc_name, swaption_ref_rate_name, &swp_rate);

  if (err) {
    goto FREE_RETURN;
  }

  spread = swp_rate - swp_cash;
  shift -= spread;

  dShiftStrike = dvector(0, iNbStrike - 1);
  dPriceAdj = dvector(0, iNbStrike - 1);
  dPriceNonAdj = dvector(0, iNbStrike - 1);

  if (!dShiftStrike || !dPriceAdj || !dPriceNonAdj) {
    err = "Memory allocation faillure in LGMSVOptionApprox";
    goto FREE_RETURN;
  }

  dShiftFwd = swp_rate + shift;

  for (i = 0; i < iNbStrike; i++) {
    dShiftStrike[i] = dStrike[i] + shift;
  }

  if (NumerParams->iIntegMethod == 2) {
    /* rescale */

    for (i = 0; i < iNbStrike; i++) {
      dShiftStrike[i] /= dShiftFwd;
    }

    dShiftFwd = 1.0;
  }

  /*	5) Call the solver */

  if (NumerParams->iIntegMethod == 0 || NumerParams->iIntegMethod == 1) {
    if (pay_rec == -1) {
      NumerParams->dIntegParam = -1 - NumerParams->dIntegParam;
    }
  }

  /* First Calculate the expectation */
  err = LGMSVCmsRateApprox_TS(yc_name, libor_ref_rate_name, libor_freq,
                              libor_basis, lLiborFixDate, lLiborStartDate,
                              lLiborEndDate, lPayDate, 1.0 / dTau, iNbPWTime,
                              dPWTime, dSigmaTS, dAlphaTS, dLambdaEpsTS, dRhoTS,
                              dTStar, NumerParams, &dLiborFwd);

  if (err) {
    goto FREE_RETURN;
  }

  /* Then calculate the non adjusted part */
  LGMSVClosedFormApprox(1.0 / dTau, iNbPWTime, /* Term Structure of g(t) */
                        dPWTime, dSigmaTS, dAlphaTS, /* Alpha of V = Eps^2 */
                        dLambdaEpsTS,   /* LambdaEps of V = Eps^2 */
                        dRhoTS, dTStar, /* Tstar in years from today */
                        dShiftFwd, iNbStrike, dShiftStrike, dExTime, coefVol,
                        1000, coefMeanRev2, coefCMS2, 0, 0, NumerParams,
                        dPriceNonAdj);

  /* Then calculate the adjusted part */
  LGMSVClosedFormApprox(1.0 / dTau, iNbPWTime, /* Term Structure of g(t) */
                        dPWTime, dSigmaTS, dAlphaTS, /* Alpha of V = Eps^2 */
                        dLambdaEpsTS,   /* LambdaEps of V = Eps^2 */
                        dRhoTS, dTStar, /* Tstar in years from today */
                        dShiftFwd, iNbStrike, dShiftStrike, dExTime, coefVol,
                        dFixTime, coefMeanRev1, coefCMS1, coefMeanRev2,
                        coefCMS2, NumerParams, dPriceAdj);

  if (iIsCMS) {
    /* change the numeraire */
    numeraire = swp_f_df(today, lPayDate, yc_name);
  } else {
    numeraire = level;
  }

  if (NumerParams->iIntegMethod == 2) {
    dShiftFwd = swp_rate + shift;

    for (i = 0; i < iNbStrike; i++) {
      if (1.0 > dShiftStrike[i]) {
        dPriceNonAdj[i] *= dShiftFwd;
        dPriceAdj[i] *= dShiftFwd;
      } else {
        dPriceNonAdj[i] = (dPriceNonAdj[i] + dShiftStrike[i] - 1.0) * dShiftFwd;
        dPriceAdj[i] = (dPriceNonAdj[i] + dShiftStrike[i] - 1.0) * dShiftFwd;
      }
    }

    for (i = 0; i < iNbStrike; i++) {
      if (pay_rec == 1) {
        dPriceNonAdj[i] *= numeraire;
        dPriceAdj[i] *= numeraire;
      } else {
        /* call put parity */
        dPriceNonAdj[i] = (dPriceNonAdj[i] + dStrike[i] - swp_rate) * numeraire;
        dPriceAdj[i] = (dPriceAdj[i] + dStrike[i] - swp_rate) * numeraire;
      }
    }
  } else {
    for (i = 0; i < iNbStrike; i++) {
      dPriceNonAdj[i] *= numeraire;
      dPriceAdj[i] *= numeraire * (dLiborFwd * theta + 1.0);
      pSwaptionPrice[i] = 1.0 / theta * (dPriceAdj[i] - dPriceNonAdj[i]);
    }
  }

FREE_RETURN:

  if (cpn_date)
    free_lvector(cpn_date, 0, ncpn - 1);
  if (cpn_time)
    free_dvector(cpn_time, 0, ncpn - 1);
  if (cpn_cvg)
    free_dvector(cpn_cvg, 0, ncpn - 1);

  if (dShiftStrike)
    free_dvector(dShiftStrike, 0, iNbStrike - 1);
  if (dPriceAdj)
    free_dvector(dPriceAdj, 0, iNbStrike - 1);
  if (dPriceNonAdj)
    free_dvector(dPriceNonAdj, 0, iNbStrike - 1);

  return err;
}

Err LGMSVCmsRateApprox(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lFixDate, long lStartDate, long lEndDate, long lPayDate,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *dFwd) {
  Err err = NULL;

  int dNbSig;
  double *dSigTime = NULL, *dSig = NULL, *dAlphaTS = NULL, *dRhoTS = NULL,
         *dLambdaEpsTS = NULL;

  double dTStar;
  double dTau;
  int one2F;

  /* Get the TS */

  err = Get_LGMSV_TermStructure(und_name, &dSigTime, &dSig, &dAlphaTS, &dRhoTS,
                                &dLambdaEpsTS, &dTStar, &dNbSig, &dTau, &one2F,
                                NULL, NULL, NULL);

  if (err) {
    goto FREE_RETURN;
  }

  err = LGMSVCmsRateApprox_TS(
      yc_name, ref_rate_name, swaption_freq, swaption_basis, lFixDate,
      lStartDate, lEndDate, lPayDate, 1.0 / dTau, dNbSig, dSigTime, dSig,
      dAlphaTS, dLambdaEpsTS, dRhoTS, dTStar, NumerParams, dFwd);

FREE_RETURN:

  if (dSigTime)
    free(dSigTime);
  if (dSig)
    free(dSig);

  if (dAlphaTS)
    free(dAlphaTS);
  if (dRhoTS)
    free(dRhoTS);
  if (dLambdaEpsTS)
    free(dLambdaEpsTS);

  return err;
}

Err LGMSVOptionApprox(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lEndDate, long lPayDate, int IsCMS,
    long lRStartDate, long lREndDate, int iIsRatioNum, int iNbStrike,
    double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice) {
  Err err = NULL;

  int dNbSig;
  double *dSigTime = NULL, *dSig = NULL, *dAlphaTS = NULL, *dRhoTS = NULL,
         *dLambdaEpsTS = NULL;

  double dTStar;

  double dTau;
  int one2F;

  /* Get the TS */

  err = Get_LGMSV_TermStructure(und_name, &dSigTime, &dSig, &dAlphaTS, &dRhoTS,
                                &dLambdaEpsTS, &dTStar, &dNbSig, &dTau, &one2F,
                                NULL, NULL, NULL);

  if (err) {
    goto FREE_RETURN;
  }

  err = LGMSVOptionApprox_TS(
      yc_name, ref_rate_name, swaption_freq, swaption_basis, lExDate,
      lStartDate, lEndDate, lPayDate, IsCMS, iNbStrike, dStrike, pay_rec,
      1.0 / dTau, dNbSig, dSigTime, dSig, dAlphaTS, dLambdaEpsTS, dRhoTS,
      dTStar, NumerParams, pSwaptionPrice);

FREE_RETURN:

  if (dSigTime)
    free(dSigTime);
  if (dSig)
    free(dSig);

  if (dAlphaTS)
    free(dAlphaTS);
  if (dRhoTS)
    free(dRhoTS);
  if (dLambdaEpsTS)
    free(dLambdaEpsTS);

  return err;
}

Err LGMSVLiborOptionApprox(
    char *und_name,               /*	Name of the underlying */
    char *yc_name,                /*	Name of the yield curve */
    char *swaption_ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,
    char *libor_ref_rate_name, /*	Name of the reference rate */
    char *libor_freq,          /*	Frequency and basis of underlying swaptions */
    char *libor_basis,

    long lExDate, long lStartDate, long lEndDate, long lPayDate, int IsCMS,
    long lLiborFixDate, long lLiborStartDate, long lLiborEndDate, int iNbStrike,
    double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice) {
  Err err = NULL;

  int dNbSig;
  double *dSigTime = NULL, *dSig = NULL, *dAlphaTS = NULL, *dRhoTS = NULL,
         *dLambdaEpsTS = NULL;

  double dTStar;

  double dTau;
  int one2F;

  /* Get the TS */

  err = Get_LGMSV_TermStructure(und_name, &dSigTime, &dSig, &dAlphaTS, &dRhoTS,
                                &dLambdaEpsTS, &dTStar, &dNbSig, &dTau, &one2F,
                                NULL, NULL, NULL);

  if (err) {
    goto FREE_RETURN;
  }

  err = LGMSVLiborOptionApprox_TS(
      yc_name, swaption_ref_rate_name, swaption_freq, swaption_basis,
      libor_ref_rate_name, libor_freq, libor_basis, lExDate, lStartDate,
      lEndDate, lPayDate, IsCMS, lLiborFixDate, lLiborStartDate, lLiborEndDate,
      iNbStrike, dStrike, pay_rec, 1.0 / dTau, dNbSig, dSigTime, dSig, dAlphaTS,
      dLambdaEpsTS, dRhoTS, dTStar, NumerParams, pSwaptionPrice);

FREE_RETURN:

  if (dSigTime)
    free(dSigTime);
  if (dSig)
    free(dSig);

  if (dAlphaTS)
    free(dAlphaTS);
  if (dRhoTS)
    free(dRhoTS);
  if (dLambdaEpsTS)
    free(dLambdaEpsTS);

  return err;
}

Err Calculate_ShiftAndVol_LGMSV_FromLGM_struct(/* product information */
                                               CALIBINSTRUMENTDLM CalibInst,
                                               CALIBCPNSCHEDULEDLM
                                                   CalibCpnSchedule,

                                               /* model information */
                                               LGMSV_MODEL model,

                                               /* output */
                                               LGMSV_HESTONINST HestonParam) {
  int i, index;

  double strike[3], price[3], vol[3];
  double price1, price2, vol1, vol2, ratio2F;
  double sumx, sumy, sumxy, sumx2, coef_lin, const_lin;
  double t2, t1;
  double swp_cash, std, cum_var, cum_var2, cum_var12;
  double ex_zeta1, ex_zeta2, ex_zeta12;
  Err err = NULL;

  /*	first calculate the atm normal std */
  if (CalibInst->dExeTime > 1.0E-08) {
    swp_cash = CalibInst->dFwdCash;

    index = Get_Index(CalibInst->dExeTime, model->dPWTime, model->iNbPWTime);

    /* Compute Zeta */
    cum_var = 0.0;

    for (i = 0; i < index + 1; i++) {
      if (i > 0) {
        t1 = model->dPWTime[i - 1];
      } else {
        t1 = 0.0;
      }

      if (i == index || 0 == index) {
        /* Last part */
        t2 = CalibInst->dExeTime;
      } else {
        t2 = model->dPWTime[i];
      }

      cum_var += model->dSigma[i] * model->dSigma[i] * (t2 - t1);
    }

    ex_zeta1 = cum_var * exp(2.0 / model->dTau * model->dTStar);

    if (model->iOne2F == 2) {
      cum_var2 = 0.0;
      cum_var12 = 0.0;

      for (i = 0; i < index + 1; i++) {
        if (i > 0) {
          t1 = model->dPWTime[i - 1];
        } else {
          t1 = 0.0;
        }

        if (i == index || 0 == index) {
          /* Last part */
          t2 = CalibInst->dExeTime;
        } else {
          t2 = model->dPWTime[i];
        }

        cum_var2 += model->dLGMAlpha[i] * model->dLGMAlpha[i] *
                    model->dSigma[i] * model->dSigma[i] * (t2 - t1);
        cum_var12 += model->dLGMAlpha[i] * model->dLGMRho[i] *
                     model->dSigma[i] * model->dSigma[i] * (t2 - t1);
      }

      ex_zeta2 = cum_var2 * exp(2.0 / model->dTau2 * model->dTStar);
      ex_zeta12 =
          cum_var12 * exp((model->dLambdaX + model->dLambdaX2) * model->dTStar);
    }

    strike[1] = swp_cash;

    if (model->iOne2F == 1) {
      price[1] = lgmsopval1F(CalibInst->iNbCoupon,
                             &(CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn]),
                             &(CalibCpnSchedule->dCpnCvg[CalibInst->iStartCpn]),
                             HestonParam->dCpnG1, ex_zeta1, HestonParam->dExG1,
                             strike[1]);
    } else {
      price[1] = lgmsopval2F(CalibInst->iNbCoupon,
                             &(CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn]),
                             &(CalibCpnSchedule->dCpnCvg[CalibInst->iStartCpn]),
                             HestonParam->dCpnG1, HestonParam->dCpnG2, ex_zeta1,
                             ex_zeta2, ex_zeta12, HestonParam->dExG1,
                             HestonParam->dExG2, strike[1]);
    }

    err = srt_f_optimpvol(price[1], swp_cash, strike[1], CalibInst->dExeTime,
                          CalibInst->dLevel, SRT_PUT, SRT_NORMAL, &(vol[1]));

    if (err) {
      goto FREE_RETURN;
    }

    std = vol[1] * sqrt(CalibInst->dExeTime);
    HestonParam->dATMNormalStd = std;

    /* update the limit strikes */
    HestonParam->dMinStrike =
        CalibInst->dFwdCash + CalibInst->dSpread - HestonParam->dMinStd * std;
    HestonParam->dMaxStrike =
        CalibInst->dFwdCash + CalibInst->dSpread + HestonParam->dMaxStd * std;

    /* Now calculate the shift */
    if (CalibInst->iNbCoupon == 2) {
      /* special case of the FRA */
      HestonParam->dShift = CalibInst->dSumDf / CalibInst->dLevel;
      HestonParam->dCoefVol =
          (exp(-(CalibCpnSchedule->dCpnTime[CalibInst->iStartCpn] -
                 model->dTStar) /
               model->dTau) -
           exp(-(CalibCpnSchedule->dCpnTime[CalibInst->iEndCpn] -
                 model->dTStar) /
               model->dTau)) *
          model->dTau;

      if (model->iOne2F == 2) {
        HestonParam->dCoefVol_2F =
            (exp(-(CalibCpnSchedule->dCpnTime[CalibInst->iStartCpn] -
                   model->dTStar) /
                 model->dTau2) -
             exp(-(CalibCpnSchedule->dCpnTime[CalibInst->iEndCpn] -
                   model->dTStar) /
                 model->dTau2)) *
            model->dTau2;
      }

      /* adjustment for the spread */
      HestonParam->dShift -= CalibInst->dSpread;

      HestonParam->iNegVol = 0;

      if (CalibInst->dFwdCash + HestonParam->dShift < 0.0) {
        HestonParam->iNegVol = 1;
      }
    } else {
      strike[0] = strike[1] - 0.5 * std;
      strike[2] = strike[1] + 0.5 * std;

      if (model->iOne2F == 1) {
        price[0] = lgmsopval1F(
            CalibInst->iNbCoupon,
            &(CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn]),
            &(CalibCpnSchedule->dCpnCvg[CalibInst->iStartCpn]),
            HestonParam->dCpnG1, ex_zeta1, HestonParam->dExG1, strike[0]);
      } else {
        price[0] = lgmsopval2F(
            CalibInst->iNbCoupon,
            &(CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn]),
            &(CalibCpnSchedule->dCpnCvg[CalibInst->iStartCpn]),
            HestonParam->dCpnG1, HestonParam->dCpnG2, ex_zeta1, ex_zeta2,
            ex_zeta12, HestonParam->dExG1, HestonParam->dExG2, strike[0]);
      }

      err = srt_f_optimpvol(price[0], swp_cash, strike[0], CalibInst->dExeTime,
                            CalibInst->dLevel, SRT_PUT, SRT_NORMAL, &(vol[0]));

      if (err) {
        goto FREE_RETURN;
      }

      if (model->iOne2F == 1) {
        price[2] = lgmsopval1F(
            CalibInst->iNbCoupon,
            &(CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn]),
            &(CalibCpnSchedule->dCpnCvg[CalibInst->iStartCpn]),
            HestonParam->dCpnG1, ex_zeta1, HestonParam->dExG1, strike[2]);
      } else {
        price[2] = lgmsopval2F(
            CalibInst->iNbCoupon,
            &(CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn]),
            &(CalibCpnSchedule->dCpnCvg[CalibInst->iStartCpn]),
            HestonParam->dCpnG1, HestonParam->dCpnG2, ex_zeta1, ex_zeta2,
            ex_zeta12, HestonParam->dExG1, HestonParam->dExG2, strike[2]);
      }

      err = srt_f_optimpvol(price[2], swp_cash, strike[2], CalibInst->dExeTime,
                            CalibInst->dLevel, SRT_PUT, SRT_NORMAL, &(vol[2]));

      if (err) {
        goto FREE_RETURN;
      }

      if (vol[2] - vol[1] < 1.0E-08 || vol[1] - vol[0] < 1.0E-08) {
        HestonParam->dShift = 100;
      } else {
        /* Linear regression */
        sumx = 0.0;
        sumx2 = 0.0;
        sumy = 0.0;
        sumxy = 0.0;

        for (i = 0; i < 3; i++) {
          sumx += strike[i];
          sumx2 += strike[i] * strike[i];
          sumy += vol[i];
          sumxy += strike[i] * vol[i];
        }

        sumx /= 3.0;
        sumx2 /= 3.0;
        sumy /= 3.0;
        sumxy /= 3.0;

        coef_lin = (sumxy - sumx * sumy) / (sumx2 - sumx * sumx);
        const_lin = sumy - coef_lin * sumx;

        HestonParam->dShift =
            (const_lin - coef_lin * swp_cash) / (2.0 * coef_lin);
      }

      /* adjustment for the spread */
      HestonParam->dShift -= CalibInst->dSpread;

      err = srt_f_optimpvol(price[1], fabs(swp_cash + HestonParam->dShift),
                            fabs(strike[1] + HestonParam->dShift),
                            CalibInst->dExeTime, CalibInst->dLevel, SRT_PUT,
                            SRT_LOGNORMAL, &HestonParam->dCoefVol);

      if (err) {
        goto FREE_RETURN;
      }

      if (model->iOne2F == 1) {
        HestonParam->dCoefVol =
            HestonParam->dCoefVol / sqrt(cum_var / CalibInst->dExeTime);
      } else {
        /* we need to find the ratio between the two brownian contributions */
        price1 = lgmsopval1F(CalibInst->iNbCoupon,
                             &(CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn]),
                             &(CalibCpnSchedule->dCpnCvg[CalibInst->iStartCpn]),
                             HestonParam->dCpnG1, ex_zeta1, HestonParam->dExG1,
                             strike[1]);
        price2 = lgmsopval1F(CalibInst->iNbCoupon,
                             &(CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn]),
                             &(CalibCpnSchedule->dCpnCvg[CalibInst->iStartCpn]),
                             HestonParam->dCpnG2, ex_zeta2, HestonParam->dExG2,
                             strike[1]);

        err = srt_f_optimpvol(price1, fabs(swp_cash + HestonParam->dShift),
                              fabs(strike[1] + HestonParam->dShift),
                              CalibInst->dExeTime, CalibInst->dLevel, SRT_PUT,
                              SRT_LOGNORMAL, &vol1);

        if (err)
          goto FREE_RETURN;

        err = srt_f_optimpvol(price2, fabs(swp_cash + HestonParam->dShift),
                              fabs(strike[1] + HestonParam->dShift),
                              CalibInst->dExeTime, CalibInst->dLevel, SRT_PUT,
                              SRT_LOGNORMAL, &vol2);

        if (err)
          goto FREE_RETURN;

        ratio2F = vol2 / vol1 * sqrt(cum_var / cum_var2);

        HestonParam->dCoefVol = HestonParam->dCoefVol /
                                sqrt((cum_var + cum_var2 * ratio2F * ratio2F +
                                      2.0 * cum_var12 * ratio2F) /
                                     CalibInst->dExeTime);
        HestonParam->dCoefVol_2F = HestonParam->dCoefVol * ratio2F;
      }

      HestonParam->iNegVol = 0;

      if (swp_cash + HestonParam->dShift < 0.0) {
        HestonParam->iNegVol = 1;
      }
    }

    HestonParam->dNewCoefCMS = HestonParam->dCoefCMS * HestonParam->dCoefVol;
    HestonParam->dNewCoefCMS2 = HestonParam->dCoefCMS2 * HestonParam->dCoefVol;

    if (model->iOne2F == 2) {
      HestonParam->dNewCoefCMS_2F =
          HestonParam->dCoefCMS_2F * HestonParam->dCoefVol_2F;
      HestonParam->dNewCoefCMS2_2F =
          HestonParam->dCoefCMS2_2F * HestonParam->dCoefVol_2F;

      HestonParam->dNewCoefCMS_cross_2F =
          HestonParam->dCoefCMS_2F * HestonParam->dCoefVol +
          HestonParam->dCoefCMS * HestonParam->dCoefVol_2F;

      HestonParam->dNewCoefCMS2_cross_2F =
          HestonParam->dCoefCMS2_2F * HestonParam->dCoefVol +
          HestonParam->dCoefCMS2 * HestonParam->dCoefVol_2F;
    }
  }

FREE_RETURN:

  return err;
}

Err Calculate_HestonEquivalent_LGMSV_FromLGM_struct(
    /*	Market information */
    /* product information */
    CALIBINSTRUMENTDLM CalibInst, CALIBCPNSCHEDULEDLM CalibCpnSchedule,

    /* model information */
    LGMSV_MODEL model,

    /* output */
    LGMSV_HESTONINST HestonParam) {
  int i;
  double expTi1, dCoefMeanRev, dCoefCMS, temp;
  Err err = NULL;

  dCoefMeanRev = 0.0;

  for (i = 1; i < CalibInst->iNbCoupon; i++) {
    expTi1 = exp(-(CalibCpnSchedule->dCpnTime[CalibInst->iStartCpn + i] -
                   model->dTStar) /
                 model->dTau);
    dCoefMeanRev +=
        CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn + i] * (1.0 - expTi1);
  }

  dCoefMeanRev *= model->dTau;
  dCoefMeanRev /= CalibInst->dSumDf;

  if (CalibInst->iIsCMS || CalibInst->iIsLiborOption) {
    temp = (1.0 - exp(-(CalibInst->dPayTime - model->dTStar) / model->dTau)) *
           model->dTau;
    dCoefCMS = dCoefMeanRev - temp;
    dCoefMeanRev = temp;

    HestonParam->dCoefMeanRev = dCoefMeanRev;
    HestonParam->dCoefCMS = dCoefCMS;

    if (CalibInst->iIsLiborOption) {
      temp =
          (exp(-(CalibInst->dLiborEndTime - model->dTStar) / model->dTau) -
           exp(-(CalibInst->dLiborStartTime - model->dTStar) / model->dTau)) *
          model->dTau;
      dCoefCMS -= temp;
      dCoefMeanRev += temp;

      HestonParam->dCoefMeanRev2 = dCoefMeanRev;
      HestonParam->dCoefCMS2 = dCoefCMS;

      /* then switch */
      temp = HestonParam->dCoefMeanRev2;
      HestonParam->dCoefMeanRev2 = HestonParam->dCoefMeanRev;
      HestonParam->dCoefMeanRev = temp;

      temp = HestonParam->dCoefCMS2;
      HestonParam->dCoefCMS2 = HestonParam->dCoefCMS;
      HestonParam->dCoefCMS = temp;
    }
  } else {
    HestonParam->dCoefCMS = 0.0;
    HestonParam->dCoefCMS2 = 0.0;
    HestonParam->dCoefMeanRev = dCoefMeanRev;
  }

  if (model->iOne2F == 2) {
    /* we need to add the adjustment due to the second factor */
    dCoefMeanRev = 0.0;

    for (i = 1; i < CalibInst->iNbCoupon; i++) {
      expTi1 = exp(-(CalibCpnSchedule->dCpnTime[CalibInst->iStartCpn + i] -
                     model->dTStar) /
                   model->dTau2);
      dCoefMeanRev +=
          CalibCpnSchedule->dCpnDf[CalibInst->iStartCpn + i] * (1.0 - expTi1);
    }

    dCoefMeanRev *= model->dTau2;
    dCoefMeanRev /= CalibInst->dSumDf;

    if (CalibInst->iIsCMS || CalibInst->iIsLiborOption) {
      temp =
          (1.0 - exp(-(CalibInst->dPayTime - model->dTStar) / model->dTau2)) *
          model->dTau2;
      dCoefCMS = dCoefMeanRev - temp;
      dCoefMeanRev = temp;

      HestonParam->dCoefMeanRev_2F = dCoefMeanRev;
      HestonParam->dCoefCMS_2F = dCoefCMS;

      if (CalibInst->iIsLiborOption) {
        temp =
            (exp(-(CalibInst->dLiborEndTime - model->dTStar) / model->dTau2) -
             exp(-(CalibInst->dLiborStartTime - model->dTStar) /
                 model->dTau2)) *
            model->dTau2;
        dCoefCMS -= temp;
        dCoefMeanRev += temp;

        HestonParam->dCoefMeanRev2_2F = dCoefMeanRev;
        HestonParam->dCoefCMS2_2F = dCoefCMS;

        /* then switch */
        temp = HestonParam->dCoefMeanRev2_2F;
        HestonParam->dCoefMeanRev2_2F = HestonParam->dCoefMeanRev_2F;
        HestonParam->dCoefMeanRev_2F = temp;

        temp = HestonParam->dCoefCMS2_2F;
        HestonParam->dCoefCMS2_2F = HestonParam->dCoefCMS_2F;
        HestonParam->dCoefCMS_2F = temp;
      }
    } else {
      HestonParam->dCoefCMS_2F = 0.0;
      HestonParam->dCoefCMS2_2F = 0.0;
      HestonParam->dCoefMeanRev_2F = dCoefMeanRev;
    }
  }

  /* Fill the precalculations for LGM */

  HestonParam->dCpnG1[0] =
      (1.0 -
       exp(-CalibCpnSchedule->dCpnTime[CalibInst->iStartCpn] / model->dTau)) *
      model->dTau;

  for (i = 1; i < CalibInst->iNbCoupon; i++) {
    HestonParam->dCpnG1[i] =
        (1.0 - exp(-CalibCpnSchedule->dCpnTime[CalibInst->iStartCpn + i] /
                   model->dTau)) *
        model->dTau;
  }

  HestonParam->dExG1 =
      (1.0 - exp(-CalibInst->dExeTime / model->dTau)) * model->dTau;

  if (model->iOne2F == 2) {
    HestonParam->dCpnG2[0] =
        (1.0 - exp(-CalibCpnSchedule->dCpnTime[CalibInst->iStartCpn] /
                   model->dTau2)) *
        model->dTau2;

    for (i = 1; i < CalibInst->iNbCoupon; i++) {
      HestonParam->dCpnG2[i] =
          (1.0 - exp(-CalibCpnSchedule->dCpnTime[CalibInst->iStartCpn + i] /
                     model->dTau2)) *
          model->dTau2;
    }

    HestonParam->dExG2 =
        (1.0 - exp(-CalibInst->dExeTime / model->dTau2)) * model->dTau2;
  }

  return err;
}

void LGMSVClosedFormApprox_struct(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    LGMSV_MODEL model,

    /* Product description */
    CALIBINSTRUMENTDLM CalibInst, LGMSV_HESTONINST HestonParam,

    /* Parameter of grids */
    LGMSV_PRICINGCONST NumerParams, LGMSV_NUMERINST InstNumer,

    /* Outputs */
    double *Price) {
  double temp_vol, temp_var, XMax;
  double temp_vol2, temp_vol_tot, correl_tot;
  double t1, t2;
  int i, i2;
  int iIntegMethod;
  int do_fwd;

  if (fabs(CalibInst->dStrike[InstNumer->iNbStrike - 1]) < 1.0E-08) {
    do_fwd = 1;
  } else {
    do_fwd = 0;
  }

  if (CalibInst->dExeTime > 1.0E-08) {
    /* Precalculations */

    Update_NumerInst(CalibInst, HestonParam, InstNumer);

    if (HestonParam->iNegVol) {
      NumerParams->NumerParams->dIntegParam =
          -1.0 - NumerParams->NumerParams->dIntegParam;
    }

    t1 = 0.0;
    temp_var = 0.0;

    if (model->iOne2F == 1) {
      /* First part */
      for (i = 0; i <= InstNumer->iSwitchIndex; i++) {
        /* First the time discretisation */
        if (i < InstNumer->iSwitchIndex) {
          t2 = model->dPWTime[i];
        } else {
          t2 = InstNumer->dNewSwitchTime;
        }

        /* Precalculation on the option i*/
        temp_vol = HestonParam->dCoefVol * model->dSigma[i];
        NumerParams->dt[i] = (t2 - t1);
        NumerParams->CoefIntT[i] = model->dLvlEps[i];
        NumerParams->Coef1ReT[i] = -0.5 * model->dAlpha[i] * model->dAlpha[i];
        NumerParams->Coef2ReT[i] =
            model->dLambdaEps[i] + HestonParam->dCoefMeanRev *
                                       model->dAlpha[i] * model->dRho[i] *
                                       model->dSigma[i];
        NumerParams->Coef2ImT[i] =
            -model->dAlpha[i] * model->dRho[i] * temp_vol;
        temp_vol *= temp_vol;
        temp_var += temp_vol * NumerParams->dt[i];
        NumerParams->Coef3ReT[i] = 0.5 * temp_vol;
        NumerParams->Coef3ImT[i] = 0.5 * temp_vol - HestonParam->dNewCoefCMS *
                                                        model->dSigma[i] *
                                                        model->dSigma[i];

        t1 = t2;
      }

      /* Second part */
      for (i = i; i <= InstNumer->iEndIndex; i++) {
        /* First the time discretisation */
        if (i < InstNumer->iEndIndex) {
          t2 = model->dPWTime[i - InstNumer->iDoSwitch];
        } else {
          t2 = CalibInst->dExeTime;
        }

        i2 = i - InstNumer->iDoSwitch;

        /* Precalculation on the option i*/
        temp_vol = HestonParam->dCoefVol * model->dSigma[i2];
        NumerParams->dt[i] = (t2 - t1);
        NumerParams->CoefIntT[i] = model->dLvlEps[i2];
        NumerParams->Coef1ReT[i] = -0.5 * model->dAlpha[i2] * model->dAlpha[i2];
        NumerParams->Coef2ReT[i] =
            model->dLambdaEps[i2] + HestonParam->dCoefMeanRev2 *
                                        model->dAlpha[i2] * model->dRho[i2] *
                                        model->dSigma[i2];
        NumerParams->Coef2ImT[i] =
            -model->dAlpha[i2] * model->dRho[i2] * temp_vol;
        temp_vol *= temp_vol;
        temp_var += temp_vol * NumerParams->dt[i];
        NumerParams->Coef3ReT[i] = 0.5 * temp_vol;
        NumerParams->Coef3ImT[i] = 0.5 * temp_vol - HestonParam->dNewCoefCMS2 *
                                                        model->dSigma[i2] *
                                                        model->dSigma[i2];

        t1 = t2;
      }

      temp_var = sqrt(temp_var);
    } else {
      /* First part */
      for (i = 0; i <= InstNumer->iSwitchIndex; i++) {
        /* First the time discretisation */
        if (i < InstNumer->iSwitchIndex) {
          t2 = model->dPWTime[i];
        } else {
          t2 = InstNumer->dNewSwitchTime;
        }

        /* Precalculation on the option i*/
        temp_vol = HestonParam->dCoefVol * model->dSigma[i];
        temp_vol2 =
            HestonParam->dCoefVol_2F * model->dSigma[i] * model->dLGMAlpha[i];
        temp_vol_tot = sqrt(temp_vol * temp_vol + temp_vol2 * temp_vol2 +
                            2.0 * model->dLGMRho[i] * temp_vol * temp_vol2);
        correl_tot = (model->dRho[i] * temp_vol + model->dRho2[i] * temp_vol2) /
                     temp_vol_tot;

        NumerParams->dt[i] = (t2 - t1);
        NumerParams->CoefIntT[i] = model->dLvlEps[i];
        NumerParams->Coef1ReT[i] = -0.5 * model->dAlpha[i] * model->dAlpha[i];
        NumerParams->Coef2ReT[i] =
            model->dLambdaEps[i] +
            model->dAlpha[i] * model->dSigma[i] *
                (HestonParam->dCoefMeanRev * model->dRho[i] +
                 HestonParam->dCoefMeanRev_2F * model->dRho2[i] *
                     model->dLGMAlpha[i]);
        NumerParams->Coef2ImT[i] =
            -model->dAlpha[i] * correl_tot * temp_vol_tot;
        temp_vol_tot *= temp_vol_tot;
        temp_var += temp_vol_tot * NumerParams->dt[i];
        NumerParams->Coef3ReT[i] = 0.5 * temp_vol_tot;
        NumerParams->Coef3ImT[i] =
            0.5 * temp_vol_tot -
            model->dSigma[i] * model->dSigma[i] *
                (HestonParam->dNewCoefCMS +
                 HestonParam->dNewCoefCMS_2F * model->dLGMAlpha[i] *
                     model->dLGMAlpha[i] +
                 model->dLGMRho[i] * model->dLGMAlpha[i] *
                     HestonParam->dNewCoefCMS_cross_2F);

        t1 = t2;
      }

      /* Second part */
      for (i = i; i <= InstNumer->iEndIndex; i++) {
        /* First the time discretisation */
        if (i < InstNumer->iEndIndex) {
          t2 = model->dPWTime[i - InstNumer->iDoSwitch];
        } else {
          t2 = CalibInst->dExeTime;
        }

        i2 = i - InstNumer->iDoSwitch;

        /* Precalculation on the option i*/
        temp_vol = HestonParam->dCoefVol * model->dSigma[i2];
        temp_vol2 =
            HestonParam->dCoefVol_2F * model->dSigma[i2] * model->dLGMAlpha[i2];
        temp_vol_tot = sqrt(temp_vol * temp_vol + temp_vol2 * temp_vol2 +
                            2.0 * model->dLGMRho[i2] * temp_vol * temp_vol2);
        correl_tot =
            (model->dRho[i2] * temp_vol + model->dRho2[i2] * temp_vol2) /
            temp_vol_tot;

        NumerParams->dt[i] = (t2 - t1);
        NumerParams->CoefIntT[i] = model->dLvlEps[i2];
        NumerParams->Coef1ReT[i] = -0.5 * model->dAlpha[i2] * model->dAlpha[i2];
        NumerParams->Coef2ReT[i] =
            model->dLambdaEps[i2] +
            model->dAlpha[i2] * model->dSigma[i2] *
                (HestonParam->dCoefMeanRev2 * model->dRho[i2] +
                 HestonParam->dCoefMeanRev2_2F * model->dRho2[i2] *
                     model->dLGMAlpha[i2]);
        NumerParams->Coef2ImT[i] =
            -model->dAlpha[i2] * correl_tot * temp_vol_tot;
        temp_vol_tot *= temp_vol_tot;
        temp_var += temp_vol_tot * NumerParams->dt[i];
        NumerParams->Coef3ReT[i] = 0.5 * temp_vol_tot;
        NumerParams->Coef3ImT[i] =
            0.5 * temp_vol_tot -
            model->dSigma[i2] * model->dSigma[i2] *
                (HestonParam->dNewCoefCMS2 +
                 HestonParam->dNewCoefCMS2_2F * model->dLGMAlpha[i2] *
                     model->dLGMAlpha[i2] +
                 model->dLGMRho[i2] * model->dLGMAlpha[i2] *
                     HestonParam->dNewCoefCMS2_cross_2F);

        t1 = t2;
      }

      temp_var = sqrt(temp_var);
    }

    /* select the method */
    if (NumerParams->NumerParams->iIntegMethod == 4) {
      /* select the method */
      if (temp_var < NumerParams->NumerParams->dVolLimit) {
        /* special case for low volatility */
        iIntegMethod = 3;
      } else {
        /* classical case */
        iIntegMethod = 1;
      }
    } else {
      iIntegMethod = NumerParams->NumerParams->iIntegMethod;
    }

    /* construct the grid */

    if (iIntegMethod == 0) {
      /* simple unifrom grid */
      XMax = (NumerParams->NumerParams->dParam1 +
              NumerParams->NumerParams->dParam2) *
             HestonParam->dCoefVol *
             model->dSigma[InstNumer->iEndIndex - InstNumer->iDoSwitch] /
             (NumerParams->iNbX - 1);
      XMax = 1.0 / XMax;

      LGMSVCalculateGride(NumerParams->iNbX, XMax, NumerParams->X);
    }

    /* calculate the integrand */
    if (iIntegMethod <= 1) {
      LGMSVSolveODEGride(
          NumerParams->iNbX, NumerParams->X, InstNumer->dLogShiftFwd,
          InstNumer->iEndIndex, NumerParams->dt, NumerParams->CoefIntT,
          NumerParams->Coef1ReT, NumerParams->Coef1ReT, NumerParams->Coef2ReT,
          NumerParams->Coef2ImT, NumerParams->Coef3ReT, NumerParams->Coef3ImT,
          NumerParams->NumerParams->dIntegParam, NumerParams->IntRe,
          NumerParams->IntIm);
    } else if (iIntegMethod == 2) {
      LGMSVSolveODEGrideTV(
          NumerParams->iNbX, NumerParams->X, InstNumer->dLogShiftFwd,
          InstNumer->iEndIndex, NumerParams->dt, NumerParams->CoefIntT,
          NumerParams->Coef1ReT, NumerParams->Coef1ReT, NumerParams->Coef2ReT,
          NumerParams->Coef2ImT, NumerParams->Coef3ReT, NumerParams->Coef3ImT,
          NumerParams->NumerParams->dIntegParam, NumerParams->IntRe,
          NumerParams->IntIm);
    } else if (iIntegMethod == 3) {
      LGMSVSolveODEGride(
          NumerParams->iNbX, NumerParams->XX, InstNumer->dLogShiftFwd,
          InstNumer->iEndIndex, NumerParams->dt, NumerParams->CoefIntT,
          NumerParams->Coef1ReT, NumerParams->Coef1ReT, NumerParams->Coef2ReT,
          NumerParams->Coef2ImT, NumerParams->Coef3ReT, NumerParams->Coef3ImT,
          NumerParams->NumerParams->dIntegParam, NumerParams->IntRe,
          NumerParams->IntIm);
    }

    /* do the integral */
    if (iIntegMethod == 0) {
      LGMSVComputeIntegral(NumerParams->iNbX, NumerParams->X,
                           NumerParams->IntRe, NumerParams->IntIm,
                           InstNumer->iNbStrike - do_fwd,
                           InstNumer->dLogShiftStrike,
                           NumerParams->NumerParams->dIntegParam, Price);
    } else if (iIntegMethod == 1) {
      LGMSVComputeIntegralLaguerre(
          NumerParams->iNbX, NumerParams->X, NumerParams->W, 0.0,
          NumerParams->IntRe, NumerParams->IntIm, InstNumer->iNbStrike - do_fwd,
          InstNumer->dLogShiftStrike, NumerParams->NumerParams->dIntegParam,
          Price);
    } else if (iIntegMethod == 2) {
      LGMSVComputeIntegralHermite(
          NumerParams->iNbX, NumerParams->X, NumerParams->W, NumerParams->IntRe,
          NumerParams->IntIm, InstNumer->iNbStrike - do_fwd,
          InstNumer->dLogShiftStrike, NumerParams->NumerParams->dIntegParam,
          Price);

      for (i = 0; i < InstNumer->iNbStrike - do_fwd; i++) {
        if (InstNumer->dLogShiftStrike[i] < 0) {
          Price[i] +=
              CalibInst->dFwdCash + CalibInst->dSpread - CalibInst->dStrike[i];
        } else {
          Price[i] -=
              CalibInst->dFwdCash + CalibInst->dSpread - CalibInst->dStrike[i];
        }
      }
    } else if (iIntegMethod == 3) {
      LGMSVComputeIntegralLegendre(
          NumerParams->iNbX, NumerParams->XX, NumerParams->WW,
          NumerParams->IntRe, NumerParams->IntIm, InstNumer->iNbStrike - do_fwd,
          InstNumer->dLogShiftStrike, NumerParams->NumerParams->dIntegParam,
          Price);
    }

    if (HestonParam->iNegVol) {
      NumerParams->NumerParams->dIntegParam =
          -1.0 - NumerParams->NumerParams->dIntegParam;
    }
  } else {
    for (i = 0; i < InstNumer->iNbStrike - do_fwd; i++) {
      if (NumerParams->NumerParams->dIntegParam < 0.0) {
        /* Payer */
        Price[i] = DMAX(CalibInst->dFwdCash + CalibInst->dSpread -
                            CalibInst->dStrike[i],
                        0.0);
      } else {
        /* Receiver */
        Price[i] = DMAX(CalibInst->dStrike[i] - CalibInst->dFwdCash -
                            CalibInst->dSpread,
                        0.0);
      }
    }
  }

  if (!CalibInst->iIsCMS) {
    for (i = 0; i < InstNumer->iNbStrike - do_fwd; i++) {
      Price[i] *= CalibInst->dLevel;
    }
  }

  if (do_fwd) {
    if (CalibInst->dExeTime > 1.0E-08) {
      LGMSVSolveODE(0.0, -1.0, InstNumer->dLogShiftFwd, InstNumer->iEndIndex,
                    NumerParams->dt, NumerParams->CoefIntT,
                    NumerParams->Coef1ReT, NumerParams->Coef1ReT,
                    NumerParams->Coef2ReT, NumerParams->Coef2ImT,
                    NumerParams->Coef3ReT, NumerParams->Coef3ImT,
                    NumerParams->IntRe, NumerParams->IntIm);

      if (HestonParam->dCoefVol > 0.0) {
        Price[InstNumer->iNbStrike - 1] =
            exp(NumerParams->IntRe[0]) - HestonParam->dShift;
      } else {
        Price[InstNumer->iNbStrike - 1] =
            -exp(NumerParams->IntRe[0]) - HestonParam->dShift;
      }
    } else {
      Price[InstNumer->iNbStrike - 1] =
          CalibInst->dFwdCash + CalibInst->dSpread;
    }
  }
}

void LGMSVFwdClosedFormApprox_struct(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    LGMSV_MODEL model,

    /* Product description */
    CALIBINSTRUMENTDLM CalibInst, LGMSV_HESTONINST HestonParam,

    /* Parameter of grids */
    LGMSV_PRICINGCONST NumerParams, LGMSV_NUMERINST InstNumer,

    /* Outputs */
    double *Price) {
  double temp_vol, temp_vol2, temp_vol_tot, correl_tot;
  int i;
  double t1, t2;

  /* Precalculations */
  if (CalibInst->dExeTime > 1.0E-08) {

    Update_NumerInst(CalibInst, HestonParam, InstNumer);

    t1 = 0.0;

    if (model->iOne2F == 1) {
      for (i = 0; i <= InstNumer->iLiborIndex; i++) {
        /* First the time discretisation */
        if (i < InstNumer->iLiborIndex) {
          t2 = model->dPWTime[i];
        } else {
          t2 = CalibInst->dLiborFixTime;
        }

        /* Precalculation on the option i*/
        temp_vol = HestonParam->dCoefVol * model->dSigma[i];
        NumerParams->dt[i] = (t2 - t1);
        NumerParams->CoefIntT[i] = model->dLvlEps[i];
        NumerParams->Coef1ReT[i] = -0.5 * model->dAlpha[i] * model->dAlpha[i];
        NumerParams->Coef2ReT[i] =
            model->dLambdaEps[i] + HestonParam->dCoefMeanRev *
                                       model->dAlpha[i] * model->dRho[i] *
                                       model->dSigma[i];
        NumerParams->Coef2ImT[i] =
            -model->dAlpha[i] * model->dRho[i] * temp_vol;
        temp_vol *= temp_vol;
        NumerParams->Coef3ReT[i] = 0.5 * temp_vol;
        NumerParams->Coef3ImT[i] = 0.5 * temp_vol - HestonParam->dNewCoefCMS *
                                                        model->dSigma[i] *
                                                        model->dSigma[i];

        t1 = t2;
      }
    } else {
      for (i = 0; i <= InstNumer->iLiborIndex; i++) {
        /* First the time discretisation */
        if (i < InstNumer->iLiborIndex) {
          t2 = model->dPWTime[i];
        } else {
          t2 = CalibInst->dLiborFixTime;
        }

        /* Precalculation on the option i*/
        temp_vol = HestonParam->dCoefVol * model->dSigma[i];
        temp_vol2 =
            HestonParam->dCoefVol_2F * model->dSigma[i] * model->dLGMAlpha[i];
        temp_vol_tot = sqrt(temp_vol * temp_vol + temp_vol2 * temp_vol2 +
                            2.0 * model->dLGMRho[i] * temp_vol * temp_vol2);
        correl_tot = (model->dRho[i] * temp_vol + model->dRho2[i] * temp_vol2) /
                     temp_vol_tot;

        NumerParams->dt[i] = (t2 - t1);
        NumerParams->CoefIntT[i] = model->dLvlEps[i];
        NumerParams->Coef1ReT[i] = -0.5 * model->dAlpha[i] * model->dAlpha[i];
        NumerParams->Coef2ReT[i] =
            model->dLambdaEps[i] +
            model->dAlpha[i] * model->dSigma[i] *
                (HestonParam->dCoefMeanRev * model->dRho[i] +
                 HestonParam->dCoefMeanRev_2F * model->dAlpha[i] *
                     model->dRho2[i]);
        NumerParams->Coef2ImT[i] =
            -model->dAlpha[i] * correl_tot * temp_vol_tot;
        temp_vol_tot *= temp_vol_tot;
        NumerParams->Coef3ReT[i] = 0.5 * temp_vol_tot;
        NumerParams->Coef3ImT[i] =
            0.5 * temp_vol_tot -
            model->dSigma[i] * model->dSigma[i] *
                (HestonParam->dNewCoefCMS +
                 HestonParam->dNewCoefCMS_2F * model->dLGMAlpha[i] *
                     model->dLGMAlpha[i] +
                 model->dLGMRho[i] * model->dLGMAlpha[i] *
                     HestonParam->dNewCoefCMS_cross_2F);

        t1 = t2;
      }
    }

    LGMSVSolveODE(0.0, -1.0, InstNumer->dLogShiftFwd, InstNumer->iLiborIndex,
                  NumerParams->dt, NumerParams->CoefIntT, NumerParams->Coef1ReT,
                  NumerParams->Coef1ReT, NumerParams->Coef2ReT,
                  NumerParams->Coef2ImT, NumerParams->Coef3ReT,
                  NumerParams->Coef3ImT, NumerParams->IntRe,
                  NumerParams->IntIm);

    if (HestonParam->dCoefVol > 0.0) {
      (*Price) = exp(NumerParams->IntRe[0]) - HestonParam->dShift;
    } else {
      (*Price) = -exp(NumerParams->IntRe[0]) - HestonParam->dShift;
    }
  } else {
    *Price = CalibInst->dFwdCash + CalibInst->dSpread;
  }
}

Err ConstructSingleInstrumentSchedule(char *sYcName, long lToday,
                                      char *cInstFreq, char *cInstBasis,
                                      long lExeDate, long lStartDate,
                                      long lTheoEndDate,
                                      CALIBCPNSCHEDULEDLM CalibCpnSchedule,
                                      CALIBEXESCHEDULEDLM CalibExeSchedule) {
  SwapDP Swap;
  SrtDateList SwapDates;
  SrtBasisCode ibasis;
  int k;
  Err err = NULL;

  /* Initialisation */
  CalibCpnSchedule->lToday = lToday;

  /* Fill the Exe Schedule */
  CalibExeSchedule->iNExe = 1;
  CalibExeSchedule->lActEndDates[0] =
      bus_date_method(lTheoEndDate, MODIFIED_SUCCEEDING);
  CalibExeSchedule->lExeDates[0] = lExeDate;
  CalibExeSchedule->dExeTimes[0] = (lExeDate - lToday) * YEARS_IN_DAY;

  /* Fill the Cpn Schedule */
  err = interp_basis(cInstBasis, &ibasis);
  if (err)
    goto FREE_RETURN;

  err =
      swp_f_initSwapDP(lStartDate, lTheoEndDate, cInstFreq, cInstBasis, &Swap);

  if (err)
    goto FREE_RETURN;

  SwapDates = SwapDP_to_DateList(&Swap, MODIFIED_SUCCEEDING);

  CalibCpnSchedule->iNCpn = SwapDates.len;

  for (k = 0; k < CalibCpnSchedule->iNCpn; k++) {
    CalibCpnSchedule->lCpnDate[k] = SwapDates.date[k];
    CalibCpnSchedule->lCpnTheoDate[k] = SwapDates.date[k];
    CalibCpnSchedule->dCpnTime[k] =
        (CalibCpnSchedule->lCpnDate[k] - lToday) * YEARS_IN_DAY;
    CalibCpnSchedule->dCpnDf[k] =
        swp_f_df(lToday, CalibCpnSchedule->lCpnDate[k], sYcName);

    if (k > 0) {
      CalibCpnSchedule->dCpnCvg[k] =
          coverage(CalibCpnSchedule->lCpnDate[k - 1],
                   CalibCpnSchedule->lCpnDate[k], ibasis);
    } else {
      CalibCpnSchedule->dCpnCvg[k] = 0.0;
    }
  }

FREE_RETURN:

  swp_f_free_in_DateList(SwapDates);

  return err;
}

Err LGMSVOptionApprox_new(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lTheoEndDate, long lPayDate, int IsCMS,
    long lRStartDate, long lREndDate, int iIsRatioNum, int iNbStrike,
    double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice) {
  Err err = NULL;

  CalibCpnScheduleDLM *CpnSchedule = NULL;
  CalibExeScheduleDLM *ExeSchedule = NULL;
  CalibInstrumentDLM *Inst = NULL;
  AllCalibInstrumentsDLM *AllInst = NULL;

  SrtCurvePtr yc_ptr;
  long today;

  LGMSV_model *model = NULL;
  LGMSV_HestonInst *HestonParam = NULL;
  LGMSV_PricingConst *PricingConst = NULL;
  LGMSV_NumerInst *InstNumer = NULL;

  double dExeTime;
  long lEndDate;
  double dDfPay;

  char Diag_tenor[20] = "DIAG";
  char *tenor = &(Diag_tenor[0]);

  int i;
  int iIntegParam_needs_reset = 0;

  /* All allocations */
  CpnSchedule = calloc(1, sizeof(CalibCpnScheduleDLM));
  ExeSchedule = calloc(1, sizeof(CalibExeScheduleDLM));
  Inst = calloc(1, sizeof(CalibInstrumentDLM));
  AllInst = calloc(1, sizeof(AllCalibInstrumentsDLM));

  model = calloc(1, sizeof(LGMSV_model));
  HestonParam = calloc(1, sizeof(LGMSV_HestonInst));
  PricingConst = calloc(1, sizeof(LGMSV_PricingConst));
  InstNumer = calloc(1, sizeof(LGMSV_NumerInst));

  if (!CpnSchedule || !ExeSchedule || !Inst || !AllInst || !model ||
      !HestonParam || !PricingConst || !InstNumer) {
    err = "Memory allocation faillure in LGMSVOptionApprox_new";
    goto FREE_RETURN;
  }

  /* Get the curve informations */
  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found in LGMSVOptionApprox_new";
    goto FREE_RETURN;
  }

  today = get_today_from_curve(yc_ptr);

  /* Get the model */
  err = Get_LGMSV_model(und_name, model);

  if (err)
    goto FREE_RETURN;

  /* Get the instrument informations */
  AllocateCalibExeSchedule(&lExDate, &dExeTime, &lStartDate, &lTheoEndDate,
                           &lEndDate, dStrike, NULL, NULL, NULL, ExeSchedule);

  /* Old Code */
  /*
  err = Construct_CalibSchedule(	yc_name  ,
                                                                  today  ,
                                                                  swaption_freq
  , swaption_basis  , 1  , &lExDate  , NULL  , &tenor  , lTheoEndDate  , dStrike
  , NULL  , NULL  , NULL  , CpnSchedule  , ExeSchedule);
  */

  err = ConstructSingleInstrumentSchedule(
      yc_name, today, swaption_freq, swaption_basis, lExDate, lStartDate,
      lTheoEndDate, CpnSchedule, ExeSchedule);

  if (err)
    goto FREE_RETURN;

  err = AllocateCalibInst(iNbStrike, Inst);

  if (err)
    goto FREE_RETURN;

  for (i = 0; i < iNbStrike; i++) {
    Inst->dStrike[i] = dStrike[i];
  }

  /* Fill the calibration instrument */

  AllInst->iNbInst = 1;
  AllInst->sCalibInst = Inst;

  err = Calculate_CalibInst(today, yc_name, NULL, NULL, NULL, swaption_freq,
                            swaption_basis, ref_rate_name, 0, NULL, CpnSchedule,
                            ExeSchedule, AllInst, NULL);

  if (err)
    goto FREE_RETURN;

  /* Set the CMS */
  Inst->iIsCMS = IsCMS;
  Inst->lPayDate = lPayDate;
  Inst->dPayTime = (Inst->lPayDate - today) * YEARS_IN_DAY;

  /* All the precalculation and memory allocation for pricing */
  err = Initialise_PricingConst(model, NumerParams, PricingConst);

  if (err)
    goto FREE_RETURN;

  if (pay_rec == -1) {
    NumerParams->dIntegParam = -1 - NumerParams->dIntegParam;
    iIntegParam_needs_reset = 1;
  }

  /* Check if time value */
  if (dExeTime < 1.0E-06) {
    if (pay_rec == 1) {
      for (i = 0; i < iNbStrike; i++) {
        pSwaptionPrice[i] = max(Inst->dFwdCash + Inst->dSpread - dStrike[i], 0);
      }
    } else {
      for (i = 0; i < iNbStrike; i++) {
        pSwaptionPrice[i] = max(dStrike[i] - Inst->dFwdCash - Inst->dSpread, 0);
      }
    }

    if (IsCMS) {
      for (i = 0; i < iNbStrike; i++) {
        pSwaptionPrice[i] *= Inst->dLevel;
      }
    } else {
      dDfPay = swp_f_df(today, lPayDate, yc_name);

      for (i = 0; i < iNbStrike; i++) {
        pSwaptionPrice[i] *= dDfPay;
      }
    }

    goto FREE_RETURN;
  }

  err = Calculate_HestonEquivalent_LGMSV_FromLGM_struct(Inst, CpnSchedule,
                                                        model, HestonParam);

  if (err)
    goto FREE_RETURN;

  err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(Inst, CpnSchedule, model,
                                                   HestonParam);

  if (err)
    goto FREE_RETURN;

  err = Initialise_NumerInst(model, Inst, HestonParam, InstNumer);

  if (err)
    goto FREE_RETURN;

  Update_NumerInst(Inst, HestonParam, InstNumer);

  /* pricing */
  LGMSVClosedFormApprox_struct(model, Inst, HestonParam, PricingConst,
                               InstNumer, pSwaptionPrice);

  if (Inst->iIsCMS) {
    dDfPay = swp_f_df(today, lPayDate, yc_name);

    for (i = 0; i < iNbStrike; i++) {
      pSwaptionPrice[i] *= dDfPay;
    }
  }

FREE_RETURN:

  if (iIntegParam_needs_reset)
    NumerParams->dIntegParam = -1 - NumerParams->dIntegParam;

  if (CpnSchedule)
    free(CpnSchedule);
  if (ExeSchedule)
    free(ExeSchedule);

  if (Inst) {
    FreeCalibInst(Inst);
    free(Inst);
  }

  if (AllInst)
    free(AllInst);

  if (model) {
    free_LGMSV_model(model);
    free(model);
  }

  if (HestonParam)
    free(HestonParam);

  if (PricingConst) {
    Free_PricingConst(PricingConst);
    free(PricingConst);
  }

  if (InstNumer) {
    Free_NumerInst(InstNumer);
    free(InstNumer);
  }

  return err;
}

Err LGMSVOptionApproxFromModel(
    LGMSV_model *sModel,

    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lTheoEndDate, long lPayDate, int IsCMS,
    long lRStartDate, long lREndDate, int iIsRatioNum, int iNbStrike,
    double *dStrike, int pay_rec, /*	pay:1 rec:-1 */

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *pSwaptionPrice) {
  Err err = NULL;

  CalibCpnScheduleDLM *CpnSchedule = NULL;
  CalibExeScheduleDLM *ExeSchedule = NULL;
  CalibInstrumentDLM *Inst = NULL;
  AllCalibInstrumentsDLM *AllInst = NULL;

  SrtCurvePtr yc_ptr;
  long today;

  LGMSV_HestonInst *HestonParam = NULL;
  LGMSV_PricingConst *PricingConst = NULL;
  LGMSV_NumerInst *InstNumer = NULL;

  double dExeTime;
  long lEndDate;
  double dDfPay;

  char Diag_tenor[20] = "DIAG";
  char *tenor = &(Diag_tenor[0]);

  int i;
  int iIntegParam_needs_reset = 0;

  /* All allocations */
  CpnSchedule = calloc(1, sizeof(CalibCpnScheduleDLM));
  ExeSchedule = calloc(1, sizeof(CalibExeScheduleDLM));
  Inst = calloc(1, sizeof(CalibInstrumentDLM));
  AllInst = calloc(1, sizeof(AllCalibInstrumentsDLM));

  HestonParam = calloc(1, sizeof(LGMSV_HestonInst));
  PricingConst = calloc(1, sizeof(LGMSV_PricingConst));
  InstNumer = calloc(1, sizeof(LGMSV_NumerInst));

  if (!CpnSchedule || !ExeSchedule || !Inst || !AllInst || !sModel ||
      !HestonParam || !PricingConst || !InstNumer) {
    err = "Memory allocation faillure in LGMSVOptionApprox_new";
    goto FREE_RETURN;
  }

  /* Get the curve informations */
  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found in LGMSVOptionApprox_new";
    goto FREE_RETURN;
  }

  today = get_today_from_curve(yc_ptr);

  /* Get the instrument informations */
  AllocateCalibExeSchedule(&lExDate, &dExeTime, &lStartDate, &lTheoEndDate,
                           &lEndDate, dStrike, NULL, NULL, NULL, ExeSchedule);

  /* Old Code */
  /*
  err = Construct_CalibSchedule(	yc_name  ,
                                                                  today  ,
                                                                  swaption_freq
  , swaption_basis  , 1  , &lExDate  , NULL  , &tenor  , lTheoEndDate  , dStrike
  , NULL  , NULL  , NULL  , CpnSchedule  , ExeSchedule);
  */

  err = ConstructSingleInstrumentSchedule(
      yc_name, today, swaption_freq, swaption_basis, lExDate, lStartDate,
      lTheoEndDate, CpnSchedule, ExeSchedule);

  if (err)
    goto FREE_RETURN;

  err = AllocateCalibInst(iNbStrike, Inst);

  if (err)
    goto FREE_RETURN;

  for (i = 0; i < iNbStrike; i++) {
    Inst->dStrike[i] = dStrike[i];
  }

  /* Fill the calibration instrument */

  AllInst->iNbInst = 1;
  AllInst->sCalibInst = Inst;

  err = Calculate_CalibInst(today, yc_name, NULL, NULL, NULL, swaption_freq,
                            swaption_basis, ref_rate_name, 0, NULL, CpnSchedule,
                            ExeSchedule, AllInst, NULL);

  if (err)
    goto FREE_RETURN;

  /* Set the CMS */
  Inst->iIsCMS = IsCMS;
  Inst->lPayDate = lPayDate;
  Inst->dPayTime = (Inst->lPayDate - today) * YEARS_IN_DAY;

  /* All the precalculation and memory allocation for pricing */
  err = Initialise_PricingConst(sModel, NumerParams, PricingConst);

  if (err)
    goto FREE_RETURN;

  if (pay_rec == -1) {
    NumerParams->dIntegParam = -1 - NumerParams->dIntegParam;
    iIntegParam_needs_reset = 1;
  }

  /* Check if time value */
  if (dExeTime < 1.0E-06) {
    if (pay_rec == 1) {
      for (i = 0; i < iNbStrike; i++) {
        pSwaptionPrice[i] = max(Inst->dFwdCash + Inst->dSpread - dStrike[i], 0);
      }
    } else {
      for (i = 0; i < iNbStrike; i++) {
        pSwaptionPrice[i] = max(dStrike[i] - Inst->dFwdCash - Inst->dSpread, 0);
      }
    }

    if (IsCMS) {
      for (i = 0; i < iNbStrike; i++) {
        pSwaptionPrice[i] *= Inst->dLevel;
      }
    } else {
      dDfPay = swp_f_df(today, lPayDate, yc_name);

      for (i = 0; i < iNbStrike; i++) {
        pSwaptionPrice[i] *= dDfPay;
      }
    }

    goto FREE_RETURN;
  }

  err = Calculate_HestonEquivalent_LGMSV_FromLGM_struct(Inst, CpnSchedule,
                                                        sModel, HestonParam);

  if (err)
    goto FREE_RETURN;

  err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(Inst, CpnSchedule, sModel,
                                                   HestonParam);

  if (err)
    goto FREE_RETURN;

  err = Initialise_NumerInst(sModel, Inst, HestonParam, InstNumer);

  if (err)
    goto FREE_RETURN;

  Update_NumerInst(Inst, HestonParam, InstNumer);

  /* pricing */
  LGMSVClosedFormApprox_struct(sModel, Inst, HestonParam, PricingConst,
                               InstNumer, pSwaptionPrice);

  if (Inst->iIsCMS) {
    dDfPay = swp_f_df(today, lPayDate, yc_name);

    for (i = 0; i < iNbStrike; i++) {
      pSwaptionPrice[i] *= dDfPay;
    }
  }

FREE_RETURN:

  if (iIntegParam_needs_reset)
    NumerParams->dIntegParam = -1 - NumerParams->dIntegParam;

  if (CpnSchedule)
    free(CpnSchedule);
  if (ExeSchedule)
    free(ExeSchedule);

  if (Inst) {
    FreeCalibInst(Inst);
    free(Inst);
  }

  if (AllInst)
    free(AllInst);

  if (HestonParam)
    free(HestonParam);

  if (PricingConst) {
    Free_PricingConst(PricingConst);
    free(PricingConst);
  }

  if (InstNumer) {
    Free_NumerInst(InstNumer);
    free(InstNumer);
  }

  return err;
}

Err LGMSVCmsRateApprox_new(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lFixDate, long lStartDate, long lTheoEndDate, long lPayDate,

    /* Parameter of grids */
    LGMSV_NUMERPARAMS NumerParams,

    /* Output */
    double *dCmsRate) {
  Err err = NULL;

  CalibCpnScheduleDLM *CpnSchedule = NULL;
  CalibExeScheduleDLM *ExeSchedule = NULL;
  CalibInstrumentDLM *Inst = NULL;
  AllCalibInstrumentsDLM *AllInst = NULL;

  SrtCurvePtr yc_ptr;
  long today;

  LGMSV_model *model = NULL;
  LGMSV_HestonInst *HestonParam = NULL;
  LGMSV_PricingConst *PricingConst = NULL;
  LGMSV_NumerInst *InstNumer = NULL;

  double dFixTime;
  long lEndDate;

  char Diag_tenor[20] = "DIAG";
  char *tenor = &(Diag_tenor[0]);
  double dStrike;
  int iNbStrike;

  /* All allocations */
  CpnSchedule = calloc(1, sizeof(CalibCpnScheduleDLM));
  ExeSchedule = calloc(1, sizeof(CalibExeScheduleDLM));
  Inst = calloc(1, sizeof(CalibInstrumentDLM));
  AllInst = calloc(1, sizeof(AllCalibInstrumentsDLM));

  model = calloc(1, sizeof(LGMSV_model));
  HestonParam = calloc(1, sizeof(LGMSV_HestonInst));
  PricingConst = calloc(1, sizeof(LGMSV_PricingConst));
  InstNumer = calloc(1, sizeof(LGMSV_NumerInst));

  if (!CpnSchedule || !ExeSchedule || !Inst || !AllInst || !model ||
      !HestonParam || !PricingConst || !InstNumer) {
    err = "Memory allocation faillure in LGMSVOptionApprox_new";
    goto FREE_RETURN;
  }

  /* Get the curve informations */
  yc_ptr = lookup_curve(yc_name);
  if (!yc_ptr) {
    err = "Yield Curve not found in LGMSVCmsRateApprox_new";
    goto FREE_RETURN;
  }

  today = get_today_from_curve(yc_ptr);

  /* Get the model */
  err = Get_LGMSV_model(und_name, model);

  if (err)
    goto FREE_RETURN;

  /* Get the instrument informations */
  AllocateCalibExeSchedule(&lFixDate, &dFixTime, &lStartDate, &lTheoEndDate,
                           &lEndDate, &dStrike, NULL, NULL, NULL, ExeSchedule);

  /*
  err = Construct_CalibSchedule(	yc_name  ,
                                                                  today  ,
                                                                  swaption_freq
  , swaption_basis  , 1  , &lFixDate  , NULL  , &tenor  , lTheoEndDate  ,
                                                                  &dStrike  ,
                                                                  NULL  ,
                                                                  NULL  ,
                                                                  NULL  ,
                                                                  CpnSchedule  ,
                                                                  ExeSchedule);
                                                                  */

  err = ConstructSingleInstrumentSchedule(
      yc_name, today, swaption_freq, swaption_basis, lFixDate, lStartDate,
      lTheoEndDate, CpnSchedule, ExeSchedule);

  if (err)
    goto FREE_RETURN;

  iNbStrike = 1;

  err = AllocateCalibInst(iNbStrike, Inst);

  if (err)
    goto FREE_RETURN;

  /* Fill the calibration instrument extra informations */

  AllInst->iNbInst = 1;
  AllInst->sCalibInst = Inst;
  Inst->dStrike[0] = 0.05;

  err = Calculate_CalibInst(today, yc_name, NULL, NULL, NULL, swaption_freq,
                            swaption_basis, ref_rate_name, 0, NULL, CpnSchedule,
                            ExeSchedule, AllInst, NULL);

  if (err)
    goto FREE_RETURN;

  /* Set the CMS */
  Inst->iIsCMS = 1;
  Inst->lPayDate = lPayDate;
  Inst->dPayTime = (Inst->lPayDate - today) * YEARS_IN_DAY;
  Inst->dLiborFixTime = dFixTime;

  /* All the precalculation and memory allocation for pricing */
  err = Initialise_PricingConst(model, NumerParams, PricingConst);

  if (err)
    goto FREE_RETURN;

  if (dFixTime < 1.0E-06) {
    *dCmsRate = Inst->dFwdCash + Inst->dSpread;
    goto FREE_RETURN;
  }

  err = Calculate_HestonEquivalent_LGMSV_FromLGM_struct(Inst, CpnSchedule,
                                                        model, HestonParam);

  if (err)
    goto FREE_RETURN;

  err = Calculate_ShiftAndVol_LGMSV_FromLGM_struct(Inst, CpnSchedule, model,
                                                   HestonParam);

  if (err)
    goto FREE_RETURN;

  Inst->iIsLiborOption = 1;

  err = Initialise_NumerInst(model, Inst, HestonParam, InstNumer);

  Inst->iIsLiborOption = 0;

  if (err)
    goto FREE_RETURN;

  Update_NumerInst(Inst, HestonParam, InstNumer);

  /* pricing */
  LGMSVFwdClosedFormApprox_struct(model, Inst, HestonParam, PricingConst,
                                  InstNumer, dCmsRate);

FREE_RETURN:

  if (CpnSchedule)
    free(CpnSchedule);
  if (ExeSchedule)
    free(ExeSchedule);

  if (Inst) {
    FreeCalibInst(Inst);
    free(Inst);
  }

  if (AllInst)
    free(AllInst);

  if (model) {
    free_LGMSV_model(model);
    free(model);
  }

  if (HestonParam)
    free(HestonParam);

  if (PricingConst) {
    Free_PricingConst(PricingConst);
    free(PricingConst);
  }

  if (InstNumer) {
    Free_NumerInst(InstNumer);
    free(InstNumer);
  }

  return err;
}