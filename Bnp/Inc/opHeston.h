#ifndef OPHESTON_H
#define OPHESTON_H

#include "utallhdr.h"
#include <UTERROR.H>
#include <UTGREEKS.H>
#include <UTTYPES.H>

/* ========================================================================== */

#define HESTON_CONST 0.06
/*
PRICING ADD_INS
*/

Err ConvertSABRParamInHestonParam(double Forward, double Maturity,
                                  double SABRBetaVol, double SABRAlpha,
                                  double SABRMeanReversion, double SABRBeta,
                                  double SABRRho, double *HESTONVol,
                                  double *HESTONAlpha,
                                  double *HESTONMeanReversion,
                                  double *HESTONBeta, double *HESTONRho);

Err HestonPrice(double Forward, double *Strike, int nStrikes, double Maturity,
                double Sigma, double Alpha, double SigmaIfty, double b,
                double Beta, double Rho, double Disc, double UpperBound,
                SrtCallPutType call_put, SrtGreekType greek,
                SRT_Boolean isVolInfFix, int IntegrType, int nSteps,
                double *result);

/* Auxiliary and mathematical functions */

double HestonFindRightNode(double ImAppr, double Strike, double NodeLeft,
                           SRT_Boolean isVolInfFix, double Maturity, double a,
                           double b, double Alpha, double Rho, double Forward,
                           double Sigma, double U, SrtGreekType greek);

Err PsiFunction(SRT_Boolean isVolInfFix, double T, double U, double V, double a,
                double b, double alpha, double rho, double x1, double x2,
                SrtGreekType greek, double *RePsi, double *ImPsi,
                double *DRePsi, double *DImPsi, double *DDRePsi,
                double *DDImPsi);

/*
VOLATILITY CONVERSION ADD_INS
*/

Err HestonATMVol(double Forward, double Maturity, double Sigma, double Alpha,
                 double SigmaIfty, double b, double Beta, double Rho,
                 double Disc, double UpperBound, int nSteps,
                 SRT_Boolean isVolInfFix, double *result);

Err HestonVol(double Forward, double *Strikes, int nStrikes, double Maturity,
              double Vol, double VolInfty, double Alpha, double Beta,
              double Rho, double MeanRev, SrtVolConversion TypeConversion,
              double UpperBound, int nSteps, int IntegerType,
              SRT_Boolean isVolInfFix, double *NewVols);

Err srt_opthestonvol(double Forward, double *Strikes, int nStrikes,
                     double Maturity, double Sigma, double Alpha, double Beta,
                     double Rho, double Lambda, SrtDiffusionType input,
                     SrtDiffusionType output, double *NumericalParams,
                     int nParams, double *OutputVols);

Err srt_opthestonvol2(double Forward, double Strike, double Maturity,
                      double Sigma, double Alpha, double Beta, double Rho,
                      double Lambda, double *NumericalParams, int nParams,
                      SrtDiffusionType input, SrtDiffusionType output,
                      double *OutputVol);

/*
CALIBRATION ADD_INS
*/

Err CalibrateHestonModel(const double dForward, const double dBeta,
                         const double dTau, const double dExpiry,
                         const double dUpperBound, const int iNumSteps,
                         const int iIntegerType, const SRT_Boolean isVolInfFix,
                         const SrtCalibrationType CalibrType,
                         const double *CalibrDblArgs,
                         const double *CalibrInitArgs, double *CalibratedArgs);

Err HestonGetShift(double Forward, double Maturity, double SigmaBetaSABR,
                   double AlphaSABR, double BetaSABR, double RhoSABR,
                   double ATMVol, double *Beta);

double SmileChi2(double Forward, double Maturity, double Sigma0,
                 double AlphaHeston, double SigmaIftyHeston,
                 double LambdaHeston, double BetaHeston, double RhoHeston,
                 double UpperBound, int nSteps, int IntegerType, int nK,
                 double *StrikeVec, double SigmaBeta, double AlphaSABR,
                 double BetaSABR, double RhoSABR, double gap,
                 SRT_Boolean isVolInfFix);

Err HestonCalibrateSigmaInit(double Forward, double Maturity, double ATMVol,
                             double AlphaHeston, double SigmaIftyHeston,
                             double LambdaHeston, double dBeta,
                             double RhoHeston, double UpperBound, int nSteps,
                             int iIntegerType, SRT_Boolean isVolInfFix,
                             double *result);

Err HestonCalibrRho(const double dForward, const double dExpiry,
                    const double *CalibrDblArgs, const double AlphaHeston,
                    const double SigmaIftyHeston, const double LambdaHeston,
                    const double dBeta, const double SigmaHeston,
                    const double dUpperBound, const int iNumSteps,
                    const int iIntegerType, const SRT_Boolean isVolInfFix,
                    double *RhoHeston);

Err HestonCalibrAlpha(const double dForward, const double dExpiry,
                      const double *CalibrDblArgs, const double SigmaHeston,
                      const double SigmaIftyHeston, const double LambdaHeston,
                      const double dBeta, const double RhoHeston,
                      const double dUpperBound, const int iNumSteps,
                      const int iIntegerType, const SRT_Boolean isVolInfFix,
                      double *AlphaHeston);

Err CalibrateHestonToSABR(double Forward, double Maturity, double SigmaBeta,
                          double AlphaSABR, double BetaSABR, double RhoSABR,
                          double *SigmaHeston, double *AlphaHeston,
                          double *SigmaIftyHeston, double *LambdaHeston,
                          double *BetaHeston, double *RhoHeston,
                          double UpperBound, int nSteps, int iIntegerType,
                          double nStdDev, SRT_Boolean isVolInfFix,
                          SrtCalibrationType CalibrType);

#endif