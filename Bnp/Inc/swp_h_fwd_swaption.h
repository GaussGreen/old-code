/* ===================================================================================
   FILENAME:      swp_h_swaption.h

   PURPOSE:       Compute a fwd swaption with the CMS methos compatible with the
   smile
   ===================================================================================
 */

#ifndef SWP_H_FWD_SWAPTION_H
#define SWP_H_FWD_SWAPTION_H

#include "swp_h_all.h"

Err PriceFwdSwaption(char *cYCname, char *cRefRname, long xlMaturityDate,
                     long xlStartDate, long xlEndDate, SrtCompounding srtFreq,
                     SrtBasisCode srtBasis, int nStrikes, double *Strikes,
                     SrtCallPutType srtCallPut, int xlStudDegree, int xlNumSim,
                     char *xlTypeVol, int xlNPts, double xlUpBound,
                     int xlnClasses, double *SigmaHeston, double *AlphaHeston,
                     double *LambdaHeston, double *BetaHeston,
                     double *RhoHeston,

                     double *SigmaHestonShort, double *AlphaHestonShort,
                     double *LambdaHestonShort, double *BetaHestonShort,
                     double *RhoHestonShort, double **CorrelLong,
                     double *dPrice, double *dPriceSwnShort,
                     double *dPriceSwnLong);

Err CalibrateCorrelation(char *cYCname, char *cRefRname, int xlMaturityDate,
                         int nDates, long *lDates, double *dCoverages,
                         double dCoverageShort, double dFwd,

                         SrtCompounding srtFreq, SrtBasisCode srtBasis,
                         int nStrikes, double *Strikes, double *BSVols,

                         SrtCallPutType srtCallPut, int xlStudDegree,
                         int xlNumSim, char *xlTypeVol, int xlNPts,
                         double xlUpBound, int xlnClasses,

                         double *SigmaHeston, double *AlphaHeston,
                         double *LambdaHeston, double *BetaHeston,
                         double *RhoHeston,

                         double **CorrelLong, double **CorrelOutput);

Err Fwd_Swaption_Get_Cum_Density(Date lToday, char *cYCname, char *cRefRname,
                                 SrtCompounding srtFreq, SrtBasisCode srtBasis,

                                 int iNPayDatesLong, Date *lPayDatesLong,
                                 double *dCoveragesLong,

                                 double *SigmaHeston, double *AlphaHeston,
                                 double *LambdaHeston, double *ShiftHeston,
                                 double *RhoHeston,

                                 int iNPts, double dUpBound, int nClass,

                                 double **dPPDx, double **dPPDy,
                                 double **dCPDy);

Err GetDensity_Lvl_Measure(int nClass, double dMat, double dFwd, double dSigmaH,
                           double dAlphaH, double dSigmaIftyH, double dMeanRevH,
                           double dShiftH, double dRhoH, double dDisc,
                           double dUpBound, int nSteps, int nTen, double dCov,
                           double **dPPDx, double **dPPDy, double **dCPDy);

Err GetSimulatedLevel(double *DfSet, int iNDates, Date *lEndDates,
                      double *dCoverages, double *dLevel);

Err GetSpanningSetDf(double **dCopulaCube, int iNDates, Date *lEndDates,
                     double *dCoverages, int iMC, double *DfSet);

Err GetFwdSwaptionPayoff(double **dCopulaCube, int iNDatesLong,
                         Date *lEndDatesLong, double *dCoveragesLong,

                         int iNDatesShort, Date *lEndDatesShort,
                         double *dCoveragesShort,

                         int nStrikes, double *xlStrike, int iMC,
                         double *dPayoff, double *dPayoffSwnShort,
                         double *dPayoffSwnLong, int flag);

Err FindDecayRate(double x1, double x2, double tol,

                  double Maturity, int iNDates, Date *lEndDates,
                  double *dCoverages, double dCoverageShort,

                  int nStrikes, double *Strikes, double *BSVols, double dFwd,

                  double **dCopulaCube,

                  double **dPPDx, double **dCPDy, int xlStudDegree,
                  double *dMean_v, int nClass, int nSim,

                  double *SigmaHestonShort, double *AlphaHestonShort,
                  double *LambdaHestonShort, double *BetaHestonShort,
                  double *RhoHestonShort,

                  double **dCorrel, double *dDecay);

Err EvalChi2Error(double dTryDecay, double Maturity, int iNDates,
                  Date *lEndDates, double *dCoverages,

                  double dCoverageShort, int nStrikes, double *Strikes,
                  double *BSVols, double dFwd,

                  double **dCopulaCube, double **dPPDx, double **dCPDy,
                  int xlStudDegree, double *dMean_v, int nClass, int nSim,
                  double *SigmaHestonShort, double *AlphaHestonShort,
                  double *LambdaHestonShort, double *BetaHestonShort,
                  double *RhoHestonShort, double **Correl, double *ChiError);

Err GetCapletPayoff(double **dCopulaCube, int iNDates, Date *lEndDates,
                    double *dCoverages, double dCoverageShort,

                    int nStrikes, double *Strikes, int iMC, double *dPayoff);

#endif