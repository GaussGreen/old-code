/*
 * $Log: bsxtic.h,v $
 * Revision 1.4  2004/03/24 13:43:04  rguillemot
 * SUMMIT Payment Lag
 *
 * Revision 1.3  2004/02/02 11:02:45  mab
 * RCS comment added
 *
 */


/*----------------------------------------------------------------------------*
 
    bsxtic.h
 
    Header for use of options analytic functions
 
*----------------------------------------------------------------------------*/
#ifndef _BSXTIC_H
#define _BSXTIC_H


#include "armglob.h"
#include "dates.h"

#include "zerocurv.h"



extern double bsOption(double spot,
                       double strike,
                       double volatility,
                       double dividend,
                       double discountRate,
                       double maturity,
                       double CallPut);

extern double bsOptionSq(double spot,
                         double strike,
                         double volatility,
                         double dividend,
                         double discountRate,
                         double maturity,
                         double CallPut);

extern double bsDelta(double spot,
                      double strike,
                      double volatility,
                      double dividend,
                      double discountRate,
                      double maturity,
                      double CallPut);

extern double bsKappa(double spot,
                      double strike,
                      double volatility,
                      double dividend,
                      double discountRate,
                      double maturity,
                      double CallPut);

extern double bsGamma(double spot,
                      double strike,
                      double volatility,
                      double dividend,
                      double discountRate,
                      double maturity,
                      double CallPut);

extern double bsVega(double spot,
                     double strike,
                     double volatility,
                     double dividend,
                     double discountRate,
                     double maturity,
                     double CallPut);

extern double bsRho(double spot,
                    double strike,
                    double volatility,
                    double dividend,
                    double discountRate,
                    double maturity,
                    double CallPut);

extern double bsTheta(double spot,
                      double strike,
                      double volatility,
                      double dividend,
                      double discountRate,
                      double maturity,
                      double CallPut);

extern double dekBarrier(double spot,
                         double strike,
                         double barrier,
                         double volatility,
                         double dividend,
                         double discountRate,
                         double maturity,
                         double Rebate,
                         double CallPut,
                         double InOut,
                         double UpDown,
                         int numOptions);

extern double bsBarrier(double spot,
                        double strike,
                        double barrier,
                        double volatility,
                        double dividend,
                        double discountRate,
                        double maturity,
                        double Rebate,
                        double CallPut,
                        double InOut,
                        double UpDown);

extern double bsDoubleBarrier(double PrixSousJacent,
                              double PrixExercice,
                              double UpBarrier,
                              double DownBarrier,
                              double discountRate,
                              double dividend,
                              double volatility,
                              double maturity,
                              double CallPut);

extern double bsBinary(double spot,
                       double strike,
                       double volatility,
                       double dividend,
                       double discountRate,
                       double maturity,
                       double CallPut,
                       double CashAsset);

extern double bsGap(double spot,
                    double strike,
                    double payoff,
                    double volatility,
                    double dividend,
                    double discountRate,
                    double maturity,
                    double CallPut);

extern double bsSuperShare(double spot,
                           double strikeInf,
                           double strikeSup,
                           double volatility,
                           double dividend,
                           double discountRate,
                           double maturity);

extern double bsForwardStart(double spot,
                             double IOMoney,
                             double volatility,
                             double dividend,
                             double discountRate,
                             double startDate,
                             double maturity,
                             double CallPut);

extern double bsPayLater(double spot,
                         double strike,
                         double volatility,
                         double dividend,
                         double discountRate,
                         double maturity,
                         double CallPut);

extern double bsLookBack(double spot,
                         double currentMinOrMax,
                         double volatility,
                         double dividend,
                         double discountRate,
                         double maturity,
                         double CallPut);

extern double bsAveragePrice(double spot,
                             double strike,
                             double numAvrgPoints,
                             double firstAvrgPoints,
                             double dt,
                             double volatility,
                             double dividend,
                             double discountRate,
                             double maturity,
                             double CallPut);

extern double bsChooser(double spot,
                        double callStrike,
                        double putStrike,
                        double volatility,
                        double dividend,
                        double discountRate,
                        double callMaturity,
                        double putMaturity,
                        double callPutDate);

extern double bsCompound(double spot,
                         double volatility,
                         double dividend,
                         double discountRate,
                         double uoStrike,
                         double uoMaturity,
                         double uoCallPut,
                         double oStrike,
                         double oMaturity,
                         double oCallPut);

extern double bsSpot_CallPutParity(double& f, double& df, double spot, 
                                   void** fixedParams);

extern double bsSpot_OptionPrice(double& bf, double& df, double spot, 
                                 void** fixedParams);



class ARM_BSNorModel;
extern double bsVolImpNor(double Price,                
                       double spot,
                       double strike,
                       double dividend,
                       double discountRate,
                       double mat,
                       double CallPut,
                       double Ytol,
					   ARM_BSNorModel* model = NULL);

extern double bsVolImp(double Price,                
                       double spot,
                       double strike,
                       double dividend,
                       double discountRate,
                       double mat,
                       double CallPut,
                       double Ytol,
					   double InitVal = 0.5,
                       int algo = 0);

extern double bsKImpDelta(double delta,
                          double spot,
                          double vol,
                          double dividend,
                          double discountRate,
                          double mat,
                          double CallPut,
                          double Ytol);

extern double bsKImpDelta2(double delta,
                           double spot,
                           double vol,
                           double dividend,
                           double discountRate,
                           double mat,
                           double CallPut,
                           double Ytol);

extern double bsKImp(double Price,
                     double spot,
                     double vol,
                     double dividend,
                     double discountRate,
                     double mat,
                     double CallPut,
                     double Ytol);

extern double CallForwardValue(double F,
                               double K,
                               double Sigma, 
                               double T);

extern double PutForwardValue(double F,
                              double K,
                              double Sigma, 
                              double T);

extern double ComputeFwdBSDelta(double fwd,
                                double strike, 
                                double vol,
                                double T,
                                int CallPut = K_CALL);

extern double ComputeSplinedSigmaATMF(ARM_Vector* deltas,
                                      ARM_Vector* sigmas,
                                      double matu,
                                      double SigmaZDS  = -1,
                                      double Precision = 1e-4,
                                      double FX_SPOT = -1.0);

extern double ComputeDeltaFwdFromDeltaWP(ARM_Date& AsOf,
                                         double matu,
                                         double sigma,
                                         double fxSpot,
                                         double deltaWithPremium,
                                         ARM_ZeroCurve* domCrv, // JPY
                                         ARM_ZeroCurve* foreignCrv); // USD

extern double CalculateImpliedStrikeFromDeltaWithPremium(ARM_Date& AsOf,
                                                  double matu,
                                                  double sigma,
                                                  double fxSpot,
                                                  double deltaWithPremium,
                                                  ARM_ZeroCurve* domCrv, // JPY
                                                  ARM_ZeroCurve* foreignCrv); // USD
/* TMP
extern double CalcFwdFXSpot(ARM_Date& AsOfDate,
							  double Spot,
							  ARM_Date& aFwdDate ,
							  ARM_ZeroCurve* NumDiscountCurve, 
							  ARM_ZeroCurve* UndDiscountCurve); 

*/


//
// New functions For FX VOL
//

extern void ARM_ComputeImpliedZDSWithoutPremiumSpot(double volZDS, double fxFwd, double spot, 
                                                    double matu,
                                                    ARM_ZeroCurve* zcDom,
                                                    ARM_ZeroCurve* zcFor,
                                                    double& strike, 
                                                    double& deltaCall, double& deltaPut);

extern void ARM_ComputeImpliedZDSWithPremiumSpot(double volZDS, double fxFwd, double spot, 
                                                 double matu,
                                                 ARM_ZeroCurve* zcDom,
                                                 ARM_ZeroCurve* zcFor,
                                                 double& strike, 
                                                 double& deltaCall, double& deltaPut);

extern void ARM_ComputeImpliedZDSWithoutPremium(double volZDS, double fxFwd, double matu,
                                                double& strike, double& deltaCall, double& deltaPut);

extern void ARM_ComputeImpliedZDS_FWP(double volZDS, double fxFwd, double matu,
                                      double& strike, double& deltaCall, double& deltaPut);

extern void ARM_ComputeImpliedATMF_FWP(double volATMF, double fxFwd, double matu,
                                       double& strike, double& deltaCall, double& deltaPut);

extern void ARM_ComputeImpliedATMF(double volATMF, double fxFwd, double matu,
                                   double& strike, double& deltaCall, double& deltaPut);


                         /*---- Implied strike calculation ----*/

extern void ARM_ComputeImpliedStrike(double vol, double fxFwd, double matu, 
                                 double target, int callPut,
                                 double& strike);

extern void ARM_ComputeImpliedStrike_FWP(double vol, double fxFwd, double matu, 
                                         double target, int callPut,
                                         double& strike,double firstStrike = -1.0,
									     bool isFirstTime = true);

extern void ARM_ComputeImpliedStrikeFromSpotDeltaWithPremium(double vol, 
                                                             double fxFwd,
                                                             double spot,
                                                             double matu, 
                                                             double target, 
                                                             int callPut,
                                                             ARM_ZeroCurve* zcDom,
                                                             ARM_ZeroCurve* zcFor,
                                                             double& strike);

extern void ARM_ComputeImpliedStrikeFromSpotDeltaWithoutPremium(double vol, 
                                                                double fxFwd,
                                                                double spot,
                                                                double matu, 
                                                                double target, 
                                                                int callPut,
                                                                ARM_ZeroCurve* zcDom,
                                                                ARM_ZeroCurve* zcFor,
                                                                double& strike);

/*------------------------*/

extern double computeDeltaFWP_withStrikeAndVol(double strike, double vol, 
                                               double fxFwd, double matu, int callPut);


extern double computeVol_withDeltaFWPasInput_aux(double strike, double fxFwd, double matu,
                                                 ARM_Vector* deltaCall, ARM_Vector* volsCall,
                                                 ARM_Vector* deltaPut, ARM_Vector* volsPut,
                                                 double pivotStrike, double guessVol,
                                                 int interpolType, int correctSplineWithLinear);



extern double computeVol_withDeltaAsInput_aux(double strike, double fxFwd, double matu,
                                              ARM_Vector* deltaCall, ARM_Vector* volsCall,
                                              ARM_Vector* deltaPut, ARM_Vector* volsPut,
                                              double pivotStrike,
                                              int interpolType);



extern double ComputeVolWithDeltaSpotAsInputWithoutPremium(double strike, 
                                                           double fxFwd,
                                                           double FXSpot,
                                                           double matu,
                                                           ARM_Vector* deltaCall, ARM_Vector* volsCall,
                                                           ARM_Vector* deltaPut, ARM_Vector* volsPut,
                                                           double pivotStrike,
                                                           int interpolType,
                                                           ARM_ZeroCurve* zcDom,
                                                           ARM_ZeroCurve* zcFor);

extern double ComputeVolWithDeltaSpotAsInputWithPremium(double strike, double fxFwd, 
                                                        double FXSpot,
                                                        double matu,
                                                        ARM_Vector* deltaCall, ARM_Vector* volsCall,
                                                        ARM_Vector* deltaPut, ARM_Vector* volsPut,
                                                        double pivotStrike, double guessVol, 
                                                        int interpolType, 
                                                        int correctSplineWithLinear,
                                                        ARM_ZeroCurve* zcDom,
                                                        ARM_ZeroCurve* zcFor);

/*------------------------*/


extern double ComputeSigmaFromDeltaSIMPLEX(double curDelta,
                                           ARM_Vector* delta, ARM_Vector* vols,
                                           int interpolType);


extern double computeDelta_withStrikeAndVol(double strike, double vol, double fxFwd, 
                                            double matu, int callPut);

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
