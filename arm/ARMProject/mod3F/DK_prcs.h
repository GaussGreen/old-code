/*
 * $Log: DK_prcs.h,v $
 * Revision 1.11  2004/05/10 12:43:29  mcampet
 *  MC from SP add PRCS3F_ConvertObjFwdVolToSpotVol
 *
 * Revision 1.10  2004/05/04 13:31:32  mcampet
 *  MC from SP modif prcs3f_ConvertObjSpotVolToFwdVol
 *
 * Revision 1.9  2004/04/29 07:22:48  mcampet
 * MC from SP add HW3F_Calibration
 *
 * Revision 1.8  2004/04/14 13:30:03  mcampet
 * MC update SP release
 * add ImpliedFwdCorrelation_VFDK_HW1To3F
 *
 * Revision 1.7  2004/03/25 16:54:47  mcampet
 * MC integration SwaptionPrice_VFDK_HW1To3F
 *
 * Revision 1.5  2003/10/21 10:56:03  mab
 * last release from Dimitri
 *
 * Revision 1.4  2003/10/08 13:33:09  mab
 * parameters added to Vol Conversion Func.
 *
 * Revision 1.3  2003/08/18 08:25:47  mab
 * Last Version!
 *
 * Revision 1.2  2003/08/01 13:03:06  mab
 * Improvements
 *
 * Revision 1.1  2003/06/30 16:35:38  mab
 * Initial revision
 *
*/

/*-----------------------------------------------------------------------------------------*/
/*                                                                                         */
/* DK_prcs.h: Header file for intefacing DK PRCS pricing function using ARM objects        */
/*                                                                                         */
/*-----------------------------------------------------------------------------------------*/
#ifndef _DK_PRCS_H
#define _DK_PRCS_H




#include "zeroint.h"
#include "volcurv.h"
#include "volint.h"
#include "linalg.h"





extern ARM_VolLInterpol* PRCS3F_ConvertObjSpotVolToFwdVol( ARM_Date& AsOfDate,
                                                           ARM_ZeroCurve* dBaseRatesNoBasisCurve,
							                               ARM_ZeroCurve* dBaseYieldCurve,
                                                           ARM_ZeroCurve* dForeignRatesNoBasisCurve,
								                           ARM_ZeroCurve* dForeignYieldCurve,
								                           ARM_VolCurve*  volSwopBase,
                                                           ARM_VolCurve*  volSwopForeign,
                                                           ARM_VolLInterpol*  volFx,
						    	                           double dMeanReversionBase, 
							      	                       double dMeanReversionForeign,
								                           ARM_VolCurve* dBaseSpotFXCorrelation,
								                           ARM_VolCurve* dForeignSpotFXCorrelation,
                                                           ARM_VolCurve* dBaseForeignCorrelation,
                                                           ARM_Vector* NoticeDates,
                                                           ARM_Vector* FXCouponResetDates,
                                                           ARM_Vector* FXCouponPaymentDates,
                                                           double CutOff             = 0.0,
                                                           double LongDatedSpotFxVol = 0.0,
                                                           int CalibSwoptBasis = 1,
                                                           ARM_Vector* ForwardVolDates = NULL);


extern ARM_Vector* PRCS3F_Lattice_HWVFDK_Pricing(ARM_Vector* dLatticeGeometryData,
		 										 double dNumTimeLinesBeforeFirstNotice,
					                             double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
												 double dNumTimeLines, 
												 ARM_Date& AsOfDate, 
												 ARM_ZeroCurve* dBaseYieldCurve,
												 ARM_ZeroCurve* dForeignYieldCurve,
												 ARM_ZeroCurve* dBaseRatesNoBasisCurve,
                                                 ARM_ZeroCurve* dForeignRatesNoBasisCurve,
												 ARM_VolCurve* volSwopBase,
												 ARM_VolCurve* volSwopForeign,
												 ARM_VolCurve* volFx,
												 ARM_Vector* dRedemptionData,
												 double dStrike,
												 double dType,
												 ARM_Date& dOptionExpiry, 
												 double dMeanReversionBase, 
												 double dMeanReversionForeign,  
												 double dSpotFX,
												 ARM_VolCurve* dBaseForeignCorrelation,
												 ARM_VolCurve* dBaseSpotFXCorrelation,
												 ARM_VolCurve* dForeignSpotFXCorrelation,
												 double dProductModelCode,
												 double dX1Limit,
												 double dX2Limit,
												 double dX3Limit,
												 double dI1Limit,
												 double dI2Limit,
												 double dI3Limit,
												 double dADPLimit,
												 double dOptimal,
												 double dTimeBoost,
												 double dDeltaFlag,
												 double smoothingValue    = 0.0,
                                                 long calcProbaSurvOrNot  = 0,
				                                 double QBaseSmile        = 0.0,
				                                 double QForeignSmile     = 0.0,
												 double CutOff            = 0,
                                                 double LongDatedSpotFxVol = 0.0,
												 double dStringModel      = 0.0,
												 ARM_Matrix* dBoosterData = NULL,
                                                 int CalibSwoptBasis = 1);

extern ARM_Vector* PRCS3F_Bootstrapping_VFDK_HW1To3F(ARM_VolCurve* volCurve,
		 										     ARM_ZeroCurve* zeroCurve,
												     ARM_Vector* NoticeDates,
                                                     ARM_Vector* StartDates,
                                                     ARM_Vector* EndDates,
                                                     ARM_Vector* HW3FParams,
                                                     ARM_Date& AsOfDate);

extern double PRCS3F_SwaptionPrice_VFDK_HW1To3F(double dSwaptionExpiryInYears,
                                                double dSwaptionTenorInYears,
                                                double dNoticePeriodInDays,
                                                double dStrike,
                                                double dCallPut,
                                                ARM_ZeroCurve* zeroCurve,
                                                ARM_Vector* NoticeDates,
                                                ARM_Vector* Sigma,
                                                ARM_Vector* HW3FParams,
                                                ARM_Date& AsOfDate);

double PRCS3F_ImpliedFwdCorrelation_VFDK_HW1To3F(   double dSwaptionExpiryInYears,
                                                    double dSwaptionTenorInYears,
                                                    double dSwaptionTenor2InYears,
                                                    double dNoticePeriodInDays,
                                                    ARM_ZeroCurve* zeroCurve,
                                                    ARM_Vector* NoticeDates,
                                                    ARM_Vector* Sigma,
                                                    ARM_Vector* HW3FParams,
                                                    ARM_Date& ObservationDate);

ARM_Vector* HW3F_CALIBRATION(   ARM_Date& AsOfDate,
                                ARM_ZeroCurve* zeroCurve,
                                ARM_Vector* HW3FParamsIn,
                                ARM_VolCurve* volCurve,
                                ARM_VolCurve* correlationCurve,
                                ARM_VolCurve* volWeightCurve,
                                ARM_VolCurve* correlationWeightCurve,
                                ARM_Vector* NoticeDates,
                                ARM_Vector* SwapStartDates,
                                ARM_Vector* SwapEndDates);


extern void PRCS3F_FromFwdVolToSpotVol(ARM_Date& AsOfDate,
                                ARM_ZeroCurve* dBaseRatesNoBasisCurve,
                                ARM_ZeroCurve* dBaseYieldCurve,
                                ARM_ZeroCurve* dForeignRatesNoBasisCurve,
                                ARM_ZeroCurve* dForeignYieldCurve,
                                ARM_VolCurve*  volSwopBase,
                                ARM_VolCurve*  volSwopForeign,
                                ARM_VolCurve*  volFx, // This is a Spot Vol!.
                                double dMeanReversionBase, 
                                double dMeanReversionForeign,
                                ARM_VolCurve* dBaseSpotFXCorrelation,
                                ARM_VolCurve* dForeignSpotFXCorrelation,
                                ARM_VolCurve* dBaseForeignCorrelation,
                                ARM_Vector* NoticeDates,
                                ARM_Vector* FXCouponResetDates,
                                ARM_Vector* FXCouponPaymentDates,
                                int    CutOff,
                                double LongDatedSpotFxVol,
                                int CalibSwoptBasis,
                                ARM_Vector* &SpotVolDates,
                                ARM_Vector* &SpotVol);

ARM_VolLInterpol* PRCS3F_ConvertObjFwdVolToSpotVol(ARM_Date& AsOfDate,
                                                   ARM_ZeroCurve* dBaseRatesNoBasisCurve,
                                                   ARM_ZeroCurve* dBaseYieldCurve,
                                                   ARM_ZeroCurve* dForeignRatesNoBasisCurve,
                                                   ARM_ZeroCurve* dForeignYieldCurve,
                                                   ARM_VolCurve*  volSwopBase,
                                                   ARM_VolCurve*  volSwopForeign,
                                                   ARM_VolLInterpol*  volFx, // This is the inputed Fwd Vol!.
                                                   double dMeanReversionBase, 
                                                   double dMeanReversionForeign,
                                                   ARM_VolCurve* dBaseSpotFXCorrelation,
                                                   ARM_VolCurve* dForeignSpotFXCorrelation,
                                                   ARM_VolCurve* dBaseForeignCorrelation,
                                                   ARM_Vector* NoticeDates,
                                                   ARM_Vector* FXCouponResetDates,
                                                   ARM_Vector* FXCouponPaymentDate,
                                                   double CutOff,
                                                   double LongDatedSpotFxVol,
                                                   int CalibSwoptBasis);





#endif
/*-----------------------------------------------------------------------------------------*/
/*---- End Of File ----*/
