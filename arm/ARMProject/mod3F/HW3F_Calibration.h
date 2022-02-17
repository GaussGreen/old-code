/*
 * $Log: HW3F_Calibration.h,v $
 * Revision 1.3  2004/05/25 13:30:18  mab
 * Added : Calibration control
 *
 * Revision 1.2  2004/04/29 07:29:04  mcampet
 * MC Updates from SP
 *
 * Revision 1.1  2004/04/22 12:08:07  mcampet
 * Initial revision
 *
 */


/*----------------------------------------------------------------------------*/

#ifdef unix
#include <ieeefp.h>
#else
#define finite _finite
#endif



#include "DKMaille.h"
#include "DKMaille2D.h"


#define MAX_3F_PARAM    100 // Max dimension of X

#define CONSTANT_X_ITER 4   // After how many MINIMUM iterations we compare
                            // X to X_PREC

#define THE_3F_PRECISION 0.01


/*----------------------------------------------------------------------------*/


class HW3F_Calibration_Params
{
    public :

        double SpotDate;
        DKMaille<double> ZCDates;
        DKMaille<double> ZCRates;
        DKMaille<double> SwaptionPaymentDates;
        DKMaille<double> WhichModelParametersToCalibrate;
        DKMaille<double> ModelParameters;
        DKMaille<double> InputArray_X_AbsVols;
        DKMaille<double> InputArray_Y_AbsVols;
        DKMaille2D<double> InputSurface_AbsVols;
        DKMaille<double> InputArray_X_FwdCorrelation_2_20;
        DKMaille<double> InputArray_FwdCorrelation_2_20;
        DKMaille2D<double> InputSurface_AbsVolsWeights;
        DKMaille<double> InputArray_FwdCorrelation_2_20Weights;
        int uiIterations;
        DKMaille<double> PreviousModelParameters;
        DKMaille<double> NoticeDatesBis;
        DKMaille<double> SwapStartDatesBis;
        DKMaille<double> SwapEndDatesBis;
        int Plantage;
        DKMaille<double> BestParams;
        double BestResult;
 
        // For managing special slow convergence cases!

		double X_PREC[200]; // Previous iteration Point

		int    iterCounter; // "Global" iterations counter
};




double HW3F_Calibration_Function_TO_Minimize(double* dLocalModelParameters);


DKMaille<double> HW3F_CALIBRATION(double dSpotDate,
                                  DKMaille<double> &dZCDates,
                                  DKMaille<double> &dZCRates,
                                  DKMaille<double> &dModelParameters,
                                  DKMaille<double> &dInputArray_X_AbsVols,
                                  DKMaille<double> &dInputArray_Y_AbsVols,
                                  DKMaille2D<double> &dInputSurface_AbsVols,
                                  DKMaille<double> &dInputArray_X_FwdCorrelation_2_20,
                                  DKMaille<double> &dInputArray_FwdCorrelation_2_20,
                                  DKMaille2D<double> &dInputSurface_AbsVolsWeights,
                                  DKMaille<double> &dInputArray_FwdCorrelation_2_20Weights,
                                  DKMaille<double> &dNoticesDates,
                                  DKMaille<double> &dSwapStartDates,
                                  DKMaille<double> &dSwapEndDates
                                  //double MaxTime
                                  );

void Fill(double dSpotDate,
          DKMaille<double> &dZCDates,
          DKMaille<double> &dZCRates,
          DKMaille<double> &dModelParameters,
          DKMaille<double> &dWhichModelParametersToCalibrate,
          DKMaille<double> &dInputArray_X_AbsVols,
          DKMaille<double> &dInputArray_Y_AbsVols,
          DKMaille2D<double> &dInputSurface_AbsVols,
          DKMaille<double> &dInputArray_X_FwdCorrelation_2_20,
          DKMaille<double> &dInputArray_FwdCorrelation_2_20,
          DKMaille2D<double> &dInputSurface_AbsVolsWeights,
          DKMaille<double> &dInputArray_FwdCorrelation_2_20Weights,
          DKMaille<double> &dNoticesDates,
          DKMaille<double> &dSwapStartDates,
          DKMaille<double> &dSwapEndDates);


double GetSwaptionData( double dXX, 
                        double dYY,
                        DKMaille<double> &dSwaptionDates,
                        DKMaille<double> &dAccrualPeriods,
                        DKMaille<double> &dBasisSwaptionDates,
                        DKMaille<double> &dBasisAccrualPeriods,
                        DKMaille<double> &dBasis,
                        double dFixedPaymentFrequency,
                        double dFloatPaymentFrequency);



/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/
