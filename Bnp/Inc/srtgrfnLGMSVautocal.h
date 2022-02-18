#ifndef _SRTGRFNLGMSVAUTOCAL_H_
#define _SRTGRFNLGMSVAUTOCAL_H_

#include "srtgrfnLGM2Fautocal.h"
#include "LGMSVCalibApprox.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//	Struct s_grfn_SV_pricing_params to hold Grfn pricing Params for a variety of algorithms. //
//	Observe how the structure contains the basic non-SV equivalent, s_grfn_pricing_params	 //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
typedef struct
{
    s_grfn_pricing_params standard_params;

    long num_steps_vol;
    long num_steps_phi;

} s_grfn_SV_pricing_params;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//	Sort Static level function to calibrate an LGM SV model having either 1 or 2 LGM factors and
// 1 SV	// 	factor, initialise a SV underlying then price a Grfn tableau using either PDE, MC or
// MC with		// 	optimal ex boundary.
//// 						Author:  Paul McCallum				Date:  13th
///Oct 2003
////
//																										//
//  Note the memory management policy in this routine. In particular the fact that we free the
//  Product  // Values when there has been an error but leave the calling routine to free in the
//  absence of errors. //
//	Also how we have the capability to reassign smile values in the same arrays as the trial
// values by	// 	repointing the pointers.
////
//																										//
//	Modified to return primary and secondary model prices.		PMc		2ndDec03
////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
Err GrfnLGMSVAutocalCaller(
    char* YieldCurveName,
    char* VolCubeName,
    Err (*pGetCashVol)(char*, double, double, double, int, char*, double*, double*),
    char*                     DefaultRefRate,
    char*                     PrimSwapFreq,
    char*                     PrimSwapBasis,
    char*                     PrimRefRate,
    int                       NumPrimCalibDates,
    long*                     ExerDatesPrim,
    int*                      DateSelecSpecifiersPrim,
    char**                    CalibTenorsPrim,
    long                      EndDatePrim,
    double**                  PrimStrikes,
    cpd_diag_calib_param*     pPrimParams,
    char*                     SecSwapFreq,
    char*                     SecSwapBasis,
    char*                     SecRefRate,
    int                       NumSecCalibDates,
    long*                     ExerDatesSec,
    int*                      DateSelecSpecifiersSec,
    char**                    CalibTenorsSec,
    long                      EndDateSec,
    double**                  SecStrikes,
    cpd_diag_calib_param*     pSecParams,
    LGMSV_CalibParams*        pCalibParams,
    int                       OneOrTwoFactor,
    double*                   pLambda,
    int*                      pNumSetsSmileParams,
    double**                  pSmileTimes,
    double**                  pAlphaSmileVals,
    double**                  pLambdaSmileVals,
    double**                  pRhoSmileVals,
    double**                  pRho2SmileVals,
    double                    TStar,
    double                    LGMAlpha,
    double                    LGMGamma,
    double                    LGMRho,
    LGMSV_NumerParams*        pNumericalParams,
    int*                      pNumSigmas,
    double**                  pSigmaTimes,
    double**                  pSigmaVals,
    diag_calib_lm_params*     pLmParams,
    long                      Today,
    char*                     UndName,
    LGMSVParam*               pAlgorithmSpecification,
    int                       NumEventDates,
    long*                     EventDates,
    long                      TableauRows,
    long                      TableauCols,
    char***                   TableauStrings,
    int**                     TableauMask,
    long                      AuxWidth,
    long*                     AuxLen,
    double**                  Aux,
    int                       is_end_of_day_fixing,
    int                       is_end_of_day_payment,
    int                       TreeOrMCOrMCoebAsInteger,
    s_grfn_SV_pricing_params* pPricingParams,
    int*                      ExerBoundaryOptimDateSpecifiers,
    MCEBParams*               pExerBoundaryOptimParams,
    int*                      pNumProducts,
    double***                 pProductVals,
    long*                     pNumRowsExerBoundary,
    cpd_calib_inst_data*      pCalibratedInstrumentData,
    double**                  pPrimCalibratedModelPrices,
    double**                  pSecCalibratedModelPrices);

#endif