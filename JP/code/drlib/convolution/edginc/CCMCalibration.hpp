//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMCalibration.hpp
//
//   Description : 
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#ifndef EDR_CCMCALIBRATION_HPP
#define EDR_CCMCALIBRATION_HPP

#include "edginc/CCMConvolution.hpp"


DRLIB_BEGIN_NAMESPACE

/**
 !! This is an implementation file !!
 This file is not supposed to be included by any other header file
 
 This file contains only private functions that are used by the CCM
 correlated recovery model in the fast convolution and the discrete
 convolution algorithms.
*/
class CONVOLUTION_DLL CCMCalibration{
public:
    /** Caluclates the weight on LGD1 as a function of M and recovery
     * parameters */
    static double recoveryWeightCalc(
        double  W1ind,          /* (I) */
        double  T1,             /* (I) */
        bool    isT1Calibrated, /* (I) */
        double  a,              /* (I) */
        double  beta_r,         /* (I) */
        double  qM,             /* (I) */
        double  M);             /* (I) */

    /* Calibration of recovery parameters for one name */
    static void recoveryModelCalibrate(
        const CCMConvolution::NameParam&  b,      /* (I) info on 1 name: survival, loss, correlation */
        double                        pind,   //(I) indep survival
        double                        pgauss, //(I) gaussian copula survival
        bool                          isDisc, //(I) do we calibrate to discretised loss 0=false
        CCMConvolution::LgdParamSP&   p);     //(O) calibrated loss param 
};

DRLIB_END_NAMESPACE
#endif

