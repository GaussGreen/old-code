//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : A generator of FtDEffectiveCurveSPs for the fee and contingent
//               legs, using the information in the instrument and model
//
// Date        : October 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_FTDEFFECTIVECURVEGEN_HPP
#define QLIB_FTDEFFECTIVECURVEGEN_HPP

#include "edginc/ICreditLossGen.hpp"
#include "edginc/MarketFactor.hpp"
#include "edginc/ConvolutionModel.hpp"
#include "edginc/IModelConfigMapper.hpp"
#include "edginc/IEffectiveCurveGen.hpp"
#include "edginc/NToDefaultLossConfig.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(GeneralisedCDO);

class CONVOLUTION_DLL FtDEffectiveCurveGen : virtual public IEffectiveCurveGen  
// JLHP need to implement ICondLossDistributionsGen too, to be able to have 
// FtDs in a tranche
{

public:
    ~FtDEffectiveCurveGen();

    FtDEffectiveCurveGen(NToDefaultLossConfigConstSP ntdLossCfg,
                         IModelConfigMapperConstSP   mapper,
                         IConditionalDefaultsModelSP condDefaultsModel);

    void getFeeAndCtgEffectiveCurve(IDiscountCurveRiskySP&   feeEffCurve, // (O)
                                    IDiscountCurveRiskySP&   ctgEffCurve, // (O)
                                    YieldCurveConstSP        discount,
                                    const IConvolutionModel* convolutionModel,
                                    const DateTime&          lastObservationDate); 

private:
    FtDEffectiveCurveGen(const FtDEffectiveCurveGen& rhs); // don't use
    FtDEffectiveCurveGen& operator=(const FtDEffectiveCurveGen& rhs); // don't use
    
    NToDefaultLossConfigConstSP ftdLossCfg;
    IModelConfigMapperConstSP   modelConfigMapper;
    IConditionalDefaultsModelSP condDefaultsModel;
};

DECLARE_REF_COUNT(FtDEffectiveCurveGen);

DRLIB_END_NAMESPACE

#endif
