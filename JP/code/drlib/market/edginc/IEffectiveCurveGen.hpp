//----------------------------------------------------------------------------
//                                                                           
// Group       : Credit Hybrids QR                                    
//                                                                           
// Description : Interface defining what a generator for 'effective 
//               curves' based engine needs to be able to do.
//                                                                           
// Date        : October 2006                                                   
//                                                                           
//----------------------------------------------------------------------------

#ifndef QLIB_IEFFECTIVECURVEGEN_HPP
#define QLIB_IEFFECTIVECURVEGEN_HPP

#include "edginc/Control.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ICreditLossGen.hpp"
#include "edginc/FORWARD_DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(IDiscountCurveRisky);
FORWARD_DECLARE(IConvolutionModel);

/** What a generator for an 'effective  curve' based engine needs to be
    able to do */
class MARKET_DLL IEffectiveCurveGen : virtual public ICreditLossGen {

public:
    virtual ~IEffectiveCurveGen();

    virtual void getFeeAndCtgEffectiveCurve(
        IDiscountCurveRiskySP&   feeEffCurve, // (O)
        IDiscountCurveRiskySP&   ctgEffCurve, // (O)
        YieldCurveConstSP        discount,
        const IConvolutionModel* convolutionModel,
        const DateTime&          lastObservationDate) = 0;

protected:
    IEffectiveCurveGen();

private:    
    IEffectiveCurveGen(const IEffectiveCurveGen& rhs); // don't use
    IEffectiveCurveGen& operator=(const IEffectiveCurveGen& rhs); // don't use
};

DECLARE_REF_COUNT(IEffectiveCurveGen);

DRLIB_END_NAMESPACE

#endif
