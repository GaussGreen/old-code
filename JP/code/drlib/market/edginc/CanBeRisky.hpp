//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//
//----------------------------------------------------------------------------

#ifndef CAN_BE_RISKY_HPP
#define CAN_BE_RISKY_HPP

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/Model.hpp"
#include "edginc/CreditCurve.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

/** this class changes an asset to use a risky growth curve rather than the standard yield curve
    primary client of this functionality is the convertible bond model */
class MARKET_DLL ICanBeRisky: virtual public IObject {
public:

    /** adds the credit spread to the asset's growth curve */
    virtual void makeRisky(ICreditCurveSP creditSpreads,
                        const  DateTime   *maturityDate) = 0;

    static CClassConstSP const TYPE; // in Object.cpp
};

// typedef for smart pointers to ISensitiveStrikes
typedef smartConstPtr<ICanBeRisky> ICanBeRiskyConstSP;
typedef smartPtr<ICanBeRisky> ICanBeRiskySP;

DRLIB_END_NAMESPACE

#endif // CAN_BE_RISKY_HPP
