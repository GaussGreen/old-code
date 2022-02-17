//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//  Description  : Handles the adjustment of delta shift size 
//                 (to make sure we don't spill over bounding strikes)
//
//
//----------------------------------------------------------------------------

#ifndef VOL_DELTA_SHIFT_SIZE_HPP
#define VOL_DELTA_SHIFT_SIZE_HPP

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IVolDeltaShiftSize: virtual public IObject {
public:

    /** handles the delat shift size adjustment */
    virtual void adjustDeltaShiftSize(ShiftSizeCollector* collector,
                                      const string assetName,
                                      double spot) const = 0;

    static CClassConstSP const TYPE; 
};

// typedef for smart pointers to ISensitiveStrikes
typedef smartConstPtr<IVolDeltaShiftSize> IVolDeltaShiftSizeConstSP;
typedef smartPtr<IVolDeltaShiftSize> IVolDeltaShiftSizeSP;

DRLIB_END_NAMESPACE

#endif // VOL_DELTA_SHIFT_SIZE_HPP
