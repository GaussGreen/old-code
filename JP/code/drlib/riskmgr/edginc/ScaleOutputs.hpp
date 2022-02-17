//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//
//----------------------------------------------------------------------------

#ifndef SCALE_OUTPUTS_HPP
#define SCALE_OUTPUTS_HPP

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/Model.hpp"
#include "edginc/Control.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

/** this class changes an asset to use a risky growth curve rather than the standard yield curve
    primary client of this functionality is the convertible bond model */
class RISKMGR_DLL IScaleOutputs: virtual public IObject {
public:

    /** adds the credit spread to the asset's growth curve */
    virtual void scaleOutputs(CControlSP control, ResultsSP unscaledResults) = 0;

    static CClassConstSP const TYPE; // in Object.cpp
};

// typedef for smart pointers to ISensitiveStrikes
typedef smartPtr<IScaleOutputs> IScaleOutputsSP;
typedef smartConstPtr<IScaleOutputs> IScaleOutputsConstSP;

DRLIB_END_NAMESPACE

#endif // SCALE_OUTPUTS_HPP
