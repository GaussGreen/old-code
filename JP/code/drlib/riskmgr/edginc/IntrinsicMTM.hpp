//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//
//----------------------------------------------------------------------------

#ifndef INTRINSIC_MTM_HPP
#define INTRINSIC_MTM_HPP

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/Model.hpp"
#include "edginc/Control.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL IIntrinsicMTM: virtual public IObject {
public:

    /** calculates MTM intrinsic value and replaces fair value by MTM intrinsic value */
    virtual void calculateIntrinsic(CControlSP control, ResultsSP resuts) = 0;

    static CClassConstSP const TYPE; 
};

typedef smartPtr<IIntrinsicMTM>      IIntrinsicMTMSP;
typedef smartConstPtr<IIntrinsicMTM> IIntrinsicMTMConstSP;

DRLIB_END_NAMESPACE

#endif // INTRINSIC_MTM_HPP
