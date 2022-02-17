/**
 * @file LocalCorrExpiry.hpp
 */

#ifndef QLIB_LocalCorrExpiry_H
#define QLIB_LocalCorrExpiry_H

#include "edginc/ExpiryWindow.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

struct LocalCorrExpiry;
typedef smartConstPtr<LocalCorrExpiry> LocalCorrExpiryConstSP;


struct RISKMGR_DLL LocalCorrExpiry : CObject {
    static CClassConstSP const TYPE;

    /** 
     * PARALLEL_SHIFT
     *      additive shift up or down (depending on shift only)
     *      shifted_squeeze(strike) = squeeze(strike) + shift
     *  
     * SKEWED_SHIFT
     *      additive shift up or down (depending on both, shift & strike)
     *      shifted_squeeze(strike) = squeeze(strike) + strike * shift 
     */

    /** read by LocalCorrSqueeze::sensShift */
    enum Operation { PARALLEL_SHIFT, SKEWED_SHIFT };

    Operation op;

    LocalCorrExpiry(Operation op = PARALLEL_SHIFT);

    static LocalCorrExpiryConstSP SP(Operation op);

    ~LocalCorrExpiry();

    typedef ExpiryWindow Qualifier;

    /** not sure why I need this ... */
    enum { discrete = 0 };
};

#ifndef QLIB_LOCALCORREXPIRY_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<LocalCorrExpiry>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<LocalCorrExpiry>);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskProperty<LocalCorrExpiry>);
#endif

DRLIB_END_NAMESPACE

#endif
