/**
 * @file LocalCorrVoid.hpp
 */

#ifndef QLIB_LocalCorrVoid_H
#define QLIB_LocalCorrVoid_H

#include "edginc/Void.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

struct LocalCorrVoid;
typedef smartConstPtr<LocalCorrVoid> LocalCorrVoidConstSP;


struct RISKMGR_DLL LocalCorrVoid : CObject {    
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

    LocalCorrVoid(Operation op = PARALLEL_SHIFT);

    static LocalCorrVoidConstSP SP(Operation op);
    
    ~LocalCorrVoid();

    typedef Void Qualifier;

    /** not sure why I need this ... */
    enum { discrete = 0 };
};

#ifndef QLIB_LOCALCORRVOID_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<LocalCorrVoid>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<LocalCorrVoid>);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskProperty<LocalCorrVoid>);
#endif

DRLIB_END_NAMESPACE

#endif
