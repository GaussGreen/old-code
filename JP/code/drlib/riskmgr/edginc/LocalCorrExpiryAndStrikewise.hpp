/**
 * @file LocalCorrExpiryAndStrikewise.hpp
 */

#ifndef QLIB_LocalCorrExpiryAndStrikewise_H
#define QLIB_LocalCorrExpiryAndStrikewise_H

#include "edginc/ExpiryAndStrike.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

struct RISKMGR_DLL LocalCorrExpiryAndStrikewise : CObject {
    static CClassConstSP const TYPE;
    LocalCorrExpiryAndStrikewise();
    ~LocalCorrExpiryAndStrikewise();

    typedef ExpiryAndStrike Qualifier;

    enum {
        discrete = 0
    };
};

#ifndef QLIB_LOCALCORREXPIRYANDSTRIKEWISE_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<LocalCorrExpiryAndStrikewise>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<LocalCorrExpiryAndStrikewise>);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskProperty<LocalCorrExpiryAndStrikewise>);
#endif

DRLIB_END_NAMESPACE

#endif
