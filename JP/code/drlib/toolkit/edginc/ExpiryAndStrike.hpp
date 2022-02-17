/**
 * @file ExpiryAndStrike.hpp
 */

#ifndef QLIB_ExpiryAndStrike_H
#define QLIB_ExpiryAndStrike_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ExpiryAndStrike)

#ifndef QLIB_EXPIRYANDSTRIKE_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<ExpiryAndStrike>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<ExpiryAndStrike>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<ExpiryAndStrikeSP _COMMA_ ExpiryAndStrike>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<ExpiryAndStrike>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<ExpiryAndStrike>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<ExpiryAndStrikeSP _COMMA_ ExpiryAndStrike>);
#endif

class TOOLKIT_DLL ExpiryAndStrike: public CObject {

    ExpiryAndStrike();
    static IObject* emptyShell();
    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    ExpiryConstSP expiry;
    double strike;

public:

    ExpiryAndStrike(ExpiryConstSP expiry, double strike);
    static ExpiryAndStrikeConstSP SP(ExpiryConstSP expiry, double strike);

    ~ExpiryAndStrike();

    static ExpiryAndStrikeArrayConstSP series(ExpiryArrayConstSP expiries,
                                              const DoubleArray& strikes);
};

DRLIB_END_NAMESPACE

#endif
