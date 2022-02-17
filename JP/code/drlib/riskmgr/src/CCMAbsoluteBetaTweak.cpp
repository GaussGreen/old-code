/**
 * @file CCMAbsoluteBetaTweak.cpp
 *
 *
 */

#include "edginc/config.hpp"
#include "edginc/CCMAbsoluteBetaTweak.hpp"
#include "edginc/RestorableWith.hpp"

DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<CCMAbsoluteBetaTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<CCMAbsoluteBetaTweak>",
        typeid(TweakableWith<CCMAbsoluteBetaTweak>),
        TweakableWith<CCMAbsoluteBetaTweak>::load);

template<> CClassConstSP const RestorableWith<CCMAbsoluteBetaTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<CCMAbsoluteBetaTweak>",
        typeid(RestorableWith<CCMAbsoluteBetaTweak>),
        RestorableWith<CCMAbsoluteBetaTweak>::load);

DRLIB_END_NAMESPACE

