/**
 * @file IRiskProperty.cpp
 */

#include "edginc/config.hpp"
#define QLIB_IQUALIFIEDRISKAXIS_CPP
#include "edginc/Atomic.hpp"
#include "edginc/Void.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/ExpiryAndStrike.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/IQualifiedRiskAxis.hpp"

DRLIB_BEGIN_NAMESPACE


template <class QUALIFIER>
IQualifiedRiskAxis<QUALIFIER>::IQualifiedRiskAxis() {}

//template <class QUALIFIER>
//IQualifiedRiskAxis<QUALIFIER>::~IQualifiedRiskAxis() {}

template <class QUALIFIER>
void IQualifiedRiskAxis<QUALIFIER>::load(CClassSP& clazz) {
    REGISTER_INTERFACE(IQualifiedRiskAxis, clazz);
    EXTENDS(IRiskAxis);
}

#ifdef _MSC_VER
# pragma warning( disable : 4661 )
#endif

#define IMPL(T)                                                         \
                                                                        \
    template class RISKMGR_DLL IQualifiedRiskAxis<T>;                   \
                                                                        \
    template <>                                                         \
    CClassConstSP const IQualifiedRiskAxis<T>::TYPE =                   \
        CClass::registerInterfaceLoadMethod(                            \
            "IQualifiedRiskAxis<" #T ">", typeid(IQualifiedRiskAxis<T>), load);

IMPL(Void)
IMPL(ExpiryWindow)
IMPL(ExpiryPair)
IMPL(ExpiryAndStrike)
IMPL(BoxedInt)

DRLIB_END_NAMESPACE
