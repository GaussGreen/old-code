/**
 * @file Void.hpp
 */

#ifndef DRLIB_Void_H
#define DRLIB_Void_H

#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
/** Class mainly exists to make life easy for templates where you may want
    a field of varying type. For templates where you don't want the parameter
    then you make the type of it Void. Note that Void fields are make
    transient automatically and hence are hidden from the client interface. */
class TOOLKIT_DLL Void: public CObject {
public:
    Void();
    virtual ~Void();
    static CClassConstSP const TYPE;
};
DECLARE(Void);

#ifndef QLIB_VOID_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<Void>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<Void>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<Void>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<Void>);
#endif

DRLIB_END_NAMESPACE

#endif
