//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : DRObject.hpp
//
//   Description : Interface which marks objects which can come through the
//                 DR interface
//
//   Author      : Mark A Robson
//
//   Date        : 24 Nov 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_DROBJECT_HPP
#define EDR_DROBJECT_HPP
#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE
class IDRObject;
typedef smartPtr<IDRObject> IDRObjectSP;
typedef smartConstPtr<IDRObject> IDRObjectConstSP;

#ifndef QLIB_ATOMIC_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<IDRObject>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<IDRObject>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<IDRObject>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<IDRObject>);
#endif

/** Interface which marks objects which can come through the DR interface  */
class TOOLKIT_DLL IDRObject: public virtual IObject{
public:
    static CClassConstSP const TYPE; // defined in Atomic.cpp
    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const = 0;
};


DRLIB_END_NAMESPACE
#endif
