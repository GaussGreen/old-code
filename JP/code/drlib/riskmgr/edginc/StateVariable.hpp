//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StateVariable.hpp
//
//   Description : Parent class for State Variables
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_STATEVARIABLE_HPP
#define EDR_STATEVARIABLE_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DateTime.hpp"
DRLIB_BEGIN_NAMESPACE

/** Interface implemented by all classes which represent state variables. 
    Note derive virtually from this class to ensure only one refCount field */
class RISKMGR_DLL IStateVariable: public virtual VirtualDestructorBase {
public:
    IStateVariable(); // in StateVariableGen.cpp
    virtual ~IStateVariable();  // in StateVariableGen.cpp

    /** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables - 
        need to see how hard it is to implement */
    virtual bool doingPast() const = 0;
    virtual void prepare(bool mm) {} // this will be overriden by the IQX* state variables; not used by other kinds of SV.
};

#ifndef QLIB_STATEVARIABLEGEN_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<IStateVariable>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<IStateVariable>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<IStateVariable>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<IStateVariable>);
#endif
typedef smartPtr<IStateVariable> IStateVariableSP;


class RISKMGR_DLL IAdvanceableStateVariable : public virtual IStateVariable {
public:
    virtual void advance() = 0;
    virtual void reset() = 0;
};

DRLIB_END_NAMESPACE

#endif
