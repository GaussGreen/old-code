//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StateVariableClient.hpp
//
//   Description : Interface that classes must implement if they use State
//                 Variables
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_STATEVARIABLECLIENT_HPP
#define EDR_STATEVARIABLECLIENT_HPP

#include "edginc/StateVariableCollector.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface that classes must implement if they use State Variables */
class RISKMGR_DLL IStateVariableClient {
public:
    IStateVariableClient(); // in StateVariableGen.cpp
    virtual ~IStateVariableClient(); // in StateVariableGen.cpp

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector. Implementations typically call
        IStateVariableCollector::append */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const = 0;
};


DRLIB_END_NAMESPACE

#endif
