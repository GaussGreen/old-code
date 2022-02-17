//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StateVariableGen.cpp
//
//   Description : State variable generator
//
//   Date        : Nov 6 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_STATEVARIABLEGEN_CPP
#include "edginc/StateVariableGen.hpp"
#include "edginc/StateVariableClient.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

IStateVariable::IStateVariable(){};
IStateVariable::~IStateVariable(){}
IStateVariableClient::IStateVariableClient() {};
IStateVariableClient::~IStateVariableClient() {};

void StateVarDBase::append(const IStateVariableGen* svGen, IStateVariableSP sv) 
{
    static const string method = "StateVarDBase::append";

#if 0
    // Ensure it does not exist
    DBase::const_iterator iter = svDBase.find(svGen);
    if(iter != svDBase.end()) {
        throw ModelException(method, 
            "State variable with key " +
			     Format::toString((int) (size_t)(svGen)) +
            " already exists. Internal error.");
    }
#endif

    // Append it to the map
    svDBase[svGen] = sv;
}

    
IStateVariableSP StateVarDBase::find(const IStateVariableGen* svGen) const 
{
    static const string method = "StateVarDBase::find";

    // Find the entry in the map
    DBase::const_iterator iter = svDBase.find(svGen);
    if(iter == svDBase.end()) {
        throw ModelException(method, 
            "Unable to find state variable with key " +
			     Format::toString((int) (size_t)(svGen)) +
            ". Internal error.");
    }
    
    // Return the state variable
    return (*iter).second;
}

void StateVarDBase::filterAdvanceables(vector<IAdvanceableStateVariable*>& store)
{
    static const string method = "StateVarDBase::filterAdvanceables";

    for (DBase::iterator iter = svDBase.begin(); iter != svDBase.end(); ++iter)
    {
        IStateVariableSP current = iter->second;
        IAdvanceableStateVariable* tgt = dynamic_cast<IAdvanceableStateVariable*>(current.get());
        if (tgt != 0)
            store.push_back(tgt);
    }
}

DRLIB_END_NAMESPACE
