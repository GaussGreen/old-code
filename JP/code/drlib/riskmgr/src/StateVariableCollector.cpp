//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StateVariableCollector.cpp
//
//   Description : Class for 'collecting' 'true' state variables (actually 
//                 the generators) - sort of a glorified array
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_STATEVARIABLECOLLECTOR_CPP
#include "edginc/StateVariableCollector.hpp"
#include "edginc/StateVariableClient.hpp"

DRLIB_BEGIN_NAMESPACE
IStateVariableCollector::IStateVariableCollector(){}
IStateVariableCollector::~IStateVariableCollector(){}

void StateVariableCollector::append(const IStateVariableGen* svGen) {
    static const string routine = "StateVariableCollector::append";
    
    try {
        // Check for null
        if(!svGen) {
            throw ModelException("Internal error: found Null state variable generator");
        }
    
        // if it's a IStateVariableClient then defer to that object
        const IStateVariableClient* svClient = 
            dynamic_cast<const IStateVariableClient*>(svGen);
        if (svClient){
            svClient->collectStateVars(StateVariableCollectorSP(this));
        } else {
            const IElemStateVariableGen* elemGen = 
                dynamic_cast<const IElemStateVariableGen*>(svGen);
            /** A generator must generate elementary StateVars or declare 
                itself to be an IStateVariableClient */
            if (!elemGen){
                throw ModelException("Internal error: "
                    "unable to figure out state variable generator type");
            }
            elemStateVarGens.push_back(elemGen);
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}
    

void StateVariableCollector::appendElementary(const IElemStateVariableGen* elemGen) {
    static const string routine = "StateVariableCollector::appendElementary";
    
    // Check for null
    if(!elemGen) {
        throw ModelException("Internal error: found Null state variable generator");
    }
    
    elemStateVarGens.push_back(elemGen);
}


const IElemStateVariableGenArray& StateVariableCollector::getElemStateVarGens() const {
    return elemStateVarGens;
}


DRLIB_END_NAMESPACE
