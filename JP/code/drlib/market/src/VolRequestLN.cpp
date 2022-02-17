//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequest.cpp
//
//   Description : Abstract vol request interface
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VOLREQUESTLN_CPP
#include "edginc/VolRequestLN.hpp"

DRLIB_BEGIN_NAMESPACE

CVolRequestLN::~CVolRequestLN(){}

CVolRequestLN::CVolRequestLN(const CClassConstSP& clazz): 
    CVolRequest(clazz), negFwdVarAllowed(false) {}

/** Returns the end date if applicable */
const DateTime& CVolRequestLN::getEndDate () const
{
    throw ModelException("CVolRequestLN::getEndDate", 
                         "End date does not exist for this vol request");
}

/** Returns true if this vol request asks for negative forward variance
    to be allowed. Default is false */
bool CVolRequestLN::negativeFwdVarAllowed() const{
    return negFwdVarAllowed;
}

/** sets whether this vol request asks for negative forward variance
    to be allowed */
void CVolRequestLN::allowNegativeFwdVar(bool allow){
    negFwdVarAllowed = allow;
}

/** throws an exception if negativeFwdVarAllowed() is true */
void CVolRequestLN::checkNegativeFwdVarDisallowed(const string& method) const{
    if (negFwdVarAllowed){
        throw ModelException(method, "Vol Request is allowing negative fwd "
                             "var but this is not supported here");
    }
}

void CVolRequestLN::load(CClassSP& clazz){
    REGISTER(CVolRequestLN, clazz);
    SUPERCLASS(CVolRequest);
    FIELD(negFwdVarAllowed, "Is negative forward variance allowed");
    FIELD_MAKE_TRANSIENT(negFwdVarAllowed); /* default false - we really don't
                                               want this in the public iface */
}

CClassConstSP const CVolRequestLN::TYPE = CClass::registerClassLoadMethod(
    "VolRequestLN", typeid(CVolRequestLN), load);

// initialise type for array of CVolRequestLN
DEFINE_TEMPLATE_TYPE_WITH_NAME("VolRequestLNArray", CVolRequestLNArray);

DRLIB_END_NAMESPACE


