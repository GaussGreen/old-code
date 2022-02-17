//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequest.hpp
//
//   Description : Abstract vol request interface
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOLREQUEST_HPP
#define EDR_VOLREQUEST_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** Captures how to interpolate a volatility, i.e., it captures both
    instrument specific data needed to calculate volatilities/variance
    together with the methodology to be used for the interpolation.
    The methodology is captured by the type of the VolRequest used */
class MARKET_DLL CVolRequest: public CObject {
public:
    static CClassConstSP const TYPE;

    // currently all methods are in the respective base class for each type
    // of vol 
    virtual ~CVolRequest();
protected:
    CVolRequest(const CClassConstSP& clazz);
private:
    CVolRequest(const CVolRequest &rhs);
    CVolRequest& operator=(const CVolRequest& rhs);
};

// smart pointers for CVolRequest
typedef smartConstPtr<CVolRequest> CVolRequestConstSP;
typedef smartPtr<CVolRequest> CVolRequestSP;
#ifndef QLIB_VOLREQUEST_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CVolRequest>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CVolRequest>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CVolRequest>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CVolRequest>);
#endif

// array of vol request
typedef array<CVolRequestSP, CVolRequest> CVolRequestArray;
#ifndef QLIB_VOLREQUEST_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<CVolRequestSP _COMMA_ CVolRequest>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<CVolRequestSP _COMMA_ CVolRequest>);
#endif

// smart pointers for CVolRequestArray
typedef smartConstPtr<CVolRequestArray> CVolRequestArrayConstSP;
typedef smartPtr<CVolRequestArray> CVolRequestArraySP;
#ifndef QLIB_VOLREQUEST_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CVolRequestArray>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CVolRequestArray>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CVolRequestArray>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CVolRequestArray>);
#endif

DRLIB_END_NAMESPACE

#endif
