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

#ifndef EDR_VOL_REQUEST_LN_HPP
#define EDR_VOL_REQUEST_LN_HPP

#include "edginc/DateTime.hpp"
#include "edginc/VolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

/** Captures how to interpolate a Log-Normal volatility, i.e., it captures both
    instrument specific data needed to calculate volatilities/variance
    together with the methodology to be used for the interpolation.
    The methodology is captured by the type of the VolRequest used */
class MARKET_DLL CVolRequestLN: public CVolRequest {
public:
    static CClassConstSP const TYPE;
    /** Scale the interpolation object */
    virtual void scale(double scaleFactor) = 0; // multiplicative scaling

    /** Returns the start date for forward starting volatility interpolation */
    virtual DateTime getStartDate() const = 0;

    /** Returns the end date if applicable. Default is to fail - the end
        date is currently required by the Future asset. To be reviewed */
    virtual const DateTime& getEndDate () const;

    /** Returns the sensitive strike */
    virtual void getSensitiveStrike(double spot,
                                    CDoubleArraySP sensitiveStrikes) const = 0;

    /** Returns true if this vol request asks for negative forward variance
        to be allowed. Default is false */
    bool negativeFwdVarAllowed() const;

    /** sets whether this vol request asks for negative forward variance
        to be allowed */
    void allowNegativeFwdVar(bool allow);

    /** throws an exception using supplied string as the method name 
        if negativeFwdVarAllowed() is true */
    void checkNegativeFwdVarDisallowed(const string& method) const;

    virtual ~CVolRequestLN();
protected:
    CVolRequestLN(const CClassConstSP& clazz);
private:
    bool negFwdVarAllowed; // transient field

    static void load(CClassSP& clazz);
    CVolRequestLN(const CVolRequestLN &rhs);
    CVolRequestLN& operator=(const CVolRequestLN& rhs);
};

// smart pointers for CVolRequestLN
typedef smartConstPtr<CVolRequestLN> CVolRequestLNConstSP;
typedef smartPtr<CVolRequestLN> CVolRequestLNSP;
// array of vol request
typedef array<CVolRequestLNSP, CVolRequestLN> CVolRequestLNArray;
#ifndef QLIB_VOLREQUESTLN_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CVolRequestLN>);
EXTERN_TEMPLATE(class MARKET_DLL array<CVolRequestLNSP _COMMA_ CVolRequestLN>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CVolRequestLN>);
INSTANTIATE_TEMPLATE(class MARKET_DLL array<CVolRequestLNSP _COMMA_ CVolRequestLN>);
#endif



// smart pointers for CVolRequestLNArray
typedef smartConstPtr<CVolRequestLNArray> CVolRequestLNArrayConstSP;
typedef smartPtr<CVolRequestLNArray> CVolRequestLNArraySP;

DRLIB_END_NAMESPACE

#endif
