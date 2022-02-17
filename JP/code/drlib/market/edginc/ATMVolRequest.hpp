//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ATMVolRequest.hpp
//
//   Description : "At the money" vol request 
//
//   Author      : Mark A Robson
//
//   Date        : 12 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ATM_VOLREQUEST_HPP
#define EDG_ATM_VOLREQUEST_HPP

#include "edginc/DateTime.hpp"
#include "edginc/VolRequestLN.hpp"

DRLIB_BEGIN_NAMESPACE

/** This class represents a volatility interpolation algorithm which
    always interpolates "at the money" ie at spot (precise meaning will
    depend on volatility representation) */
class MARKET_DLL ATMVolRequest: public CVolRequestLN {
public:
    static CClassConstSP const TYPE;
    friend class ATMVolRequestHelper;
    ATMVolRequest();

    /** Scale the interpolation object */
    virtual void scale(double scaleFactor); // multiplicative scaling

    /** Returns the start date for forward starting volatility 
        interpolation. This current returns DateTime() - we need to review
        whether we need this method for all vol interps once the XCB
        stuff is working */
    virtual DateTime getStartDate() const;

    /** Returns whether the vol request is for a forward starting 
        vol interpolation */
    virtual bool isForwardStarting() const;

    /** Returns the sensitive strike */
    virtual void getSensitiveStrike(double spot,
                                    CDoubleArraySP sensitiveStrikes) const;

    /** Validates that this vol request is configured correctly for forward
        starting/spot starting. If fwdStartExpected is true and the request
        corresponds to an absolute strike, for example, an exception will
        be thrown. Similarly if fwdStartExpected is false and the request
        holds a % strike. */
    void validateFwdStart(bool fwdStartExpected) const;

protected:
private:
    ATMVolRequest(const ATMVolRequest &rhs);
    ATMVolRequest& operator=(const ATMVolRequest& rhs);
};

typedef smartConstPtr<ATMVolRequest> ATMVolRequestConstSP;
typedef smartPtr<ATMVolRequest> ATMVolRequestSP;
#ifndef QLIB_ATMVOLREQUEST_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<ATMVolRequest>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<ATMVolRequest>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<ATMVolRequest>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<ATMVolRequest>);
#endif

DRLIB_END_NAMESPACE

#endif
