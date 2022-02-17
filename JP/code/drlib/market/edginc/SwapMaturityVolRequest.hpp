//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SwapMaturityVolRequest.hpp
//
//   Description : Interest rate vol volrequest where interp at given
//                 swap maturity
//
//   Author      : Andrew J Swain
//
//   Date        : 7 November 2001
//
//
//----------------------------------------------------------------------------

#ifndef _SWAPMATURITYVOLREQUEST_HPP
#define _SWAPMATURITYVOLREQUEST_HPP

#include "edginc/Expiry.hpp"
#include "edginc/VolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interest rate vol vol request supporting differing style of
    calibration type */
class MARKET_DLL SwapMaturityVolRequest: public CVolRequest {
public:
    static CClassConstSP const TYPE;

    static const string CALIB_CMS; // "CMS"
    static const string CALIB_FIX; // "FIX"

    virtual ~SwapMaturityVolRequest();

    /** Interest rate vol vol request where interp is done as CALIB_CMS */
    SwapMaturityVolRequest(const Expiry* maturity);

    /** Interest rate vol vol request where interp is done as specified.
        Currently CALIB_CMS or CALIB_FIX */
    SwapMaturityVolRequest(const Expiry* maturity, const string& calibType);
    //// Returns the expiry used to 'interpolate' with
    ExpiryConstSP swapMaturity() const;
    //// Returns the style of the calibration (currently CALIB_CMS or CALIB_FIX)
    const string& calibType() const; 
private:
    friend class SwapMaturityVolRequestHelper;
    SwapMaturityVolRequest();
    SwapMaturityVolRequest(const SwapMaturityVolRequest& rhs);
    SwapMaturityVolRequest& operator=(const SwapMaturityVolRequest& rhs);

    // fields
    ExpirySP maturity;
    string   calib;
};

typedef smartConstPtr<SwapMaturityVolRequest> SwapMaturityVolRequestConstSP;
typedef smartPtr<SwapMaturityVolRequest> SwapMaturityVolRequestSP;

DRLIB_END_NAMESPACE

#endif
