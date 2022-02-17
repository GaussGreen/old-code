//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SwapMaturityVolRequest.cpp
//
//   Description : Interest rate vol vol request where interp at given
//                 swap maturity
//
//   Author      : Andrew J Swain
//
//   Date        : 7 November 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE
const string SwapMaturityVolRequest::CALIB_CMS("CMS");
const string SwapMaturityVolRequest::CALIB_FIX("FIX");

SwapMaturityVolRequest::~SwapMaturityVolRequest(){}

/** Interest rate vol vol request where interp is done as CALIB_CMS */
SwapMaturityVolRequest::SwapMaturityVolRequest(const Expiry* maturity): 
    CVolRequest(TYPE), maturity(copy(maturity)), calib(CALIB_CMS) {}

/** Interest rate vol vol request where interp is done as specified.
    Currently CALIB_CMS or CALIB_FIX */
SwapMaturityVolRequest::SwapMaturityVolRequest(
    const Expiry* maturity, const string& calibType): 
    CVolRequest(TYPE), maturity(copy(maturity)), calib(calibType) {
    if (calibType != CALIB_CMS && calibType != CALIB_FIX){
        throw ModelException("SwapMaturityVolRequest::SwapMaturityVolRequest",
                             "Calibration type "+calibType+" not supported");
    }
}

//// Returns the expiry used to 'interpolate' with
ExpiryConstSP SwapMaturityVolRequest::swapMaturity() const {
    return maturity;
}

//// Returns the style of the calibration (currently CALIB_CMS or CALIB_FIX)
const string& SwapMaturityVolRequest::calibType() const{
    return calib;
}

/* for reflection */
SwapMaturityVolRequest::SwapMaturityVolRequest(): CVolRequest(TYPE) {}

class SwapMaturityVolRequestHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SwapMaturityVolRequest, clazz);
        SUPERCLASS(CVolRequest);
        EMPTY_SHELL_METHOD(defaultSwapMaturityVolRequest);
        FIELD(maturity, "maturity");
        FIELD(calib, "Calibration type");
        Addin::registerConstructor("SWAP_MATURITY_VOL_REQUEST",
                                   Addin::MARKET,
                                   "Creates a handle to a swap maturity"
                                   " vol request",
                                   SwapMaturityVolRequest::TYPE);
    }

    static IObject* defaultSwapMaturityVolRequest(){
        return new SwapMaturityVolRequest();
    }
};

CClassConstSP const SwapMaturityVolRequest::TYPE = 
CClass::registerClassLoadMethod(
    "SwapMaturityVolRequest", typeid(SwapMaturityVolRequest), 
    SwapMaturityVolRequestHelper::load);

DRLIB_END_NAMESPACE
