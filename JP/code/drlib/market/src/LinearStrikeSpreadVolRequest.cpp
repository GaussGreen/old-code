//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LinearStrikeSpreadVolRequest.cpp
//
//   Description : Linear Strike with spread vol request 
//   Author      : Andrew J Swain
//
//   Date        : 14 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LinearStrikeSpreadVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

LinearStrikeSpreadVolRequest::LinearStrikeSpreadVolRequest(
    double          interpLevel,
    const DateTime& startDate,
    const DateTime& endDate,
    double          spread): 
    VolRequestLNStrike(TYPE), 
    interpLevel(interpLevel < 0.0? 0.0: interpLevel), startDate(startDate), 
    endDate(endDate), spread(spread) {}


/** Scale the interpolation object */
void LinearStrikeSpreadVolRequest::scale(double scaleFactor){
    // don't bother, it's a moneyness
}


/** Returns the start date for forward starting volatility interpolation. */
DateTime LinearStrikeSpreadVolRequest::getStartDate() const{
    return startDate;
}

/** Returns the start date for forward starting volatility interpolation. */
const DateTime& LinearStrikeSpreadVolRequest::getEndDate() const{
    return endDate;
}

/** Validates that this vol request is configured correctly for forward
    starting/spot starting. If fwdStartExpected is true and the request
    corresponds to an absolute strike, for example, an exception will
    be thrown. Similarly if fwdStartExpected is false and the request
    holds a % strike. */
void LinearStrikeSpreadVolRequest::validateFwdStart(
    bool fwdStartExpected) const{
    if (!fwdStartExpected){
        throw ModelException("LinearStrikeSpreadVolRequest::getInterpLevel",
                             "Inconsistent forward start/spot starting data.\n"
                             "Instrument is doing forward starting vol "
                             "interpolation\nbut VolRequest is"
                             " flagged as not forward starting");
    }
}

/** Returns the start date for forward starting volatility interpolation. */
double LinearStrikeSpreadVolRequest::getInterpLevel(bool expectPerc) const{
    validateFwdStart(expectPerc);
    return interpLevel;
}

/** return the spread */
double LinearStrikeSpreadVolRequest::getSpread() const {
    return spread;
}

/** Returns the sensitive strike */
void LinearStrikeSpreadVolRequest::sensitiveStrikes(
    double               spot,
    const DoubleArraySP& strikes) const
{
    double sensStrike = getInterpLevel(true);

    strikes->push_back(spot);
    strikes->push_back(sensStrike * spot);
}

/** sets the strike of this request using the supplied strike */
void LinearStrikeSpreadVolRequest::setStrike(double strike)  {
    interpLevel = strike;
}

/** Returns the strike as a percentage. If the strike is held internally
    as a level, the supplied spot is used to turn into a percentage */
double LinearStrikeSpreadVolRequest::getPercStrike(double spot) const{
    return interpLevel; // always a %
}

/* for reflection */
LinearStrikeSpreadVolRequest::LinearStrikeSpreadVolRequest(): 
    VolRequestLNStrike(TYPE), interpLevel(0){}

class LinearStrikeSpreadVolRequestHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LinearStrikeSpreadVolRequest, clazz);
        SUPERCLASS(VolRequestLNStrike);
        EMPTY_SHELL_METHOD(defaultLinearStrikeSpreadVolRequest);
        FIELD(interpLevel,   "Interpolation level");
        FIELD(startDate,     "Start of vol interpolation");
        FIELD(endDate,       "End of vol interpolation");
        FIELD(spread,        "Fwd start vol spread");
    }

    static IObject* defaultLinearStrikeSpreadVolRequest(){
        return new LinearStrikeSpreadVolRequest();
    }
};

CClassConstSP const LinearStrikeSpreadVolRequest::TYPE = 
CClass::registerClassLoadMethod(
    "LinearStrikeSpreadVolRequest", typeid(LinearStrikeSpreadVolRequest), 
    LinearStrikeSpreadVolRequestHelper::load);


DRLIB_END_NAMESPACE

