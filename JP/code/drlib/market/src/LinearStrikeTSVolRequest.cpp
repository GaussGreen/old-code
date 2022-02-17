//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LinearStrikeTSVolRequest.cpp
//
//   Description : Linear Strike vol request which incorporates a spread
//                 term structure for forward starting options
//
//   Author      : Stephen Hope
//
//   Date        : 23 April 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

LinearStrikeTSVolRequest::LinearStrikeTSVolRequest(
    double          interpLevel,
    const DateTime& startDate,
    const DateTime& endDate,
    bool            isFwdStarting): 
    VolRequestLNStrike(TYPE), 
    interpLevel(interpLevel < 0.0? 0.0: interpLevel), startDate(startDate), 
    endDate(endDate), isFwdStarting(isFwdStarting){}


/** Scale the interpolation object */
void LinearStrikeTSVolRequest::scale(double scaleFactor){
    if (!isFwdStarting)
    {
        interpLevel *= scaleFactor;
    }
}


/** Returns the start date for forward starting volatility interpolation. */
DateTime LinearStrikeTSVolRequest::getStartDate() const{
    return startDate;
}

/** Returns the start date for forward starting volatility interpolation. */
const DateTime& LinearStrikeTSVolRequest::getEndDate() const{
    return endDate;
}

/** Validates that this vol request is configured correctly for forward
    starting/spot starting. If fwdStartExpected is true and the request
    corresponds to an absolute strike, for example, an exception will
    be thrown. Similarly if fwdStartExpected is false and the request
    holds a % strike. */
void LinearStrikeTSVolRequest::validateFwdStart(bool fwdStartExpected) const{
    if (fwdStartExpected != isFwdStarting){
        static const string method("LinearStrikeTSVolRequest::getInterpLevel");
        string m("Inconsistent forward start/spot starting data.\n");
        if (isFwdStarting){
            m.append("Instrument is not doing forward starting vol "
                     "interpolation\nbut VolRequest is"
                     " flagged as forward starting");
        } else {
            m.append("Instrument is doing forward starting vol "
                     "interpolation\nbut VolRequest is"
                     " flagged as not forward starting");
        }
        throw ModelException(method, m);
    }
}

/** Returns the interpolation level. The expectPerc must match up with
    the isFwdStarting flag */
double LinearStrikeTSVolRequest::getInterpLevel(bool expectPerc) const{
    validateFwdStart(expectPerc);
    return interpLevel;
}

/** Returns the sensitive strike */
void LinearStrikeTSVolRequest::sensitiveStrikes(
    double               spot,
    const DoubleArraySP& strikes) const
{
    double sensitiveStrike = interpLevel;

    if (isFwdStarting) {
        strikes->push_back(spot);
        sensitiveStrike *= spot;
    }

    strikes->push_back(sensitiveStrike);
}

/** sets the strike of this request using the supplied strike */
void LinearStrikeTSVolRequest::setStrike(double strike)  {
    interpLevel = strike;
}
                                        
/** Returns the strike as a percentage. If the strike is held internally
    as a level, the supplied spot is used to turn into a percentage */
double LinearStrikeTSVolRequest::getPercStrike(double spot) const{
    return (isFwdStarting? interpLevel: interpLevel/spot);
}


/* for reflection */
LinearStrikeTSVolRequest::LinearStrikeTSVolRequest(): 
    VolRequestLNStrike(TYPE), interpLevel(0){}

class LinearStrikeTSVolRequestHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LinearStrikeTSVolRequest, clazz);
        SUPERCLASS(VolRequestLNStrike);
        EMPTY_SHELL_METHOD(defaultLinearStrikeTSVolRequest);
        FIELD(interpLevel, "Interpolation level");
        FIELD(startDate, "Start date");
        FIELD(endDate, "End date");
        FIELD(isFwdStarting, "is vol interpolation fwd starting");

        Addin::registerConstructor("LINEAR_STRIKE_TS_VOL_REQUEST",
                                   Addin::MARKET,
                                   "Creates a handle to a linear strike vol request with term structure",
                                   LinearStrikeTSVolRequest::TYPE);
    }

    static IObject* defaultLinearStrikeTSVolRequest(){
        return new LinearStrikeTSVolRequest();
    }
};

CClassConstSP const LinearStrikeTSVolRequest::TYPE = 
CClass::registerClassLoadMethod(
    "LinearStrikeTSVolRequest", typeid(LinearStrikeTSVolRequest), 
    LinearStrikeTSVolRequestHelper::load);


DRLIB_END_NAMESPACE
