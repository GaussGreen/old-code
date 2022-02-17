//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LocVolRequest.cpp
//
//   Description : Local vol request 
//
//   Author      : Regis Guichard
//
//   Date        : 24 April 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

#define SPEEDUP_METHOD_NONE     "none"

LocVolRequest::LocVolRequest(
    DateTime        startDate,
    bool            fwdStarting,
    bool            useUnscaledTimeTweak,
    bool            useNextStepDerivs,
    bool            useMidPoint,
    double          strikeTweakUnscaled,
    double          timeTweakUnscaled,
    double          probDensRatioMin): CVolRequestDVF(TYPE), 
    startDate(startDate),
    fwdStarting(fwdStarting),
    strikeTweakUnscaled(strikeTweakUnscaled),
    timeTweakUnscaled(timeTweakUnscaled),
    probDensRatioMin(probDensRatioMin),
    useUnscaledTimeTweak(useUnscaledTimeTweak),
    useNextStepDerivs(useNextStepDerivs),
    useMidPoint(useMidPoint) {}

LocVolRequest::LocVolRequest(
    DateTime        startDate,
    bool            fwdStarting,
    bool            useUnscaledTimeTweak,
    bool            useNextStepDerivs,
    bool            useMidPoint,
    double          strikeTweakUnscaled,
    double          timeTweakUnscaled,
    double          probDensRatioMin,
    string          speedup): CVolRequestDVF(TYPE), 
    startDate(startDate),
    fwdStarting(fwdStarting),
    strikeTweakUnscaled(strikeTweakUnscaled),
    timeTweakUnscaled(timeTweakUnscaled),
    probDensRatioMin(probDensRatioMin),
    useUnscaledTimeTweak(useUnscaledTimeTweak),
    useNextStepDerivs(useNextStepDerivs),
    useMidPoint(useMidPoint),
    speedup(speedup) {}

/** Returns the start date for forward starting volatility interpolation. */ 
DateTime LocVolRequest::getStartDate() const{
    return startDate;
}

bool LocVolRequest::isFwdStarting() const{

    return fwdStarting;
}

double LocVolRequest::getStrikeTweakUnscaled() const{
    return strikeTweakUnscaled;
}

double LocVolRequest::getTimeTweakUnscaled() const{
    return timeTweakUnscaled;
}

double LocVolRequest::getProbDensRatioMin() const{
    return probDensRatioMin;
}
bool LocVolRequest::getUseUnscaledTimeTweak() const{
    return useUnscaledTimeTweak;
}
bool LocVolRequest::getUseNextStepDerivs() const{
    return useNextStepDerivs;
}
bool LocVolRequest::getUseMidPoint() const{
    return useMidPoint;
}
string LocVolRequest::getSpeedup() const{
    return speedup;
}

/* for reflection */
LocVolRequest::LocVolRequest(): 
    CVolRequestDVF(TYPE),
    fwdStarting(false),
    strikeTweakUnscaled(0.005),
    timeTweakUnscaled(0.001),
    probDensRatioMin(0.01),
    useUnscaledTimeTweak(true),
    useNextStepDerivs(false),
    useMidPoint(true),
    speedup(SPEEDUP_METHOD_NONE){}

class LocVolRequestHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LocVolRequest, clazz);
        SUPERCLASS(CVolRequestDVF);
        EMPTY_SHELL_METHOD(defaultLocVolRequest);
        FIELD(startDate, "Start date");
        FIELD(fwdStarting, "is vol interpolation fwd starting");
        FIELD(strikeTweakUnscaled, "strikeTweakUnscaled");
        FIELD(timeTweakUnscaled, "timeTweakUnscaled");
        FIELD(probDensRatioMin, "probDensRatioMin");
        FIELD(useUnscaledTimeTweak, "useUnscaledTimeTweak");
        FIELD(useNextStepDerivs, "useNextStepDerivs");
        FIELD(useMidPoint, "useMidPoint");
        FIELD(speedup, "speedup");
        FIELD_MAKE_OPTIONAL(speedup); 


        Addin::registerConstructor("LOC_VOL_REQUEST",
                                   Addin::MARKET,
                                   "Creates a handle to a linear strike vol request",
                                   LocVolRequest::TYPE);
    }

    static IObject* defaultLocVolRequest(){
        return new LocVolRequest();
    }
};

CClassConstSP const LocVolRequest::TYPE = 
CClass::registerClassLoadMethod("LocVolRequest", typeid(LocVolRequest), LocVolRequestHelper::load);

DRLIB_END_NAMESPACE








