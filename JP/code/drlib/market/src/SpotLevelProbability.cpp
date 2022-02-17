//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotLevelProbability.cpp
//
//   Description : spot level scenario - for probability supplied, use the
//                 vol of the asset to determine the spot level that encloses
//                 this confidence interval
//
//   Author      : Andrew J Swain
//
//   Date        : 16 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SpotLevelProbability.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE
SpotLevelProbability::Shift::~Shift(){} // empty

/** what the asset should set its spot to */
double SpotLevelProbability::spotLevel(const DateTime& today, 
                                       const Asset* asset) const {
    static const string method = "SpotLevelProbability::spotLevel";
    try {
        double stdDev = N1Inverse(probability);
        double spot   = asset->getSpot();
        double strike = spot * moneyness;

        // where we should interpolate vol - no fwd start, straight strike
        LinearStrikeVolRequest volRequest(strike, today, today, false);

        CVolProcessedBSSP vol(asset->getProcessedVol(&volRequest));
        
        // when we should interpolate
        MaturityTimePeriod period(maturity, DateTime::timeConvert(time));
        DateTime endDate = period.toDate(today);

        double variance = vol->CalcVar(today, endDate);

        double newSpot = spot*exp(-0.5*variance + stdDev*sqrt(variance));

        return newSpot;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** validation code to be called after object construction */
void SpotLevelProbability::validatePop2Object() {
    static const string method("SpotLevelProbability::validatePop2Object");

    if (!Maths::isPositive(probability) || probability >= 1.0) {
        throw ModelException(method,
                             "probability ("+
                             Format::toString(probability) + 
                             ") must lie in (0, 1)");
    }

    if (Maths::isNegative(moneyness)) {
        throw ModelException(method,
                             "moneyness ("+
                             Format::toString(moneyness) + 
                             ") must be >= 0.0");
    }      
}


/** for reflection */
SpotLevelProbability::SpotLevelProbability(): 
    Perturbation(TYPE), probability(0.0), moneyness(0.0) {
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP SpotLevelProbability::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool SpotLevelProbability::nameMatches(const OutputName&         name,
                                       IObjectConstSP            obj){
    // cast obj to SpotLevelProbability::Shift and then invoke name method
    const Shift& spotLevelObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(spotLevelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void SpotLevelProbability::appendName(OutputNameArray&          namesList,
                                      IObjectConstSP            obj){
    // cast obj to SpotLevelProbability::Shift and then invoke name method
    const Shift& spotLevelObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(spotLevelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool SpotLevelProbability::shift(IObjectSP obj) {
    // cast obj to SpotLevelProbability::Shift and then invoke shift method
    Shift& spotLevelObj =
        dynamic_cast<Shift&>(*obj);
    return spotLevelObj.sensShift(this);
}

class SpotLevelProbabilityHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SpotLevelProbability, clazz);
        SUPERCLASS(Perturbation);
        EMPTY_SHELL_METHOD(defaultSpotLevelProbability);
        FIELD(probability, "probability");
        FIELD(moneyness, "moneyness");
        FIELD(maturity, "maturity");
        FIELD(time, "time");
        FIELD(toTweak, "ignored - do not use");
        FIELD_MAKE_OPTIONAL(toTweak);
    }

    static IObject* defaultSpotLevelProbability(){
        return new SpotLevelProbability();
    }
};

CClassConstSP const SpotLevelProbability::TYPE = CClass::registerClassLoadMethod(
    "SpotLevelProbability", typeid(SpotLevelProbability), SpotLevelProbabilityHelper::load);

CClassConstSP const SpotLevelProbability::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "SpotLevelProbability::Shift", typeid(SpotLevelProbability::Shift), 0);


DRLIB_END_NAMESPACE
