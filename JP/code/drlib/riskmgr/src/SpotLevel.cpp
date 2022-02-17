//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotLevel.hpp
//
//   Description : spot level scenario - set spot to supplied value
//
//   Author      : Andrew J Swain
//
//   Date        : 25 April 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SpotLevel.hpp"

DRLIB_BEGIN_NAMESPACE
SpotLevel::Shift::~Shift(){} // empty

SpotLevel::SpotLevel(CClassConstSP clazz, double spot):
    ScalarPerturbation(clazz, spot), floorNegFwdPrice(false){
    validatePop2Object();
}
    

/** constructor with explicit spot level */
SpotLevel::SpotLevel(double spot):
    ScalarPerturbation(TYPE, spot), floorNegFwdPrice(false){
    validatePop2Object();
}

/** constructor with explicit spot level and floorNegFwdPrice*/
SpotLevel::SpotLevel(double spot, bool floorNegFwdPrice):
    ScalarPerturbation(TYPE, spot), floorNegFwdPrice(floorNegFwdPrice){
    validatePop2Object();
}

/** for reflection */
SpotLevel::SpotLevel(): ScalarPerturbation(TYPE), floorNegFwdPrice(false){}

/** Validation */
void SpotLevel::validatePop2Object(){
    // Testing positivity with Maths::isPositive() caused
    // testing/creditinp/neg_fwd.xml to fail since it generates a 2.0773e-176
    // entry in spotGrid in VanillaCreditSupport::preProcess() and does
    // AssetUnd::setSpot with it ... so I cowardly did plain "> 0", just
    // in case there are any production implications.

    if (!(getShiftSize() > 0)) {
        throw ModelException(
            "SpotLevel::validatePop2Object",
            "SpotLevel shiftSize is illegally set to " +
                Format::toString(shiftSize) +
                ": it must be greater than zero (NB it specifies the absolute "
                "level to which the spot is set, not a shift)");
    }
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP SpotLevel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool SpotLevel::nameMatches(const OutputName&         name,
                            IObjectConstSP          obj){
    // cast obj to SpotLevel::Shift and then invoke name method
    const Shift& spotLevelObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(spotLevelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void SpotLevel::appendName(OutputNameArray&          namesList,
                       IObjectConstSP          obj){
    // cast obj to SpotLevel::Shift and then invoke name method
    const Shift& spotLevelObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(spotLevelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool SpotLevel::shift(IObjectSP obj) {
    // cast obj to SpotLevel::Shift and then invoke shift method
    Shift& spotLevelObj =
        dynamic_cast<Shift&>(*obj);
    return spotLevelObj.sensShift(this);
}

// given today's date, when is spot defined
DateTime SpotLevel::spotDate(const DateTime& today) const {
    return today;
}


bool SpotLevel::zeroNegFwd() const {
    return floorNegFwdPrice;
}

class SpotLevelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SpotLevel, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultSpotLevel);
        FIELD(floorNegFwdPrice, "what to use for setWillZeroNegFwd");
        FIELD_MAKE_OPTIONAL(floorNegFwdPrice);
    }

    static IObject* defaultSpotLevel(){
        return new SpotLevel();
    }
};

CClassConstSP const SpotLevel::TYPE = CClass::registerClassLoadMethod(
    "SpotLevel", typeid(SpotLevel), SpotLevelHelper::load);

CClassConstSP const SpotLevel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "SpotLevel::Shift", typeid(SpotLevel::Shift), 0);

DRLIB_END_NAMESPACE
