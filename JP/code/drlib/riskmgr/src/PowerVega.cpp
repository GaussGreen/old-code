//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PowerVega.cpp
//
//   Description : vol shift scenario - add benchmark shift to vol where
//                 shift every point (across strikes) of each maturity by 
//                 Shift * T^-Power where T is the time (in years) between now 
//                 and the relevant maturity, all of them at the same time.
//
//   Author      : Andrew J Swain
//
//   Date        : 28 August 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PowerVega.hpp"

DRLIB_BEGIN_NAMESPACE

PowerVega::Shift::Shift(){} // empty
PowerVega::Shift::~Shift(){} // empty

/** constructor with explicit shift */
PowerVega::PowerVega(double power, double shift):
    ScalarPerturbation(TYPE, shift), power(power) {
}

/** for reflection */
PowerVega::PowerVega(): ScalarPerturbation(TYPE){}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP PowerVega::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool PowerVega::nameMatches(const OutputName&         name,
                           IObjectConstSP          obj){
    // cast obj to PowerVega::Shift and then invoke name method
    const Shift& volShiftObj = dynamic_cast<const Shift&>(*obj);
    return name.equals(volShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void PowerVega::appendName(OutputNameArray&          namesList,
                          IObjectConstSP          obj){
    // cast obj to PowerVega::Shift and then invoke name method
    const Shift& volShiftObj = dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(volShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool PowerVega::shift(IObjectSP obj) {
    // cast obj to PowerVega::Shift and then invoke shift method
    Shift& volShiftObj = dynamic_cast<Shift&>(*obj);
    return volShiftObj.sensShift(this);
}

/** return the power shift (shift * T^-p) */
double PowerVega::powerShift(double years) const {
    return (shiftSize / pow(years, power));
}

class PowerVegaHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PowerVega, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultPowerVega);
        FIELD(power, "power");
    }

    static IObject* defaultPowerVega(){
        return new PowerVega();
    }
};

CClassConstSP const PowerVega::TYPE = CClass::registerClassLoadMethod(
    "PowerVega", typeid(PowerVega), PowerVegaHelper::load);

CClassConstSP const PowerVega::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "PowerVega::Shift", typeid(PowerVega::Shift), 0);


DRLIB_END_NAMESPACE
