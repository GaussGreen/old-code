//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolBenchmarkShift.cpp
//
//   Description : vol shift scenario - add Benchmark shift to vol
//
//   Author      : Andrew J Swain
//
//   Date        : 5 June 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBenchmarkShift.hpp"

DRLIB_BEGIN_NAMESPACE
VolBenchmarkShift::Shift::Shift(){} // empty
VolBenchmarkShift::Shift::~Shift(){} // empty

/** constructor with explicit shift */
VolBenchmarkShift::VolBenchmarkShift(ExpirySP expiry, double shift):
    ScalarPerturbation(TYPE, shift), expiry(expiry), isLastShift(true) {}

/** for reflection */
VolBenchmarkShift::VolBenchmarkShift(): 
    ScalarPerturbation(TYPE), isLastShift(true) {}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VolBenchmarkShift::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VolBenchmarkShift::nameMatches(const OutputName&         name,
                                    IObjectConstSP          obj){
    // cast obj to VolBenchmarkShift::Shift and then invoke name method
    const Shift& volShiftObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(volShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VolBenchmarkShift::appendName(OutputNameArray&          namesList,
                                   IObjectConstSP          obj){
    // cast obj to VolBenchmarkShift::Shift and then invoke name method
    const Shift& volShiftObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(volShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VolBenchmarkShift::shift(IObjectSP obj) {
    // cast obj to VolBenchmarkShift::Shift and then invoke shift method
    Shift& volShiftObj =
        dynamic_cast<Shift&>(*obj);
    return volShiftObj.sensShift(this);
}

/** does this match the given expiry ? */
bool VolBenchmarkShift::expiryEquals(const Expiry* expiry) const {
    return this->expiry->equals(expiry);
}

/** is this the last shift for given market data ? */
bool VolBenchmarkShift::lastShift() const {
    return isLastShift;
}

/** store all expiries on vol to be shifted */
void VolBenchmarkShift::cacheExpiries(const ExpiryArrayConstSP& expiries) {
    allExpiries = expiries;
}

ExpiryConstSP VolBenchmarkShift::getExpiry() const {
    return expiry;
}

ExpiryArrayConstSP VolBenchmarkShift::getAllExpiries() const {
    return allExpiries;
}

class VolBenchmarkShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolBenchmarkShift, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultVolBenchmarkShift);
        FIELD(expiry, "expiry");
        FIELD(isLastShift, "isLastShift");
        FIELD_MAKE_OPTIONAL(isLastShift);
        FIELD(allExpiries, "transient");
        FIELD_MAKE_TRANSIENT(allExpiries);
    }

    static IObject* defaultVolBenchmarkShift(){
        return new VolBenchmarkShift();
    }
};

CClassConstSP const VolBenchmarkShift::TYPE = CClass::registerClassLoadMethod(
    "VolBenchmarkShift", typeid(VolBenchmarkShift), VolBenchmarkShiftHelper::load);

CClassConstSP const VolBenchmarkShift::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VolBenchmarkShift::Shift", typeid(VolBenchmarkShift::Shift), 0);

DRLIB_END_NAMESPACE
