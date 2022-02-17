//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRelativeShift.cpp
//
//   Description : Scenario shift where a relative shift (e.g. a 10% shift
//                 to 20% vol gives 22% vol) is applied for each benchmark.
//                 Shifts are defined for ranges of benchmarks (e.g. <= 1Y,
//                 1Y -> 5Y, 5Y ->  30Y etc).
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

VolRelativeShift::IShift::~IShift(){} // empty
VolRelativeShift::IShift::IShift(){} // empty

/** for reflection */
VolRelativeShift::VolRelativeShift(): MultiExpiryStepShift(TYPE), 
    moneyness(1.0), gotVol(false) {}

// what's the shift for a given date ? */
double VolRelativeShift::shiftSize(const DateTime& today,  // to anchor expiries
                                   const DateTime& shiftDate) const {
    static const string method("VolRelativeShift::shiftSize");
    try {
        if (!gotVol) {
        throw ModelException(method, "moneyness vol not set");            
        }

        // shift is ABSOLUTE based on a RELATIVE proportion of the 
        // moneyness vols
        double relative = MultiExpiryStepShift::shiftSize(today, shiftDate);
        double basevol  = vol->CalcVol(today, shiftDate);
        double absolute = basevol * relative;

        return absolute;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VolRelativeShift::shiftInterface() const{
    return IShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VolRelativeShift::nameMatches(const OutputName&     name,
                                   IObjectConstSP        obj){
    // cast obj to VolRelativeShift::Shift and then invoke name method
    const IShift& volShiftObj = dynamic_cast<const IShift&>(*obj);
    return name.equals(volShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VolRelativeShift::appendName(OutputNameArray&       namesList,
                                  IObjectConstSP         obj){
    // cast obj to VolRelativeShift::Shift and then invoke name method
    const IShift& volShiftObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(volShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VolRelativeShift::shift(IObjectSP obj) {
    // cast obj to VolRelativeShift::Shift and then invoke shift method
    IShift& volShiftObj = dynamic_cast<IShift&>(*obj);
    return volShiftObj.sensShift(this);
}

void VolRelativeShift::validatePop2Object() {
    try {
        MultiExpiryShift::validatePop2Object();
    }
    catch (exception& e) {
        throw ModelException(e, "VolRelativeShift::validatePop2Object");
    }
}
  

/** Stores the spot value of the associated asset */
void VolRelativeShift::setSpot(const DateTime& today, const Asset* asset) {
    static const string method("VolRelativeShift::setSpot");
    try {

        // Clone the asset before requesting processed vol, otherwise the
        // vol surface which is shifted *is* the surface that is simultaneously
        // interrogated for moneyness vols
        AssetConstSP myAsset(copy(asset));

        double spot = myAsset->getSpot();
        // need to grab the moneyness vols here as might not have 
        // access later on
        LinearStrikeVolRequest request(spot*moneyness, today, today, false);

        vol = CVolProcessedBSSP(myAsset->getProcessedVol(&request));

        gotVol = true;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Override clone method to copy our extra data over */
IObject* VolRelativeShift::clone() const{
    // first clone all the registered fields
    IObject*  copy = CObject::clone();
    VolRelativeShift* vrs = dynamic_cast<VolRelativeShift*>(copy);
    if (!vrs){
        throw ModelException("VolRelativeShift::clone"); // shouldn't happen
    }
    
    vrs->vol = vol;
    return copy;
}

class VolRelativeShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Relative shifts to vol by benchmark");
        REGISTER(VolRelativeShift, clazz);
        SUPERCLASS(MultiExpiryStepShift);
        EMPTY_SHELL_METHOD(defaultVolRelativeShift);
        FIELD(moneyness, "moneyness");
        FIELD_MAKE_OPTIONAL(moneyness);
        FIELD_NO_DESC(gotVol);
        FIELD_MAKE_TRANSIENT(gotVol);
    }

    static IObject* defaultVolRelativeShift(){
        return new VolRelativeShift();
    }
};

CClassConstSP const VolRelativeShift::TYPE = CClass::registerClassLoadMethod(
    "VolRelativeShift", typeid(VolRelativeShift), VolRelativeShiftHelper::load);

CClassConstSP const VolRelativeShift::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VolRelativeShift::IShift", typeid(VolRelativeShift::IShift), 0);


DRLIB_END_NAMESPACE
