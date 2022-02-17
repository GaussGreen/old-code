//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MuSpecial.cpp
//
//   Description : Controls calculation of MU_S
//
//   Author      : Stephen Hope
//
//   Date        : 14 March 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MuSpecial.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
MuSpecial::IShift::~IShift(){} // empty
MuSpecial::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for MU_S */
const string MuSpecial::NAME = "MU_S";
const double MuSpecial::DEFAULT_SHIFT = 0.1;

/** constructor with explicit shift size */
MuSpecial::MuSpecial(double shiftSize, const ExpiryArray* expiries):
    BucketShift(TYPE, NAME, shiftSize, expiries){
}

/** for reflection */
MuSpecial::MuSpecial(): 
    BucketShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double MuSpecial::divisor() const{
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP MuSpecial::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP MuSpecial::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this sensitivity's
    Shift interface */
bool MuSpecial::nameMatches(const OutputName&         name,
                            IObjectConstSP          obj){
    // cast obj to MuSpecial::Shift and then invoke name method
    const IShift& MuSpecialObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(MuSpecialObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    sensitivity's Shift interface */
void MuSpecial::appendName(OutputNameArray&          namesList,
                           IObjectConstSP          obj){
    // cast obj to MuSpecial::Shift and then invoke name method
    const IShift& MuSpecialObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(MuSpecialObj.sensName(this)));
    namesList.push_back(outputName);
}

bool MuSpecial::shift(IObjectSP obj) {
    // cast obj to MuSpecial::IShift and then invoke IShift method
    IShift& muObj = dynamic_cast<IShift&>(*obj);
    return muObj.sensShift(this);
}

void MuSpecial::restore(IObjectSP obj) {
    // cast obj to MuSpecial::IShift and then invoke restore method
    IRestorableShift& muObj = 
        dynamic_cast<IRestorableShift&>(*obj);
    muObj.sensRestore(this);
}


class MuSpecialHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. err none */
    class Factory: public SensitivityFactory::ICreation {};
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MuSpecial, clazz);
        SUPERCLASS(BucketShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultMuSpecial);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(MuSpecial::NAME, 
                                    new Factory(), 
                                    new MuSpecial(),
                                    MuSpecial::IShift::TYPE);
    }

    static IObject* defaultMuSpecial(){
        return new MuSpecial();
    }
};

CClassConstSP const MuSpecial::TYPE = CClass::registerClassLoadMethod(
    "MuSpecial", typeid(MuSpecial), MuSpecialHelper::load);

CClassConstSP const MuSpecial::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "MuSpecial::IShift", typeid(MuSpecial::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(MuSpecial::IRestorableShift, clazz);
    EXTENDS(MuSpecial::IShift);
}
    
CClassConstSP const MuSpecial::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "MuSpecial::IRestorableShift",
    typeid(MuSpecial::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE

