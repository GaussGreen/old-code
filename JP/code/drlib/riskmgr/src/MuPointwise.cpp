//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MuPointwise.cpp
//
//   Description : Controls calculation of MU_POINTWISE
//
//   Author      : Stephen Hope
//
//   Date        : 20 March 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MuPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
MuPointwise::IShift::~IShift(){} // empty
MuPointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for MU_POINTWISE */
const string MuPointwise::NAME = "MU_POINTWISE";
const double MuPointwise::DEFAULT_SHIFT = 0.1;

/** constructor with explicit shift size */
MuPointwise::MuPointwise(double shiftSize, const ExpiryArray* expiries):
    BucketShift(TYPE, NAME, shiftSize, expiries){
}

/** for reflection */
MuPointwise::MuPointwise(): 
    BucketShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double MuPointwise::divisor() const{
    return getShiftSize() * 100.0;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP MuPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP MuPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this sensitivity's
    Shift interface */
bool MuPointwise::nameMatches(const OutputName&         name,
                            IObjectConstSP          obj){
    // cast obj to MuPointwise::Shift and then invoke name method
    const IShift& MuPointwiseObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(MuPointwiseObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    sensitivity's Shift interface */
void MuPointwise::appendName(OutputNameArray&          namesList,
                           IObjectConstSP          obj){
    // cast obj to MuPointwise::Shift and then invoke name method
    const IShift& MuPointwiseObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(MuPointwiseObj.sensName(this)));
    namesList.push_back(outputName);
}

bool MuPointwise::shift(IObjectSP obj) {
    // cast obj to MuPointwise::IShift and then invoke IShift method
    IShift& muObj =
        dynamic_cast<IShift&>(*obj);
    return muObj.sensShift(this);
}

void MuPointwise::restore(IObjectSP obj) {
    // cast obj to MuPointwise::IShift and then invoke restore method
    IRestorableShift& muObj =
        dynamic_cast<IRestorableShift&>(*obj);
    muObj.sensRestore(this);
}


class MuPointwiseHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. err none */
    class Factory: public SensitivityFactory::ICreation {};
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MuPointwise, clazz);
        SUPERCLASS(BucketShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultMuPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(MuPointwise::NAME, 
                                    new Factory(), 
                                    new MuPointwise(),
                                    MuPointwise::IShift::TYPE);
    }

    static IObject* defaultMuPointwise(){
        return new MuPointwise();
    }
};

CClassConstSP const MuPointwise::TYPE = CClass::registerClassLoadMethod(
    "MuPointwise", typeid(MuPointwise), MuPointwiseHelper::load);

CClassConstSP const MuPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "MuPointwise::IShift", typeid(MuPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(MuPointwise::IRestorableShift, clazz);
    EXTENDS(MuPointwise::IShift);
}
    
CClassConstSP const MuPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "MuPointwise::IRestorableShift",
    typeid(MuPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE

