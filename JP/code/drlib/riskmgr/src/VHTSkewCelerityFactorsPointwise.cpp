//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VHTSkewCelerityFactorsPointwise.cpp
//
//   Description : Controls calculation of VHT skew celerity factors pointwise
//
//   Author      : Stewart Nesbitt
//
//   Date        : Sep 2005
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VHTSkewCelerityFactorsPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
VHTSkewCelerityFactorsPointwise::IShift::~IShift(){} // empty
VHTSkewCelerityFactorsPointwise::IRestorableShift::~IRestorableShift(){} // empty

// Limit of 40 chars for sens name...

const string VHTSkewCelerityFactorsPointwise::NAME = "VHT_SKEW_CELERITY_FACTORS_POINTWISE";
const double VHTSkewCelerityFactorsPointwise::DEFAULT_SHIFT = 0.01;

/** constructor with explicit shift size */
VHTSkewCelerityFactorsPointwise::VHTSkewCelerityFactorsPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** for reflection */
VHTSkewCelerityFactorsPointwise::VHTSkewCelerityFactorsPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VHTSkewCelerityFactorsPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VHTSkewCelerityFactorsPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VHTSkewCelerityFactorsPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VHTSkewCelerityFactorsPointwise::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to VHTSkewCelerityFactorsPointwise::Shift and then invoke name method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(myObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VHTSkewCelerityFactorsPointwise::appendName(OutputNameArray&          namesList,
                               IObjectConstSP          obj){
    // cast obj to VHTSkewCelerityFactorsPointwise::Shift and then invoke name method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(myObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VHTSkewCelerityFactorsPointwise::shift(IObjectSP obj) {
    // cast obj to VHTSkewCelerityFactorsPointwise::Shift and then invoke shift method
    IShift& myObj =
        dynamic_cast<IShift&>(*obj);
    return myObj.sensShift(this);
}

void VHTSkewCelerityFactorsPointwise::restore(IObjectSP obj) {
    // cast obj to VHTSkewCelerityFactorsPointwise::Shift and then invoke restore method
    IRestorableShift& myObj =
        dynamic_cast<IRestorableShift&>(*obj);
    myObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the VHTSkewCelerityFactorsPointwise.Shift interface */
IObjectConstSP VHTSkewCelerityFactorsPointwise::qualifier(IObjectConstSP obj){
    // cast obj to VHTSkewCelerityFactorsPointwise::Shift and then invoke qualifier method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    return myObj.sensExpiries(this);
}

class VHTSkewCelerityFactorsPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VHTSkewCelerityFactorsPointwise(VHTSkewCelerityFactorsPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VHTSkewCelerityFactorsPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VHTSkewCelerityFactorsPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVHTSkewCelerityFactorsPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(VHTSkewCelerityFactorsPointwise::NAME, 
                                    new Factory(), 
                                    new VHTSkewCelerityFactorsPointwise(VHTSkewCelerityFactorsPointwise::DEFAULT_SHIFT),
                                    VHTSkewCelerityFactorsPointwise::IShift::TYPE);
    }

    static IObject* defaultVHTSkewCelerityFactorsPointwise(){
        return new VHTSkewCelerityFactorsPointwise();
    }
};

CClassConstSP const VHTSkewCelerityFactorsPointwise::TYPE = CClass::registerClassLoadMethod(
    "VHTSkewCelerityFactorsPointwise", typeid(VHTSkewCelerityFactorsPointwise), VHTSkewCelerityFactorsPointwiseHelper::load);

CClassConstSP const VHTSkewCelerityFactorsPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VHTSkewCelerityFactorsPointwise::IShift", typeid(VHTSkewCelerityFactorsPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VHTSkewCelerityFactorsPointwise::IRestorableShift, clazz);
    EXTENDS(VHTSkewCelerityFactorsPointwise::IShift);
}
    
CClassConstSP const VHTSkewCelerityFactorsPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VHTSkewCelerityFactorsPointwise::IRestorableShift",
    typeid(VHTSkewCelerityFactorsPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
