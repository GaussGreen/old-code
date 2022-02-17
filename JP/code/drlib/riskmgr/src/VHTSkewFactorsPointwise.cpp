//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VHTSkewFactorsPointwise.cpp
//
//   Description : Controls calculation of VHT skew factors pointwise
//
//   Author      : Stewart Nesbitt
//
//   Date        : Sep 2005
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VHTSkewFactorsPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
VHTSkewFactorsPointwise::IShift::~IShift(){} // empty
VHTSkewFactorsPointwise::IRestorableShift::~IRestorableShift(){} // empty

const string VHTSkewFactorsPointwise::NAME = "VHT_SKEW_FACTORS_POINTWISE";
const double VHTSkewFactorsPointwise::DEFAULT_SHIFT = 0.01;

/** constructor with explicit shift size */
VHTSkewFactorsPointwise::VHTSkewFactorsPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** for reflection */
VHTSkewFactorsPointwise::VHTSkewFactorsPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VHTSkewFactorsPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VHTSkewFactorsPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VHTSkewFactorsPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VHTSkewFactorsPointwise::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to VHTSkewFactorsPointwise::Shift and then invoke name method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(myObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VHTSkewFactorsPointwise::appendName(OutputNameArray&          namesList,
                               IObjectConstSP          obj){
    // cast obj to VHTSkewFactorsPointwise::Shift and then invoke name method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(myObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VHTSkewFactorsPointwise::shift(IObjectSP obj) {
    // cast obj to VHTSkewFactorsPointwise::Shift and then invoke shift method
    IShift& myObj =
        dynamic_cast<IShift&>(*obj);
    return myObj.sensShift(this);
}

void VHTSkewFactorsPointwise::restore(IObjectSP obj) {
    // cast obj to VHTSkewFactorsPointwise::Shift and then invoke restore method
    IRestorableShift& myObj =
        dynamic_cast<IRestorableShift&>(*obj);
    myObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the VHTSkewFactorsPointwise.Shift interface */
IObjectConstSP VHTSkewFactorsPointwise::qualifier(IObjectConstSP obj){
    // cast obj to VHTSkewFactorsPointwise::Shift and then invoke qualifier method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    return myObj.sensExpiries(this);
}

class VHTSkewFactorsPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VHTSkewFactorsPointwise(VHTSkewFactorsPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VHTSkewFactorsPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VHTSkewFactorsPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVHTSkewFactorsPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(VHTSkewFactorsPointwise::NAME, 
                                    new Factory(), 
                                    new VHTSkewFactorsPointwise(VHTSkewFactorsPointwise::DEFAULT_SHIFT),
                                    VHTSkewFactorsPointwise::IShift::TYPE);
    }

    static IObject* defaultVHTSkewFactorsPointwise(){
        return new VHTSkewFactorsPointwise();
    }
};

CClassConstSP const VHTSkewFactorsPointwise::TYPE = CClass::registerClassLoadMethod(
    "VHTSkewFactorsPointwise", typeid(VHTSkewFactorsPointwise), VHTSkewFactorsPointwiseHelper::load);

CClassConstSP const VHTSkewFactorsPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VHTSkewFactorsPointwise::IShift", typeid(VHTSkewFactorsPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VHTSkewFactorsPointwise::IRestorableShift, clazz);
    EXTENDS(VHTSkewFactorsPointwise::IShift);
}
    
CClassConstSP const VHTSkewFactorsPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VHTSkewFactorsPointwise::IRestorableShift",
    typeid(VHTSkewFactorsPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
