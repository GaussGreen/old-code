//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VHTConvexityFactorsPointwise.cpp
//
//   Description : Controls calculation of VHT convexity factors pointwise
//
//   Author      : Stewart Nesbitt
//
//   Date        : Sep 2005
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VHTConvexityFactorsPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
VHTConvexityFactorsPointwise::IShift::~IShift(){} // empty
VHTConvexityFactorsPointwise::IRestorableShift::~IRestorableShift(){} // empty

const string VHTConvexityFactorsPointwise::NAME = "VHT_CONVXT_FACTORS_POINTWISE";
const double VHTConvexityFactorsPointwise::DEFAULT_SHIFT = 0.01;

/** constructor with explicit shift size */
VHTConvexityFactorsPointwise::VHTConvexityFactorsPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** for reflection */
VHTConvexityFactorsPointwise::VHTConvexityFactorsPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VHTConvexityFactorsPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VHTConvexityFactorsPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VHTConvexityFactorsPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VHTConvexityFactorsPointwise::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to VHTConvexityFactorsPointwise::Shift and then invoke name method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(myObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VHTConvexityFactorsPointwise::appendName(OutputNameArray&          namesList,
                               IObjectConstSP          obj){
    // cast obj to VHTConvexityFactorsPointwise::Shift and then invoke name method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(myObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VHTConvexityFactorsPointwise::shift(IObjectSP obj) {
    // cast obj to VHTConvexityFactorsPointwise::Shift and then invoke shift method
    IShift& myObj =
        dynamic_cast<IShift&>(*obj);
    return myObj.sensShift(this);
}

void VHTConvexityFactorsPointwise::restore(IObjectSP obj) {
    // cast obj to VHTConvexityFactorsPointwise::Shift and then invoke restore method
    IRestorableShift& myObj =
        dynamic_cast<IRestorableShift&>(*obj);
    myObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the VHTConvexityFactorsPointwise.Shift interface */
IObjectConstSP VHTConvexityFactorsPointwise::qualifier(IObjectConstSP obj){
    // cast obj to VHTConvexityFactorsPointwise::Shift and then invoke qualifier method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    return myObj.sensExpiries(this);
}

class VHTConvexityFactorsPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VHTConvexityFactorsPointwise(VHTConvexityFactorsPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VHTConvexityFactorsPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VHTConvexityFactorsPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVHTConvexityFactorsPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(VHTConvexityFactorsPointwise::NAME, 
                                    new Factory(), 
                                    new VHTConvexityFactorsPointwise(VHTConvexityFactorsPointwise::DEFAULT_SHIFT),
                                    VHTConvexityFactorsPointwise::IShift::TYPE);
    }

    static IObject* defaultVHTConvexityFactorsPointwise(){
        return new VHTConvexityFactorsPointwise();
    }
};

CClassConstSP const VHTConvexityFactorsPointwise::TYPE = CClass::registerClassLoadMethod(
    "VHTConvexityFactorsPointwise", typeid(VHTConvexityFactorsPointwise), VHTConvexityFactorsPointwiseHelper::load);

CClassConstSP const VHTConvexityFactorsPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VHTConvexityFactorsPointwise::IShift", typeid(VHTConvexityFactorsPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VHTConvexityFactorsPointwise::IRestorableShift, clazz);
    EXTENDS(VHTConvexityFactorsPointwise::IShift);
}
    
CClassConstSP const VHTConvexityFactorsPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VHTConvexityFactorsPointwise::IRestorableShift",
    typeid(VHTConvexityFactorsPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
