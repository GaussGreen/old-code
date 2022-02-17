//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VHTConvexityCelerityFactorsPointwise.cpp
//
//   Description : Controls calculation of VHT convexity celerity factors pointwise
//
//   Author      : Stewart Nesbitt
//
//   Date        : Sep 2005
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VHTConvexityCelerityFactorsPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
VHTConvexityCelerityFactorsPointwise::IShift::~IShift(){} // empty
VHTConvexityCelerityFactorsPointwise::IRestorableShift::~IRestorableShift(){} // empty

const string VHTConvexityCelerityFactorsPointwise::NAME = "VHT_CONVXT_CELERITY_FACTORS_POINTWISE";
const double VHTConvexityCelerityFactorsPointwise::DEFAULT_SHIFT = 0.01;

/** constructor with explicit shift size */
VHTConvexityCelerityFactorsPointwise::VHTConvexityCelerityFactorsPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** for reflection */
VHTConvexityCelerityFactorsPointwise::VHTConvexityCelerityFactorsPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VHTConvexityCelerityFactorsPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VHTConvexityCelerityFactorsPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VHTConvexityCelerityFactorsPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VHTConvexityCelerityFactorsPointwise::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to VHTConvexityCelerityFactorsPointwise::Shift and then invoke name method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(myObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VHTConvexityCelerityFactorsPointwise::appendName(OutputNameArray&          namesList,
                               IObjectConstSP          obj){
    // cast obj to VHTConvexityCelerityFactorsPointwise::Shift and then invoke name method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(myObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VHTConvexityCelerityFactorsPointwise::shift(IObjectSP obj) {
    // cast obj to VHTConvexityCelerityFactorsPointwise::Shift and then invoke shift method
    IShift& myObj =
        dynamic_cast<IShift&>(*obj);
    return myObj.sensShift(this);
}

void VHTConvexityCelerityFactorsPointwise::restore(IObjectSP obj) {
    // cast obj to VHTConvexityCelerityFactorsPointwise::Shift and then invoke restore method
    IRestorableShift& myObj =
        dynamic_cast<IRestorableShift&>(*obj);
    myObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the VHTConvexityCelerityFactorsPointwise.Shift interface */
IObjectConstSP VHTConvexityCelerityFactorsPointwise::qualifier(IObjectConstSP obj){
    // cast obj to VHTConvexityCelerityFactorsPointwise::Shift and then invoke qualifier method
    const IShift& myObj =
        dynamic_cast<const IShift&>(*obj);
    return myObj.sensExpiries(this);
}

class VHTConvexityCelerityFactorsPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VHTConvexityCelerityFactorsPointwise(VHTConvexityCelerityFactorsPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VHTConvexityCelerityFactorsPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VHTConvexityCelerityFactorsPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVHTConvexityCelerityFactorsPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(VHTConvexityCelerityFactorsPointwise::NAME, 
                                    new Factory(), 
                                    new VHTConvexityCelerityFactorsPointwise(VHTConvexityCelerityFactorsPointwise::DEFAULT_SHIFT),
                                    VHTConvexityCelerityFactorsPointwise::IShift::TYPE);
    }

    static IObject* defaultVHTConvexityCelerityFactorsPointwise(){
        return new VHTConvexityCelerityFactorsPointwise();
    }
};

CClassConstSP const VHTConvexityCelerityFactorsPointwise::TYPE = CClass::registerClassLoadMethod(
    "VHTConvexityCelerityFactorsPointwise", typeid(VHTConvexityCelerityFactorsPointwise), VHTConvexityCelerityFactorsPointwiseHelper::load);

CClassConstSP const VHTConvexityCelerityFactorsPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VHTConvexityCelerityFactorsPointwise::IShift", typeid(VHTConvexityCelerityFactorsPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VHTConvexityCelerityFactorsPointwise::IRestorableShift, clazz);
    EXTENDS(VHTConvexityCelerityFactorsPointwise::IShift);
}
    
CClassConstSP const VHTConvexityCelerityFactorsPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VHTConvexityCelerityFactorsPointwise::IRestorableShift",
    typeid(VHTConvexityCelerityFactorsPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
