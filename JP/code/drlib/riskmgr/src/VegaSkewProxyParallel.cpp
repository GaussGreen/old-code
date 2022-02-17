//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaSkewProxyParallel.cpp
//
//   Description : Controls calculation of fund proxy vega skew
//
//   Author      : Andrew J Swain
//
//   Date        : 12 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VegaSkewProxyParallel.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
VegaSkewProxyParallel::IShift::~IShift(){} // empty
VegaSkewProxyParallel::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega parallel */
const string VegaSkewProxyParallel::NAME = "VEGA_SKEW_PROXY_PARALLEL";
const double VegaSkewProxyParallel::DEFAULT_SHIFT = VegaSkewParallel::DEFAULT_SHIFT;

/** constructor with explicit shift size */
VegaSkewProxyParallel::VegaSkewProxyParallel(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){}

/** for reflection */
VegaSkewProxyParallel::VegaSkewProxyParallel(): 
    ScalarShift(TYPE, NAME) {}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VegaSkewProxyParallel::divisor() const{
    // add check for 0 shift size - somewhere 
    return 1000.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VegaSkewProxyParallel::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VegaSkewProxyParallel::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VegaSkewProxyParallel::nameMatches(const OutputName&    name,
                                   IObjectConstSP          obj){
    // cast obj to VegaSkewProxyParallel::Shift and then invoke name method
    const IShift& vegaSkewObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaSkewObj.sensNameMatches(this, name);
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VegaSkewProxyParallel::appendName(OutputNameArray&     namesList,
                                  IObjectConstSP          obj){
    // cast obj to VegaSkewProxyParallel::Shift and then invoke name method
    const IShift& vegaSkewObj =
        dynamic_cast<const IShift&>(*obj);
    vegaSkewObj.sensAppendName(this, namesList);
}

bool VegaSkewProxyParallel::shift(IObjectSP obj) {
    // cast obj to VegaSkewProxyParallel::IShift and then invoke shift method
    IShift& vegaSkewObj =
        dynamic_cast<IShift&>(*obj);
    return vegaSkewObj.sensShift(this);
}

void VegaSkewProxyParallel::restore(IObjectSP obj) {
    // cast obj to VegaSkewProxyParallel::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

class VegaSkewProxyParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VegaSkewProxyParallel(VegaSkewProxyParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VegaSkewProxyParallel(shiftSize);
        }
    };
public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VegaSkewProxyParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVegaSkewProxyParallel);
        // register how to build our sensitivity
        SensitivityFactory::addSens(VegaSkewProxyParallel::NAME, 
                                    new Factory(), 
                                    new VegaSkewProxyParallel(VegaSkewProxyParallel::DEFAULT_SHIFT),
                                    VegaSkewProxyParallel::IShift::TYPE);
    }

    static IObject* defaultVegaSkewProxyParallel(){
        return new VegaSkewProxyParallel();
    }
};

CClassConstSP const VegaSkewProxyParallel::TYPE = CClass::registerClassLoadMethod(
    "VegaSkewProxyParallel", typeid(VegaSkewProxyParallel), 
    VegaSkewProxyParallelHelper::load);

CClassConstSP const VegaSkewProxyParallel::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VegaSkewProxyParallel::IShift", typeid(VegaSkewProxyParallel::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VegaSkewProxyParallel::IRestorableShift, clazz);
    EXTENDS(VegaSkewProxyParallel::IShift);
}
    
CClassConstSP const VegaSkewProxyParallel::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VegaSkewProxyParallel::IRestorableShift", 
    typeid(VegaSkewProxyParallel::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
