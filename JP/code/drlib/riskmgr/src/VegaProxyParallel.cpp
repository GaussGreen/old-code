//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaProxyParallel.cpp
//
//   Description : Fund proxy vega sensitivity
//
//   Author      : Andrew J Swain
//
//   Date        : 8 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/VegaProxyParallel.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
VegaProxyParallel::Shift::~Shift(){} // empty
VegaProxyParallel::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for VegaProxy parallel */
// use same name as regular vega
const string VegaProxyParallel::NAME = "VEGA_PROXY_PARALLEL";
const double VegaProxyParallel::DEFAULT_SHIFT = VegaParallel::DEFAULT_SHIFT;

/** constructor with explicit shift size */
VegaProxyParallel::VegaProxyParallel(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
VegaProxyParallel::VegaProxyParallel(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VegaProxyParallel::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VegaProxyParallel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VegaProxyParallel::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    VegaProxyParallel.Shift interface */
bool VegaProxyParallel::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to VegaProxyParallel::Shift and then invoke name method
    const Shift& vegaObj =
        dynamic_cast<const Shift&>(*obj);
    return vegaObj.sensNameMatches(this, name);
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void VegaProxyParallel::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to VegaProxyParallel::Shift and then invoke name method
    const Shift& vegaObj =
        dynamic_cast<const Shift&>(*obj);
    vegaObj.sensAppendName(this, namesList);
}

bool VegaProxyParallel::shift(IObjectSP obj) {
    // cast obj to VegaProxyParallel::Shift and then invoke shift method
    Shift& vegaObj =
        dynamic_cast<Shift&>(*obj);
    return vegaObj.sensShift(this);
}

void VegaProxyParallel::restore(IObjectSP obj) {
    // cast obj to VegaProxyParallel::Shift and then invoke restore method
    RestorableShift& vegaObj =
        dynamic_cast<RestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

class VegaProxyParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VegaProxyParallel(VegaProxyParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VegaProxyParallel(shiftSize);
        }
    };
        
public:
    /** Invoked when this class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VegaProxyParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVegaProxyParallel);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(VegaProxyParallel::NAME, 
                                    new Factory(), 
                                    new VegaProxyParallel(VegaProxyParallel::DEFAULT_SHIFT),
                                    VegaProxyParallel::Shift::TYPE);
    }

    static IObject* defaultVegaProxyParallel(){
        return new VegaProxyParallel();
    }
};

CClassConstSP const VegaProxyParallel::TYPE = CClass::registerClassLoadMethod(
    "VegaProxyParallel", typeid(VegaProxyParallel), VegaProxyParallelHelper::load);

CClassConstSP const VegaProxyParallel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VegaProxyParallel::Shift", typeid(VegaProxyParallel::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VegaProxyParallel::RestorableShift, clazz);
    EXTENDS(VegaProxyParallel::Shift);
}
    
CClassConstSP const VegaProxyParallel::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VegaProxyParallel::RestorableShift", typeid(VegaProxyParallel::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
