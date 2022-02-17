//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RootTimeVegaProxy.hpp
//
//   Description : Controls calculation of proxy Root TIme Vega
//
//   Author      : Andrew J Swain
//
//   Date        : 13 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RootTimeVegaProxy.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
RootTimeVegaProxy::IShift::~IShift(){} // empty
RootTimeVegaProxy::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega parallel */
const string RootTimeVegaProxy::NAME = "RT_VEGA_PROXY";
const double RootTimeVegaProxy::DEFAULT_SHIFT = RootTimeVega::DEFAULT_SHIFT;

/** constructor with explicit shift size */
RootTimeVegaProxy::RootTimeVegaProxy(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

/** for reflection */
RootTimeVegaProxy::RootTimeVegaProxy(): 
    ScalarShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double RootTimeVegaProxy::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP RootTimeVegaProxy::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP RootTimeVegaProxy::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool RootTimeVegaProxy::nameMatches(const OutputName&         name,
                        IObjectConstSP          obj){
    // cast obj to RootTimeVegaProxy::Shift and then invoke name method
    const IShift& rootTimeVegaProxyObj =
        dynamic_cast<const IShift&>(*obj);
    return rootTimeVegaProxyObj.sensNameMatches(this, name);    
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void RootTimeVegaProxy::appendName(OutputNameArray&          namesList,
                       IObjectConstSP          obj){
    // cast obj to RootTimeVegaProxy::Shift and then invoke name method
    const IShift& rootTimeVegaProxyObj =
        dynamic_cast<const IShift&>(*obj);
    rootTimeVegaProxyObj.sensAppendName(this, namesList);
}

bool RootTimeVegaProxy::shift(IObjectSP obj) {
    // cast obj to RootTimeVegaProxy::IShift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

void RootTimeVegaProxy::restore(IObjectSP obj) {
    // cast obj to RootTimeVegaProxy::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

class RootTimeVegaProxyHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new RootTimeVegaProxy(RootTimeVegaProxy::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new RootTimeVegaProxy(shiftSize);
        }
    };
public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RootTimeVegaProxy, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultRootTimeVegaProxy);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(RootTimeVegaProxy::NAME, 
                                    new Factory(), 
                                    new RootTimeVegaProxy(RootTimeVegaProxy::DEFAULT_SHIFT),
                                    RootTimeVegaProxy::IShift::TYPE);
    }

    static IObject* defaultRootTimeVegaProxy(){
        return new RootTimeVegaProxy();
    }
};

CClassConstSP const RootTimeVegaProxy::TYPE = CClass::registerClassLoadMethod(
    "RootTimeVegaProxy", typeid(RootTimeVegaProxy), RootTimeVegaProxyHelper::load);

CClassConstSP const RootTimeVegaProxy::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "RootTimeVegaProxy::IShift", typeid(RootTimeVegaProxy::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(RootTimeVegaProxy::IRestorableShift, clazz);
    EXTENDS(RootTimeVegaProxy::IShift);
}
    
CClassConstSP const RootTimeVegaProxy::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "RootTimeVegaProxy::IRestorableShift", typeid(RootTimeVegaProxy::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
