//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaSkewProxyPointwise.cpp
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
#include "edginc/VegaSkewProxyPointwise.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
VegaSkewProxyPointwise::IShift::~IShift(){} // empty
VegaSkewProxyPointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega pointwise */
const string VegaSkewProxyPointwise::NAME = "VEGA_SKEW_PROXY_POINTWISE";
const double VegaSkewProxyPointwise::DEFAULT_SHIFT = VegaSkewPointwise::DEFAULT_SHIFT;

/** constructor with explicit shift size */
VegaSkewProxyPointwise::VegaSkewProxyPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize) {}

/** for reflection */
VegaSkewProxyPointwise::VegaSkewProxyPointwise(): 
    VectorShift(TYPE, NAME) {}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VegaSkewProxyPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return 1000.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VegaSkewProxyPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VegaSkewProxyPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VegaSkewProxyPointwise::nameMatches(const OutputName&         name,
                                         IObjectConstSP          obj){
    // cast obj to VegaSkewProxyPointwise::Shift and then invoke name method
    const IShift& vegaSkewObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaSkewObj.sensNameMatches(this, name);
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VegaSkewProxyPointwise::appendName(OutputNameArray&          namesList,
                                        IObjectConstSP          obj){
    // cast obj to VegaSkewProxyPointwise::Shift and then invoke name method
    const IShift& vegaSkewObj =
        dynamic_cast<const IShift&>(*obj);
    vegaSkewObj.sensAppendName(this, namesList);
}

bool VegaSkewProxyPointwise::shift(IObjectSP obj) {
    // cast obj to VegaSkewProxyPointwise::Shift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}
    
void VegaSkewProxyPointwise::restore(IObjectSP obj) {
    // cast obj to VegaSkewProxyPointwise::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the VegaSkewProxyPointwise.Shift interface */
IObjectConstSP VegaSkewProxyPointwise::qualifier(IObjectConstSP obj){
    // cast obj to VegaSkewProxyPointwise::Shift and then invoke qualifier method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaObj.sensExpiries(this);
}

class VegaSkewProxyPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VegaSkewProxyPointwise(VegaSkewProxyPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VegaSkewProxyPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VegaSkewProxyPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVegaSkewProxyPointwise);
        // register how to build our sensitivity
        SensitivityFactory::addSens(VegaSkewProxyPointwise::NAME, 
                                    new Factory(), 
                                    new VegaSkewProxyPointwise(VegaSkewProxyPointwise::DEFAULT_SHIFT),
                                    VegaSkewProxyPointwise::IShift::TYPE);
    }

    static IObject* defaultVegaSkewProxyPointwise(){
        return new VegaSkewProxyPointwise();
    }
};

CClassConstSP const VegaSkewProxyPointwise::TYPE = CClass::registerClassLoadMethod(
    "VegaSkewProxyPointwise", typeid(VegaSkewProxyPointwise),
    VegaSkewProxyPointwiseHelper::load);

CClassConstSP const VegaSkewProxyPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VegaSkewProxyPointwise::IShift", typeid(VegaSkewProxyPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VegaSkewProxyPointwise::IRestorableShift, clazz);
    EXTENDS(VegaSkewProxyPointwise::IShift);
}
    
CClassConstSP const VegaSkewProxyPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VegaSkewProxyPointwise::IRestorableShift",
    typeid(VegaSkewProxyPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
