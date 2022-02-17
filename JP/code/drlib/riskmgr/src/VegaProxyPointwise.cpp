//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaProxyPointwise.cpp
//
//   Description : Controls calculation of fund proxy vega pointwise
//
//   Author      : Andrew J Swain
//
//   Date        : 12 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/VegaProxyPointwise.hpp"
#include "edginc/VegaPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
typedef smartPtr<VegaProxyPointwise::IRestorableShift> 
VegaProxyPointwiseRestorableShiftSP;

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
VegaProxyPointwise::IShift::~IShift(){} // empty
VegaProxyPointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for proxy vega pointwise */
const string VegaProxyPointwise::NAME = "VEGA_PROXY_POINTWISE";
const double VegaProxyPointwise::DEFAULT_SHIFT = VegaPointwise::DEFAULT_SHIFT;

/** constructor with explicit shift size */
VegaProxyPointwise::VegaProxyPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** constructor with explicit shift size and name override (override allows
    a VEGA_POINTWISE calculation to be stored under, eg, VEGA_MATRIX) */
VegaProxyPointwise::VegaProxyPointwise(const string& overrideName,
                                       double        shiftSize):
    VectorShift(TYPE, overrideName, shiftSize){
}

/** for reflection */
VegaProxyPointwise::VegaProxyPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VegaProxyPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VegaProxyPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VegaProxyPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VegaProxyPointwise::nameMatches(const OutputName&         name,
                                     IObjectConstSP          obj){
    // cast obj to VegaProxyPointwise::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaObj.sensNameMatches(this, name);
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VegaProxyPointwise::appendName(OutputNameArray&          namesList,
                                    IObjectConstSP          obj){
    // cast obj to VegaProxyPointwise::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    vegaObj.sensAppendName(this, namesList);
}

bool VegaProxyPointwise::shift(IObjectSP obj) {
    // cast obj to VegaProxyPointwise::Shift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

void VegaProxyPointwise::restore(IObjectSP obj) {
    // cast obj to VegaProxyPointwise::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the VegaProxyPointwise.Shift interface */
IObjectConstSP VegaProxyPointwise::qualifier(IObjectConstSP obj){
    // cast obj to VegaProxyPointwise::Shift and then invoke qualifier method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaObj.sensExpiries(this);
}

 
class VegaProxyPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VegaProxyPointwise(VegaProxyPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VegaProxyPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VegaProxyPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVegaProxyPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(VegaProxyPointwise::NAME, 
                                    new Factory(), 
                                    new VegaProxyPointwise(VegaProxyPointwise::DEFAULT_SHIFT),
                                    VegaProxyPointwise::IShift::TYPE);
    }

    static IObject* defaultVegaProxyPointwise(){
        return new VegaProxyPointwise();
    }
};

CClassConstSP const VegaProxyPointwise::TYPE = CClass::registerClassLoadMethod(
    "VegaProxyPointwise", typeid(VegaProxyPointwise), VegaProxyPointwiseHelper::load);

CClassConstSP const VegaProxyPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VegaProxyPointwise::IShift", typeid(VegaProxyPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VegaProxyPointwise::IRestorableShift, clazz);
    EXTENDS(VegaProxyPointwise::IShift);
}
    
CClassConstSP const VegaProxyPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VegaProxyPointwise::IRestorableShift",
    typeid(VegaProxyPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
