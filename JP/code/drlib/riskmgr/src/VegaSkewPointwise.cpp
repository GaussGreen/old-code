//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaSkewPointwise.cpp
//
//   Description : Controls calculation of Vega skew pointwise
//
//   Author      : Mark A Robson
//
//   Date        : 5 March 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/VegaSkewProxyPointwise.hpp"

DRLIB_BEGIN_NAMESPACE
VegaSkewPointwise::IShift::IShift(){} // empty
VegaSkewPointwise::IShift::~IShift(){} // empty
VegaSkewPointwise::IRestorableShift::IRestorableShift(){} // empty
VegaSkewPointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega pointwise */
const string VegaSkewPointwise::NAME = "VEGA_SKEW_POINTWISE";
const double VegaSkewPointwise::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
VegaSkewPointwise::VegaSkewPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize), spot(0.0), spotSet(false) {
}

/** for reflection */
VegaSkewPointwise::VegaSkewPointwise(): 
    VectorShift(TYPE, NAME), spot(0.0), spotSet(false) {
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VegaSkewPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return 1000.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VegaSkewPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VegaSkewPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VegaSkewPointwise::nameMatches(const OutputName&         name,
                                   IObjectConstSP          obj){
    // cast obj to VegaSkewPointwise::Shift and then invoke name 
    // matched method. Used to prod owners of vol surfaces to update
    // vega skew settings
    const IShift& vegaSkewObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaSkewObj.sensNameMatches(this, name);
}

bool VegaSkewPointwise::IShift::sensNameMatches(VegaSkewPointwise* shift, const OutputName& name) const{
    // Default implementation of IShift::sensNameMatches that matches with sensName
    return name.equals(sensName(shift));
}

void VegaSkewPointwise::IShift::sensAppendName(VegaSkewPointwise* shift, OutputNameArray& namesList) const{
    // Default implementation of IShift::sensNameMatches that adds only sensName
    OutputNameSP outputName(new OutputName(sensName(shift)));
    namesList.push_back(outputName);
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VegaSkewPointwise::appendName(OutputNameArray&          namesList,
                                  IObjectConstSP          obj){
    // cast obj to VegaSkewPointwise::Shift and then invoke name method
    const IShift& vegaSkewObj =
        dynamic_cast<const IShift&>(*obj);
    vegaSkewObj.sensAppendName(this, namesList);
}

bool VegaSkewPointwise::shift(IObjectSP obj) {
    // cast obj to VegaSkewPointwise::Shift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

void VegaSkewPointwise::restore(IObjectSP obj) {
    // cast obj to VegaSkewPointwise::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the VegaSkewPointwise.Shift interface */
IObjectConstSP VegaSkewPointwise::qualifier(IObjectConstSP obj){
    // cast obj to VegaSkewPointwise::Shift and then invoke qualifier method
    const IShift& vegaSkewObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaSkewObj.sensExpiries(this);
}

/** Returns spot value of associated asset. Fails if spot has not been
    set */
double VegaSkewPointwise::getSpot() const{
    if (!spotSet){
        throw ModelException("VegaSkewPointwise::getSpot", 
                             "Spot level not set");
    }
    return spot;
}

/** Stores the spot value of the associated asset. Allows vol to 
    implement skew tweak */
void VegaSkewPointwise::setSpot(double spotValue){
    spot = spotValue;
    spotSet = true;
}

/** return skew shift for given strike i.e.
    -shift * log(strike/spot)/log(SKEW_NORMALISATION) */
double VegaSkewPointwise::skewShift(double strike) const {
    return VegaSkewParallel::skewShift(getShiftSize(), strike, getSpot()); 
}

VegaSkewPointwise* VegaSkewPointwise::fromProxy(VegaSkewProxyPointwise* proxy) {
    VegaSkewPointwiseSP vs(new VegaSkewPointwise(proxy->getShiftSize()));
    vs->expiry  = proxy->getExpiry();
    vs->setMarketDataName(proxy->getMarketDataName());   
    return vs.release();
}
 
class VegaSkewPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VegaSkewPointwise(VegaSkewPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VegaSkewPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VegaSkewPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVegaSkewPointwise);
        FIELD(spot, "spot");
        FIELD_MAKE_TRANSIENT(spot);
        FIELD(spotSet, "spotSet");
        FIELD_MAKE_TRANSIENT(spotSet);
        // register how to build our sensitivity
        SensitivityFactory::addSens(VegaSkewPointwise::NAME, 
                                    new Factory(), 
                                    new VegaSkewPointwise(VegaSkewPointwise::DEFAULT_SHIFT),
                                    VegaSkewPointwise::IShift::TYPE);
    }

    static IObject* defaultVegaSkewPointwise(){
        return new VegaSkewPointwise();
    }
};

CClassConstSP const VegaSkewPointwise::TYPE = CClass::registerClassLoadMethod(
    "VegaSkewPointwise", typeid(VegaSkewPointwise),
    VegaSkewPointwiseHelper::load);

CClassConstSP const VegaSkewPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VegaSkewPointwise::IShift", typeid(VegaSkewPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VegaSkewPointwise::IRestorableShift, clazz);
    EXTENDS(VegaSkewPointwise::IShift);
}
    
CClassConstSP const VegaSkewPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VegaSkewPointwise::IRestorableShift",
    typeid(VegaSkewPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
