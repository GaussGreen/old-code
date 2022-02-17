//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaSkewParallel.cpp
//
//   Description : Controls calculation of Vega Skew Parallel
//
//   Author      : Mark A Robson
//
//   Date        : 7 March 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE
VegaSkewParallel::IShift::IShift(){} // empty
VegaSkewParallel::IShift::~IShift(){} // empty
VegaSkewParallel::IRestorableShift::IRestorableShift(){} // empty
VegaSkewParallel::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega parallel */
const string VegaSkewParallel::NAME = "VEGA_SKEW_PARALLEL";
const double VegaSkewParallel::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
VegaSkewParallel::VegaSkewParallel(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize), spot(0.0), spotSet(false){
}

/** for reflection */
VegaSkewParallel::VegaSkewParallel(): 
    ScalarShift(TYPE, NAME), spot(0.0), spotSet(false) {
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VegaSkewParallel::divisor() const{
    // add check for 0 shift size - somewhere 
    return 1000.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VegaSkewParallel::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VegaSkewParallel::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VegaSkewParallel::nameMatches(const OutputName&         name,
                                   IObjectConstSP          obj){
    // cast obj to VegaSkewParallel::Shift and then invoke name 
    // matched method. Used to prod owners of vol surfaces to update
    // vega skew settings
    const IShift& vegaSkewObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaSkewObj.sensNameMatches(this, name);
}


bool VegaSkewParallel::IShift::sensNameMatches(VegaSkewParallel* shift, const OutputName& name) const{
    // Default implementation of IShift::sensNameMatches that matches with sensName
    return name.equals(sensName(shift));
}

void VegaSkewParallel::IShift::sensAppendName(VegaSkewParallel* shift, OutputNameArray& namesList) const{
    // Default implementation of IShift::sensNameMatches that adds only sensName
    OutputNameSP outputName(new OutputName(sensName(shift)));
    namesList.push_back(outputName);
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VegaSkewParallel::appendName(OutputNameArray&          namesList,
                                  IObjectConstSP          obj){
    // cast obj to VegaSkewParallel::Shift and then invoke name method
    const IShift& vegaSkewObj =
        dynamic_cast<const IShift&>(*obj);
    vegaSkewObj.sensAppendName(this, namesList);
}

bool VegaSkewParallel::shift(IObjectSP obj) {
    // cast obj to VegaSkewParallel::IShift and then invoke shift method
    IShift& vegaSkewObj =
        dynamic_cast<IShift&>(*obj);
    return vegaSkewObj.sensShift(this);
}

void VegaSkewParallel::restore(IObjectSP obj) {
    // cast obj to VegaSkewParallel::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

const double VegaSkewParallel::IShift::SKEW_NORMALISATION = 0.9;
const double VegaSkewParallel::IShift::LOG_SKEW_NORMALISATION = log(0.9);

/** Returns spot value of associated asset. Fails if spot has not been
    set */
double VegaSkewParallel::getSpot() const{
    if (!spotSet){
        throw ModelException("VegaSkewParallel::getSpot", 
                             "Spot level not set");
    }
    return spot;
}

/** Stores the spot value of the associated asset. Allows vol to 
    implement skew tweak */
void VegaSkewParallel::setSpot(double spotValue){
    spot = spotValue;
    spotSet = true;
}

/** return skew shift for given strike i.e.
    shift * log(strike/spot)/log(SKEW_NORMALISATION) */
double VegaSkewParallel::skewShift(double strike) const {
    return VegaSkewParallel::skewShift(getShiftSize(), strike, getSpot());
}

/** static version of skewShift so can reuse in other skew greeks */
// define skew shift to be logarithmic in region > 50% spot, then
// linear outside these boundaries
double VegaSkewParallel::skewShift(double shift, double strike, double spot) {
    static const double lowerLog = 0.5;
    static const double logHalf  = log(0.5);

    double skew;
    double loStrike = lowerLog*spot;

    if (strike < loStrike) {
        // lower linear region
        // slope is 1/K
        // skew = skew(50%) + (strike-loStrike)/loStrike
        skew = shift*(logHalf + (strike-loStrike)/loStrike)/
            VegaSkewParallel::IShift::LOG_SKEW_NORMALISATION;
    }
    else {
        // in log region
        skew = shift * log(strike/spot)/
            VegaSkewParallel::IShift::LOG_SKEW_NORMALISATION;
    }
    return skew;
}

class VegaSkewParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VegaSkewParallel(VegaSkewParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VegaSkewParallel(shiftSize);
        }
    };
public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VegaSkewParallel, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVegaSkewParallel);
        FIELD(spot, "spot");
        FIELD_MAKE_TRANSIENT(spot);
        FIELD(spotSet, "spotSet");
        FIELD_MAKE_TRANSIENT(spotSet);
        // register how to build our sensitivity
        SensitivityFactory::addSens(VegaSkewParallel::NAME, 
                                    new Factory(), 
                                    new VegaSkewParallel(VegaSkewParallel::DEFAULT_SHIFT),
                                    VegaSkewParallel::IShift::TYPE);
    }

    static IObject* defaultVegaSkewParallel(){
        return new VegaSkewParallel();
    }
};

CClassConstSP const VegaSkewParallel::TYPE = CClass::registerClassLoadMethod(
    "VegaSkewParallel", typeid(VegaSkewParallel), 
    VegaSkewParallelHelper::load);

CClassConstSP const VegaSkewParallel::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VegaSkewParallel::IShift", typeid(VegaSkewParallel::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VegaSkewParallel::IRestorableShift, clazz);
    EXTENDS(VegaSkewParallel::IShift);
}
    
CClassConstSP const VegaSkewParallel::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VegaSkewParallel::IRestorableShift", 
    typeid(VegaSkewParallel::IRestorableShift), 
    restorableShiftLoad);


class SkewShifter: public CObject {
    static CClassConstSP const TYPE;

    double shift;
    double strike;
    double spot;

    // compute skew shift for a given strike
    static double skew(SkewShifter* params){
        return VegaSkewParallel::skewShift(params->shift, 
                                           params->strike, 
                                           params->spot);
    }

    /** for reflection */
    SkewShifter():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(SkewShifter, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSkewShifter);
        FIELD(shift, "shift");
        FIELD(strike, "strike");
        FIELD(spot, "spot");

        Addin::registerClassDoubleMethod("SKEW_SHIFT",
                                         Addin::RISK,
                                         "compute skew shift for a given strike",
                                         TYPE,
                                         (Addin::DoubleMethod*)skew);
    }

    static IObject* defaultSkewShifter(){
        return new SkewShifter();
    }    
};

CClassConstSP const SkewShifter::TYPE = CClass::registerClassLoadMethod(
    "SkewShifter", typeid(SkewShifter), load);

DRLIB_END_NAMESPACE
