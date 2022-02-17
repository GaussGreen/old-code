//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaDDE.cpp
//
//   Description : Delta sensitivity for DDE
//
//   Author      : Qing Hou
//
//   Date        : 16 Mar 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DeltaDDE.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Results.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
DeltaDDE::Shift::~Shift(){} // empty
DeltaDDE::RestorableShift::~RestorableShift(){} // empty

const string DeltaDDE::NAME = "DELTA_DDE";
const string DeltaDDE::SECOND_ORDER_NAME = "GAMMA_DDE";
const double DeltaDDE::DEFAULT_SHIFT = 0.005;
const double DeltaDDE::MINIMUM_SHIFT = 0.0001;

/** constructor with explicit shift size */
DeltaDDE::DeltaDDE() : ScalarShift(TYPE, NAME) {};

DeltaDDE::DeltaDDE(double shiftSize) : ScalarShift(TYPE, NAME, shiftSize) {};

/** Scales delta and gamma numbers in results object */
void DeltaDDE::scaleResult(Results*     results,     // (M)
                        double       scaleFactor) const{
    // scale all results in packetName
    results->scale(outputName(), scaleFactor);
    results->scale(output2ndOrderName(), scaleFactor);
}


/** Adds delta and gamma numbers in results object */
void DeltaDDE::addResult(Results*           results,     // (M)
                      const Results*     resultsToAdd,
                      double             scaleFactor) const{
    // add all results in packetName
    results->add(outputName(), resultsToAdd, scaleFactor);
    results->add(output2ndOrderName(), resultsToAdd, scaleFactor);
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double DeltaDDE::divisor() const{
    static const char routine[] = "DeltaDDE::divisor";
    double initSpot;
    double shiftSize;
    try{
        // find initial value of spot - stored in sens control
        try {
            initSpot = getInitialValue();            
        }
        catch (ModelException& e){
            e.addMsg("Initial value of spot price has not been set in the " 
                     "delta calculation");
            throw e;
        }

        if (Maths::isZero(initSpot)) {
            throw ModelException(routine, "spot price is zero");
        }
        // then just scale by the shift size
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)){
            throw ModelException(routine, "Shift size is zero");
        }
    } catch (ModelException& e){
        e.addMsg(routine);
        throw e;
    }
    return (shiftSize * initSpot);
}

// own special calculate
void DeltaDDE::calculate(TweakGroup*      tweakGroup,
                         CResults*        results) {
    try {
        // can't use regular delta method as it conflicts with BasketDelta
        calculateTwoSidedDeriv(output2ndOrderName(), tweakGroup, results);
    } catch (exception& e){
        throw ModelException(e,  "DeltaDDE::calculate");
    }
}


/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP DeltaDDE::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP DeltaDDE::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool DeltaDDE::nameMatches(const OutputName&         name,
                        IObjectConstSP            obj){
    // cast obj to DeltaDDE::Shift and then invoke name method
    const Shift& deltaObj = 
        dynamic_cast<const Shift&>(*obj);
    return name.equals(deltaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void DeltaDDE::appendName(OutputNameArray&          namesList,
                       IObjectConstSP            obj){
    // cast obj to DeltaDDE::Shift and then invoke name method
    const Shift& deltaObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(deltaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool DeltaDDE::shift(IObjectSP obj) {
    // cast obj to DeltaDDE::Shift and then invoke shift method
    Shift& deltaObj = dynamic_cast<Shift&>(*obj);
    return deltaObj.sensShift(this);
}

    
void DeltaDDE::restore(IObjectSP obj) {
    // cast obj to DeltaDDE::Shift and then invoke restore method
    RestorableShift& deltaObj = 
        dynamic_cast<RestorableShift&>(*obj);
    deltaObj.sensRestore(this);
}

/** Returns the shift which has been made for the current pricing
    call */
ScalarShiftArray DeltaDDE::getComponentShifts() const{
    // need array of ScalarShiftConstSP to make this work - but this breaks
    // the array template - would need const array template...
    DeltaDDE* deltaShift = const_cast<DeltaDDE*>(this);
    return ScalarShiftArray(1, ScalarShiftSP::attachToRef(deltaShift));
}

class DeltaDDEHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new DeltaDDE(Delta::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new DeltaDDE(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaDDE, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultDeltaDDE);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(DeltaDDE::NAME, 
                                    new Factory(), 
                                    new DeltaDDE(Delta::DEFAULT_SHIFT),
                                    DeltaDDE::Shift::TYPE);
    }

    static IObject* defaultDeltaDDE(){
        return new DeltaDDE();
    }
};

CClassConstSP const DeltaDDE::TYPE = CClass::registerClassLoadMethod(
    "DeltaDDE", typeid(DeltaDDE), DeltaDDEHelper::load);

CClassConstSP const DeltaDDE::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "DeltaDDE::Shift", typeid(DeltaDDE::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(DeltaDDE::RestorableShift, clazz);
    EXTENDS(DeltaDDE::Shift);
}
    
CClassConstSP const DeltaDDE::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "DeltaDDE::RestorableShift", typeid(DeltaDDE::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
