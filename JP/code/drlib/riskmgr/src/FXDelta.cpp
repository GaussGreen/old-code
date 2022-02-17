//----------------------------------------------------------------------------
//
//   Group       : EDR
//
//   Filename    : FXDelta.hpp
//
//   Description : FX Delta sensitivity
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/FXDelta.hpp"
#include "edginc/Results.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
FXDelta::Shift::~Shift(){} // empty
FXDelta::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for vega parallel */
const string FXDelta::NAME = "FX_DELTA";
const string FXDelta::SECOND_ORDER_NAME = "FX_GAMMA";
const double FXDelta::DEFAULT_SHIFT = 0.005;

/** constructor with explicit shift size */
FXDelta::FXDelta(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize){
}

FXDelta::FXDelta(double     shiftSize,
                 IModel*    model,
                 Control*   control):
ScalarShift(TYPE, NAME, shiftSize) {
    this->algorithm = model;
    this->control   = control;
}

/** for reflection */
FXDelta::FXDelta(): 
    ScalarShift(TYPE, NAME){
}

/** Scales delta and gamma numbers in results object */
void FXDelta::scaleResult(Results*     results,     // (M)
                          double       scaleFactor) const{
    // scale all results in packetName
    results->scale(NAME, scaleFactor);
    results->scale(SECOND_ORDER_NAME, scaleFactor);
}


/** Adds delta and gamma numbers in results object */
void FXDelta::addResult(Results*           results,     // (M)
                        const Results*     resultsToAdd,
                        double             scaleFactor) const{
    // add all results in packetName
    results->add(NAME, resultsToAdd, scaleFactor);
    results->add(SECOND_ORDER_NAME, resultsToAdd, scaleFactor);
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double FXDelta::divisor() const{
    static const char routine[] = "FXDelta::divisor";
    double initSpot;
    double shiftSize;
    try{
        // find initial value of spot - stored in sens control
        try{
            initSpot = getInitialValue();
        } catch (ModelException& e){
            e.addMsg("Initial value of spot price has not been set in the " 
                     "fx delta calculation");
            throw e;
        }
        // then just scale by the shift size
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)){
            throw ModelException(routine, "FX delta shift size is zero");
        }
    } catch (ModelException& e){
        e.addMsg(routine);
        throw e;
    }
    return (shiftSize * initSpot);
}

/** override default one sided shift with a two sided one */
void FXDelta::calculate(TweakGroup*  tweakGroup,
                        CResults*    results){
    try{
        calculateTwoSidedDeriv(SECOND_ORDER_NAME, tweakGroup, results);
    } catch (exception& e){
        throw ModelException(e,  "FXDelta::calculate");
    }
}


/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP FXDelta::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP FXDelta::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool FXDelta::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to FXDelta::Shift and then invoke name method
    const Shift& fxDeltaObj =
        dynamic_cast<const Shift&>(*obj);
    return name.equals(fxDeltaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void FXDelta::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to FXDelta::Shift and then invoke name method
    const Shift& fxDeltaObj =
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(fxDeltaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool FXDelta::shift(IObjectSP obj) {
    // cast obj to Delta::Shift and then invoke shift method
    Shift& fxDeltaObj =
        dynamic_cast<Shift&>(*obj);
    return fxDeltaObj.sensShift(this);
}

void FXDelta::restore(IObjectSP obj) {
    // cast obj to Delta::Shift and then invoke restore method
    RestorableShift& fxDeltaObj =
        dynamic_cast<RestorableShift&>(*obj);
    fxDeltaObj.sensRestore(this);
}

/** Returns the shift which has been made for the current pricing
    call */
ScalarShiftArray FXDelta::getComponentShifts() const{
    // need array of ScalarShiftConstSP to make this work - but this breaks
    // the array template - would need const array template...
    FXDelta* deltaShift = const_cast<FXDelta*>(this);
    return ScalarShiftArray(1, ScalarShiftSP::attachToRef(deltaShift));
}

class FXDeltaHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new FXDelta(FXDelta::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new FXDelta(shiftSize);
        }
    };
       
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FXDelta, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultFXDelta);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(FXDelta::NAME, 
                                    new Factory(), 
                                    new FXDelta(FXDelta::DEFAULT_SHIFT),
                                    FXDelta::Shift::TYPE);
    }

    static IObject* defaultFXDelta(){
        return new FXDelta();
    }
};

CClassConstSP const FXDelta::TYPE = CClass::registerClassLoadMethod(
    "FXDelta", typeid(FXDelta), FXDeltaHelper::load);

CClassConstSP const FXDelta::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "FXDelta::Shift", typeid(FXDelta::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(FXDelta::RestorableShift, clazz);
    EXTENDS(FXDelta::Shift);
}
    
CClassConstSP const FXDelta::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "FXDelta::RestorableShift", typeid(FXDelta::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
