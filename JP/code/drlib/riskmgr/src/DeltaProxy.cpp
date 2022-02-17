//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaProxy.cpp
//
//   Description : Fund proxy delta sensitivity
//
//   Author      : Andrew J Swain
//
//   Date        : 7 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DeltaProxy.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Results.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
DeltaProxy::Shift::~Shift(){} // empty
DeltaProxy::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for delta proxy */
const string DeltaProxy::NAME = "DELTA_PROXY";
const string DeltaProxy::SECOND_ORDER_NAME = "GAMMA_PROXY";
const double DeltaProxy::DEFAULT_SHIFT = Delta::DEFAULT_SHIFT;
const double DeltaProxy::MINIMUM_SHIFT = Delta::MINIMUM_SHIFT;

/** for reflection */
DeltaProxy::DeltaProxy(): ScalarShift(TYPE, NAME) {}

/** Scales delta and gamma numbers in results object */
void DeltaProxy::scaleResult(Results*     results,     // (M)
                             double       scaleFactor) const{
    // scale all results in packetName
    results->scale(NAME, scaleFactor);
    results->scale(SECOND_ORDER_NAME, scaleFactor);
}


/** Adds delta and gamma numbers in results object */
void DeltaProxy::addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const{
    // add all results in packetName
    results->add(NAME, resultsToAdd, scaleFactor);
    results->add(SECOND_ORDER_NAME, resultsToAdd, scaleFactor);
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double DeltaProxy::divisor() const{
    static const char routine[] = "DeltaProxy::divisor";
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

/** override default one sided shift with a two sided one */
void DeltaProxy::calculate(TweakGroup*      tweakGroup,
                           CResults*        results){
    try {
        calculateTwoSidedDeriv(SECOND_ORDER_NAME, tweakGroup, results);
    } catch (exception& e){
        throw ModelException(e, "DeltaProxy::calculate");
    }
}


/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP DeltaProxy::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP DeltaProxy::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool DeltaProxy::nameMatches(const OutputName&     name,
                             IObjectConstSP      obj) {
    // cast obj to DeltaProxy::Shift and then invoke name method
    const Shift& deltaObj = 
        dynamic_cast<const Shift&>(*obj);
    return deltaObj.sensNameMatches(this, name);
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void DeltaProxy::appendName(OutputNameArray&      namesList,
                            IObjectConstSP      obj){
    // cast obj to DeltaProxy::Shift and then invoke name method
    const Shift& deltaObj = 
        dynamic_cast<const Shift&>(*obj);
    deltaObj.sensAppendName(this, namesList);
}

bool DeltaProxy::shift(IObjectSP obj) {
    // cast obj to DeltaProxy::Shift and then invoke shift method
    Shift& deltaObj = 
        dynamic_cast<Shift&>(*obj);
    return deltaObj.sensShift(this);
}

void DeltaProxy::restore(IObjectSP obj) {
    // cast obj to DeltaProxy::Shift and then invoke restore method
    RestorableShift& deltaObj = 
        dynamic_cast<RestorableShift&>(*obj);
    deltaObj.sensRestore(this);
}

/** Returns the shift which has been made for the current pricing
    call */
ScalarShiftArray DeltaProxy::getComponentShifts() const{
    // need array of ScalarShiftConstSP to make this work - but this breaks
    // the array template - would need const array template...
    DeltaProxy* deltaShift = const_cast<DeltaProxy*>(this);
    return ScalarShiftArray(1, ScalarShiftSP::attachToRef(deltaShift));
}

/** constructor with explicit shift size */
DeltaProxy::DeltaProxy(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize) {
}
  
/** make a delta shift using own shift */
DeltaSP DeltaProxy::makeDeltaShift() const {
    static const string method = "DeltaProxy::makeDeltaShift";
    try {
        DeltaSP delta(new Delta(getShiftSize()));

        if (hasOverrideNames()) {
            delta->storeOverrideNames(toTweak);
        }
        delta->setMarketDataName(getMarketDataName());
        return delta;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

class DeltaProxyHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new DeltaProxy(DeltaProxy::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new DeltaProxy(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaProxy, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultDeltaProxy);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(DeltaProxy::NAME, 
                                    new Factory(), 
                                    new DeltaProxy(DeltaProxy::DEFAULT_SHIFT),
                                    DeltaProxy::Shift::TYPE);
    }

    static IObject* defaultDeltaProxy(){
        return new DeltaProxy();
    }
};

CClassConstSP const DeltaProxy::TYPE = CClass::registerClassLoadMethod(
    "DeltaProxy", typeid(DeltaProxy), DeltaProxyHelper::load);

CClassConstSP const DeltaProxy::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "DeltaProxy::Shift", typeid(DeltaProxy::Shift), 0);

static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(DeltaProxy::RestorableShift, clazz);
    EXTENDS(DeltaProxy::Shift);
}
    
CClassConstSP const DeltaProxy::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "DeltaProxy::RestorableShift", typeid(DeltaProxy::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
