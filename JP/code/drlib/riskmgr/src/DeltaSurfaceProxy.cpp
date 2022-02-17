//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaSurfaceProxy.cpp
//
//   Description : Surface delta - ATM vol rides surface, smile remains constant
//                 for fund proxies
//
//   Author      : Andrew J Swain
//
//   Date        : 8 March 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DeltaSurfaceProxy.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DeltaProxy.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

DeltaSurfaceProxy::Shift::~Shift(){} // empty
DeltaSurfaceProxy::RestorableShift::~RestorableShift(){} // empty

const string DeltaSurfaceProxy::NAME = "DELTA_SURFACE_PROXY";
const string DeltaSurfaceProxy::SECOND_ORDER_NAME = "GAMMA_SURFACE_PROXY";
const double DeltaSurfaceProxy::DEFAULT_SHIFT = 0.005;
const double DeltaSurfaceProxy::MINIMUM_SHIFT = 0.0001;

/** constructor with explicit shift size */
DeltaSurfaceProxy::DeltaSurfaceProxy(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize) {
}

/** for reflection */
DeltaSurfaceProxy::DeltaSurfaceProxy(): 
    ScalarShift(TYPE, NAME) {
}

/** Scales DeltaSurfaceProxy and gamma numbers in results object */
void DeltaSurfaceProxy::scaleResult(Results* results,     // (M)
                                    double   scaleFactor) const{
    // scale all results in packetName
    results->scale(NAME, scaleFactor);
    results->scale(SECOND_ORDER_NAME, scaleFactor);
}


/** Adds DeltaSurfaceProxy and gamma numbers in results object */
void DeltaSurfaceProxy::addResult(Results*       results,     // (M)
                                  const Results* resultsToAdd,
                                  double         scaleFactor) const{
    // add all results in packetName
    results->add(NAME, resultsToAdd, scaleFactor);
    results->add(SECOND_ORDER_NAME, resultsToAdd, scaleFactor);
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double DeltaSurfaceProxy::divisor() const{
    static const char routine[] = "DeltaSurfaceProxy::divisor";
    double initSpot;
    double shiftSize;
    try{
        // find initial value of spot - stored in sens control
        try {
            initSpot = getInitialValue();            
        }
        catch (ModelException& e){
            e.addMsg("Initial value of spot price has not been set in the " 
                     "DeltaSurfaceProxy calculation");
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
void DeltaSurfaceProxy::calculate(TweakGroup* tweakGroup, Results* results){
    try {
        // surface delta PROXY is actually what we get from shifting MINUS delta PROXY
        // see if delta is asked for, and if so, use its shift size
        double deltaShift = control->getDeltaShiftSize();
        double shiftSize  = getShiftSize();
        double shiftToUse = Maths::isZero(deltaShift) ? shiftSize : deltaShift;

        SensitivitySP delta(new DeltaProxy(shiftToUse));
        SensitivityArraySP sens(new SensitivityArray(1, delta));
        OutputRequestArraySP req(new OutputRequestArray(0));

        CControlSP ctrl(new Control(sens, req, false, ""));

        // create new tweak group as we want to use a clean model before
        // we price for regular delta
        IModelSP   model(copy(tweakGroup->getModel()));
        TweakGroup newTweakGroup(tweakGroup->getInstrumentSP(), model);

        ctrl->calculate(newTweakGroup.getModel(), 
                        newTweakGroup.getInstrument(),
                        results);

        // now get raw delta surface
        Results rawResult;

        // need the price for a 2-sided derivative
        rawResult.storePrice(results->retrievePrice(), results->getCcyName());

        calculateTwoSidedDeriv(SECOND_ORDER_NAME, tweakGroup, &rawResult);

        // now subtract the two for what we call surface delta
        const string deltaPacket = delta->getPacketName();

        OutputNameArraySP deltaNames = OutputName::trim(results->packetContents(deltaPacket));
        if (deltaNames->empty()) {
            results->storeNotApplicable(this);            
            results->storeNotApplicable(SECOND_ORDER_NAME);            
        }

        for (int j=0; j < (int)deltaNames->size(); j++) {
            OutputNameConstSP name = (*deltaNames)[j];

            IObjectConstSP deltObj(results->retrieveGreek(deltaPacket,name));
            IObjectConstSP rawObj(rawResult.retrieveGreek(NAME,name));

            bool deltaOK = CDouble::TYPE->isInstance(deltObj.get());
            bool surfdOK = CDouble::TYPE->isInstance(rawObj.get());

            if (deltaOK && surfdOK) {
                double delt = results->retrieveScalarGreek(deltaPacket,name);
                double raw  = rawResult.retrieveScalarGreek(NAME, name);
                results->storeScalarGreek(raw-delt, NAME, name);
                
                double gamma    = results->retrieveScalarGreek(DeltaProxy::SECOND_ORDER_NAME,name);
                double rawGamma = rawResult.retrieveScalarGreek(SECOND_ORDER_NAME, name);
                results->storeScalarGreek(rawGamma-gamma, SECOND_ORDER_NAME, name);
            }
            else {
                // handle NotApplicable & Untweakable
                IObjectSP o(copy(deltaOK ? rawObj.get() : deltObj.get()));
                results->storeGreek(o, NAME, name);
                results->storeGreek(o, SECOND_ORDER_NAME, name);
            }
        }
    } 
    catch (exception& e){
        results->storeGreek(IObjectSP(new Untweakable(e)),
                            NAME, OutputNameSP(new OutputName("")));
        results->storeGreek(IObjectSP(new Untweakable(e)),
                             SECOND_ORDER_NAME, OutputNameSP(new OutputName("")));
    }
}


/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP DeltaSurfaceProxy::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP DeltaSurfaceProxy::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool DeltaSurfaceProxy::nameMatches(const OutputName&     name,
                                    IObjectConstSP        obj){
    // cast obj to DeltaSurfaceProxy::Shift and then invoke name method
    const Shift& proxyObj = dynamic_cast<const Shift&>(*obj);
    return proxyObj.sensNameMatches(this, name);
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void DeltaSurfaceProxy::appendName(OutputNameArray&      namesList,
                                   IObjectConstSP        obj){
    const IShift& proxyObj = dynamic_cast<const IShift&>(*obj);
    proxyObj.sensAppendName(this, namesList);
}

bool DeltaSurfaceProxy::shift(IObjectSP obj) {
    // cast obj to DeltaSurfaceProxy::Shift and then invoke shift method
    Shift& proxyObj = dynamic_cast<Shift&>(*obj);
    return proxyObj.sensShift(this);
}

    
void DeltaSurfaceProxy::restore(IObjectSP obj) {
    // cast obj to DeltaSurfaceProxy::Shift and then invoke restore method
    RestorableShift& proxyObj = dynamic_cast<RestorableShift&>(*obj);
    proxyObj.sensRestore(this);
}

/** Returns the shift which has been made for the current pricing
    call */
ScalarShiftArray DeltaSurfaceProxy::getComponentShifts() const{
    // need array of ScalarShiftConstSP to make this work - but this breaks
    // the array template - would need const array template...
    DeltaSurfaceProxy* DeltaSurfaceProxyShift = const_cast<DeltaSurfaceProxy*>(this);
    return ScalarShiftArray(1, ScalarShiftSP::attachToRef(DeltaSurfaceProxyShift));
}

/** make a surface delta shift using own shift */
DeltaSurface* DeltaSurfaceProxy::makeDeltaSurfaceShift() const {
    static const string method = "DeltaSurfaceProxy::makeDeltaSurfaceShift";
    try {
        DeltaSurfaceSP delta(new DeltaSurface(getShiftSize()));

        if (hasOverrideNames()) {
            delta->storeOverrideNames(toTweak);
        }
        delta->setMarketDataName(getMarketDataName());
        return delta.release();
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}



class DeltaSurfaceProxyHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new DeltaSurfaceProxy(DeltaSurfaceProxy::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new DeltaSurfaceProxy(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaSurfaceProxy, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultDeltaSurfaceProxy);
        // register how to build our sensitivity
        SensitivityFactory::addSens(DeltaSurfaceProxy::NAME, 
                                    new Factory(), 
                                    new DeltaSurfaceProxy(DeltaSurfaceProxy::DEFAULT_SHIFT),
                                    DeltaSurfaceProxy::Shift::TYPE);
    }

    static IObject* defaultDeltaSurfaceProxy(){
        return new DeltaSurfaceProxy();
    }
};

CClassConstSP const DeltaSurfaceProxy::TYPE = CClass::registerClassLoadMethod(
    "DeltaSurfaceProxy", typeid(DeltaSurfaceProxy), DeltaSurfaceProxyHelper::load);

CClassConstSP const DeltaSurfaceProxy::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "DeltaSurfaceProxy::Shift", typeid(DeltaSurfaceProxy::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(DeltaSurfaceProxy::RestorableShift, clazz);
    EXTENDS(DeltaSurfaceProxy::Shift);
}
    
CClassConstSP const DeltaSurfaceProxy::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "DeltaSurfaceProxy::RestorableShift", typeid(DeltaSurfaceProxy::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
