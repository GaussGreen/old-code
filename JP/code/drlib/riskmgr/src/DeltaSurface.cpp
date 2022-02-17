//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaSurface.cpp
//
//   Description : Surface delta - ATM vol rides surface, smile remains constant
//
//   Author      : Andrew J Swain
//
//   Date        : 24 February 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Delta.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"

DRLIB_BEGIN_NAMESPACE

DeltaSurface::Shift::Shift(){} // empty
DeltaSurface::Shift::~Shift(){} // empty
DeltaSurface::RestorableShift::~RestorableShift(){} // empty

DeltaSurface::MyString::MyString(const string& s): string(s){}
DeltaSurface::MyString::MyString(): string(){}

const string DeltaSurface::NAME = "DELTA_SURFACE";
const string DeltaSurface::SECOND_ORDER_NAME = "GAMMA_SURFACE";
const double DeltaSurface::DEFAULT_SHIFT = 0.005;
const double DeltaSurface::MINIMUM_SHIFT = 0.0001;

/** constructor with explicit shift size */
DeltaSurface::DeltaSurface(double shiftSize):
    ScalarShift(TYPE, NAME, shiftSize), spot(0.0), spotSet(false) {
}

/** constructor with explicit shift size */
DeltaSurface::DeltaSurface(double   shiftSize,
                           IModel*  model,
                           Control* control):
    ScalarShift(TYPE, NAME, shiftSize), spot(0.0), spotSet(false) {
    this->algorithm = model;
    this->control   = control;
}


/** for reflection */
DeltaSurface::DeltaSurface(): 
    ScalarShift(TYPE, NAME), spot(0.0), spotSet(false) {
}

/** Scales DeltaSurface and gamma numbers in results object */
void DeltaSurface::scaleResult(Results*     results,     // (M)
                        double       scaleFactor) const{
    // scale all results in packetName
    results->scale(NAME, scaleFactor);
    results->scale(SECOND_ORDER_NAME, scaleFactor);
}


/** Adds DeltaSurface and gamma numbers in results object */
void DeltaSurface::addResult(Results*           results,     // (M)
                             const Results*     resultsToAdd,
                             double             scaleFactor) const{
    // add all results in packetName
    results->add(NAME, resultsToAdd, scaleFactor);
    results->add(SECOND_ORDER_NAME, resultsToAdd, scaleFactor);
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double DeltaSurface::divisor() const{
    static const char routine[] = "DeltaSurface::divisor";
    double initSpot;
    double shiftSize;
    try{
        // find initial value of spot - stored in sens control
        try {
            initSpot = getInitialValue();            
        }
        catch (ModelException& e){
            e.addMsg("Initial value of spot price has not been set in the " 
                     "DeltaSurface calculation");
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
void DeltaSurface::calculate(TweakGroup* tweakGroup, Results* results){
    try{
        // surface delta is actually what we get from shifting MINUS delta
        // see if delta is asked for, and if so, use its shift size
        double deltaShift = control->getDeltaShiftSize();
        double shiftSize  = getShiftSize();
        double shiftToUse = Maths::isZero(deltaShift) ? shiftSize : deltaShift;

        SensitivitySP delta(new Delta(shiftToUse));
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
            if (!name->isEmpty()) {
                IObjectConstSP deltObj(results->retrieveGreek(deltaPacket,name));
                IObjectConstSP rawObj(rawResult.retrieveGreek(NAME,name));
            
                bool deltaOK = CDouble::TYPE->isInstance(deltObj.get());
                bool surfdOK = CDouble::TYPE->isInstance(rawObj.get());

                if (deltaOK && surfdOK) {
                    double delt = results->retrieveScalarGreek(deltaPacket,name);
                    double raw  = rawResult.retrieveScalarGreek(NAME, name);
                    results->storeScalarGreek(raw-delt, NAME, name);
                
                    double gamma    = results->retrieveScalarGreek(Delta::SECOND_ORDER_NAME,name);
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
CClassConstSP DeltaSurface::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP DeltaSurface::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool DeltaSurface::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to DeltaSurface::Shift and then invoke name method
    const Shift& DeltaSurfaceObj = 
        dynamic_cast<const Shift&>(*obj);

    string objName = DeltaSurfaceObj.sensName(this);

    // see if we match on either asset or vol name
    return (name.equals(objName) || name.equals(getEqName(objName)));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void DeltaSurface::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to DeltaSurface::Shift and then invoke name method
    const Shift& DeltaSurfaceObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(DeltaSurfaceObj.sensName(this)));
    namesList.push_back(outputName);
}

bool DeltaSurface::shift(IObjectSP obj) {
    // cast obj to DeltaSurface::Shift and then invoke shift method
    Shift& DeltaSurfaceObj = dynamic_cast<Shift&>(*obj);
    return DeltaSurfaceObj.sensShift(this);
}

void DeltaSurface::restore(IObjectSP obj) {
    // cast obj to DeltaSurface::Shift and then invoke restore method
    RestorableShift& DeltaSurfaceObj = 
        dynamic_cast<RestorableShift&>(*obj);
    DeltaSurfaceObj.sensRestore(this);
}

/** Returns the shift which has been made for the current pricing
    call */
ScalarShiftArray DeltaSurface::getComponentShifts() const{
    // need array of ScalarShiftConstSP to make this work - but this breaks
    // the array template - would need const array template...
    DeltaSurface* deltaShift = const_cast<DeltaSurface*>(this);
    return ScalarShiftArray(1, ScalarShiftSP::attachToRef(deltaShift));
}

/** Returns spot value of associated asset. Fails if spot has not been
    set */
double DeltaSurface::getSpot() const{
    if (!spotSet){
        throw ModelException("DeltaSurface::getSpot", "Spot level not set");
    }
    return spot;
}

/** Stores the spot value of the associated asset. Allows vol to 
    implement skew tweak */
void DeltaSurface::setSpot(double spotValue){
    spot = spotValue;
    spotSet = true;
}

void DeltaSurface::setSpot(double        spotValue,   // asset spot
                           const string& eqName,      // asset name
                           const string& volName)     // asset's vol name
{
    setSpot(spotValue);
    volNameMap[volName] = eqName;
}

string DeltaSurface::getEqName(const string& volName) {
    TweakNameMap::const_iterator iter = volNameMap.find(volName);
    if (iter == volNameMap.end()) {
        return "";
    }
    return (iter->second); 
}

class DeltaSurfaceHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new DeltaSurface(DeltaSurface::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new DeltaSurface(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaSurface, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultDeltaSurface);
        FIELD(spot, "spot");
        FIELD_MAKE_TRANSIENT(spot);
        FIELD(spotSet, "spotSet");
        FIELD_MAKE_TRANSIENT(spotSet);
        // register how to build our sensitivity
        SensitivityFactory::addSens(DeltaSurface::NAME, 
                                    new Factory(), 
                                    new DeltaSurface(DeltaSurface::DEFAULT_SHIFT),
                                    DeltaSurface::Shift::TYPE);
    }

    static IObject* defaultDeltaSurface(){
        return new DeltaSurface();
    }
};

CClassConstSP const DeltaSurface::TYPE = CClass::registerClassLoadMethod(
    "DeltaSurface", typeid(DeltaSurface), DeltaSurfaceHelper::load);

CClassConstSP const DeltaSurface::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "DeltaSurface::Shift", typeid(DeltaSurface::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(DeltaSurface::RestorableShift, clazz);
    EXTENDS(DeltaSurface::Shift);
}
    
CClassConstSP const DeltaSurface::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "DeltaSurface::RestorableShift", typeid(DeltaSurface::RestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
