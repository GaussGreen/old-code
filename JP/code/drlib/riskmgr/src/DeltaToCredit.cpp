//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaToCredit.cpp
//
//   Description : DeltaToCredit sensitivity
//
//   Author      : André Segger
//
//   Date        : 15 Sep 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DeltaToCredit.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/E2CModel.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
DeltaToCredit::Shift::~Shift(){} // empty
DeltaToCredit::RestorableShift::~RestorableShift(){} // empty

/** Sens Control for vega parallel */
const string DeltaToCredit::NAME = "DELTA_TO_CREDIT";
const string DeltaToCredit::SECOND_ORDER_NAME = "GAMMA_TO_CREDIT";
const double DeltaToCredit::DEFAULT_SHIFT = 0.005;
const double DeltaToCredit::MINIMUM_SHIFT = 0.0001;

/** constructor with explicit shift size */
DeltaToCredit::DeltaToCredit(double     shiftSize):
    ScalarShift(TYPE, NAME, shiftSize) {
}

/** Constructor with explicit name and shift size, useful for things
    like DeltaNextDay, e.g. */
DeltaToCredit::DeltaToCredit(double     shiftSize,
                             string     newName) :
    ScalarShift(TYPE, newName, shiftSize) {}

/** constructor with explicit shift size */
DeltaToCredit::DeltaToCredit(double     shiftSize,
                             IModel*    model,
                             Control*   control):
    ScalarShift(TYPE, NAME, shiftSize) {
    this->algorithm = model;
    this->control   = control;
}


/** for reflection */
DeltaToCredit::DeltaToCredit(): 
    ScalarShift(TYPE, NAME) {
}

/** Questionable whether we should derive from this class - CrossGamma
    is currently */
DeltaToCredit::DeltaToCredit(const CClassConstSP& clazz,
                             const string&        outputName,
                             const double&        shiftSize): 
    ScalarShift(clazz, outputName, shiftSize) {}

/** Scales delta and gamma numbers in results object */
void DeltaToCredit::scaleResult(Results*     results,     // (M)
                                double       scaleFactor) const {
    // scale all results in packetName
    results->scale(NAME, scaleFactor);
    results->scale(SECOND_ORDER_NAME, scaleFactor);
}


/** Adds delta and gamma numbers in results object */
void DeltaToCredit::addResult(Results*           results,     // (M)
                              const Results*     resultsToAdd,
                              double             scaleFactor) const{
    // add all results in packetName
    results->add(NAME, resultsToAdd, scaleFactor);
    results->add(SECOND_ORDER_NAME, resultsToAdd, scaleFactor);
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double DeltaToCredit::divisor() const{
    static const char routine[] = "DeltaToCredit::divisor";
    double initSpot;
    double shiftSize;
    try{
        // find initial value of spot - stored in sens control
        try {
            initSpot = getInitialValue();            
        }
        catch (ModelException& e){
            e.addMsg("Initial value of spot price has not been set in the " 
                     "DeltaToCredit calculation");
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
void DeltaToCredit::calculate(TweakGroup*      tweakGroup,
                              CResults*        results){
    try{
        // switch to use the correct way of pricing
        if (IE2CModel::TYPE->isInstance(tweakGroup->getModel())) {
            IE2CModel* model = dynamic_cast<IE2CModel*>(tweakGroup->getModel());
            model->setParSpreadPricing(false);
        }

        // get current valuation
        double curVal = results->retrievePrice();

        // get base price
        double basePrice = calcSensPrice(tweakGroup);
        // and use as an override for the
        // purposes of calculating this sensitivity
        results->storePrice(basePrice,results->getCcyName());

        calculateTwoSidedDeriv(SECOND_ORDER_NAME, tweakGroup, results);

        //now restore the price
        results->storePrice(curVal,results->getCcyName());

        // switch to use the correct way of pricing
        if (IE2CModel::TYPE->isInstance(tweakGroup->getModel())) {
            IE2CModel* model = dynamic_cast<IE2CModel*>(tweakGroup->getModel());
            model->setParSpreadPricing(true);
        }
    } catch (exception& e){
        // switch back pricing method for the other greeks
        if (IE2CModel::TYPE->isInstance(tweakGroup->getModel())) {
            IE2CModel* model = dynamic_cast<IE2CModel*>(tweakGroup->getModel());
            model->setParSpreadPricing(true);
        }
        throw ModelException(&e,  "DeltaToCredit::calculate");
    }
}


/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP DeltaToCredit::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP DeltaToCredit::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool DeltaToCredit::nameMatches(const OutputName&         name,
                                IObjectConstSP          obj){
    // cast obj to Delta::Shift and then invoke name method
    const Shift& deltaObj = 
        dynamic_cast<const Shift&>(*obj);
    return name.equals(deltaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void DeltaToCredit::appendName(OutputNameArray&          namesList,
                               IObjectConstSP          obj){
    // cast obj to Delta::Shift and then invoke name method
    const Shift& deltaObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(deltaObj.sensName(this)));

    if (outputName->toString() != "") {
        namesList.push_back(outputName);
    }
}

bool DeltaToCredit::shift(IObjectSP obj) {
    // cast obj to Delta::Shift and then invoke shift method
    Shift& deltaObj = dynamic_cast<Shift&>(*obj);
    return deltaObj.sensShift(this);
}

void DeltaToCredit::restore(IObjectSP obj) {
    // cast obj to Delta::Shift and then invoke restore method
    RestorableShift& deltaObj = 
        dynamic_cast<RestorableShift&>(*obj);
    deltaObj.sensRestore(this);
}

/** Returns the shift which has been made for the current pricing
    call */
ScalarShiftArray DeltaToCredit::getComponentShifts() const{
    // need array of ScalarShiftConstSP to make this work - but this breaks
    // the array template - would need const array template...
    DeltaToCredit* deltaShift = const_cast<DeltaToCredit*>(this);
    return ScalarShiftArray(1, ScalarShiftSP::attachToRef(deltaShift));
}

class DeltaToCreditHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new DeltaToCredit(DeltaToCredit::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new DeltaToCredit(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DeltaToCredit, clazz);
        SUPERCLASS(ScalarShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultDelta);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(DeltaToCredit::NAME, 
                                    new Factory(), 
                                    new DeltaToCredit(DeltaToCredit::DEFAULT_SHIFT),
                                    DeltaToCredit::Shift::TYPE);
    }

    static IObject* defaultDelta(){
        return new DeltaToCredit();
    }
};

CClassConstSP const DeltaToCredit::TYPE = CClass::registerClassLoadMethod(
    "DeltaToCredit", typeid(DeltaToCredit), DeltaToCreditHelper::load);

CClassConstSP const DeltaToCredit::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "DeltaToCredit::Shift", typeid(DeltaToCredit::Shift), 0);

static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(DeltaToCredit::RestorableShift, clazz);
    EXTENDS(DeltaToCredit::Shift);
}
    
CClassConstSP const DeltaToCredit::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "DeltaToCredit::RestorableShift", typeid(DeltaToCredit::RestorableShift),
    restorableShiftLoad);

DRLIB_END_NAMESPACE
