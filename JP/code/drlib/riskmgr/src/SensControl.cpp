
#include "edginc/config.hpp"
#define QLIB_SENSCONTROL_CPP
#include "edginc/ScalarTweak.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/VectorTweak.hpp"

DRLIB_BEGIN_NAMESPACE

IScalarTweak::IScalarTweak(){} // declaration in ScalarTweak.hpp
IScalarTweak::~IScalarTweak(){}
IVectorTweak::IVectorTweak(){}  // declaration in VectorTweak.hpp
IVectorTweak::~IVectorTweak(){}  // declaration in VectorTweak.hpp

SensControl::~SensControl(){}

/** identifies the name used storing associated results in the output */
const string& SensControl::getSensOutputName() const{
    return outputName;
}

/** stores an the intitial value of the quantity being shifted. This is
    very useful in calculating greeks like delta where the original value
    of what was being shifted is needed. It avoid the need for a 
    separate suite of interfaces allowing the retrieval of spot prices,
    correlations, fx spots etc. */
void SensControl::setInitialValue(double val){
    if (!initialValSet){
        if (initialValue == 0.0){
            initialValSet = true;
            initialValue = val;
        }
    } else {
        // make sure the value is consistent
        if (!Maths::equals(val, initialValue)){
            initialValSet = false;
            initialValue = 1.0; // any value other than zero will do
        }
    }
}

/** gets the intitial value of the quantify being shifted. See the 
    comments under {@link #initialValueSet(double) initialValueSet} */
double SensControl::getInitialValue() const{
    static const string method("SensControl::getInitialValue");
    if (!initialValSet){
        if (initialValue == 0.0){
            throw ModelException(method, "Initial value not set");
        }
        //it is useful to note how the exception below is triggered
        //-----------------------------------------------------------------
        //when a tweak is applied, it is done to objects by name and type
        //the object hierarchy (instrument, model & fetched market data)
        //is traversed for all matches (on name and type).
        //For the tweak to be consistent, the data being shifed in each
        //matched object must be consistent.
        //The initialValSet flag is managed in the "reset" method and the
        //"setInitialValue" method
        //Typically, "reset" will be called before the shift is applied to
        //the matching objects
        //Then each object is shifted in turn (at some point calling
        //"setInitialValue") and therefore inconsistent initial values cause
        //initialValSet to be false in that method and consequently trigger
        //this exception later on
        throw ModelException(method, "Inconsistent values set for initial"
                             " value");
    }
    return initialValue;
}

/** Used by objects implementing restorable tweaks to store data needed
    to restore themselves after the tweak. Stores given data on
    the SensControl */
void SensControl::setRestoreData(const IObjectSP& restoreData){
    this->restoreData = restoreData;
}

/** Used by objects implementing restorable tweaks to store data needed
    to restore themselves after the tweak. Retrieves the data previously
    stored on the SensControl. The data stored is then cleared */
IObjectSP SensControl::getRestoreData(){
    if (!restoreData){
        return restoreData;
    } else {
        IObjectSP data(restoreData);
        restoreData = IObjectSP(   );
        return data;
    }
}

/** Implemented by restorableShiftInterface()->isInstance(obj) */
bool SensControl::restorableShift(IObjectConstSP obj) const{
    return restorableShiftInterface()->isInstance(obj);
}

/** shifts instrument/market data by given shift, reprices, and
    restores data to original. Returns shifted price */
double SensControl::shiftAndPrice(TweakGroup*      tweakGroup,
                                  double           origPrice){
    // Dummy try/catch to avoid vc6.opt crashes
    try {
        return shiftAndPrice(tweakGroup, origPrice, false);
    } catch (...) { throw; }
}

/** shifts instrument/market data by given shift, reprices, and
    restores data to original. Returns shifted price. If forceReprice
    is true then price will be recalculated regardless of whether
    instrument is actually shifted during the tweak */
double SensControl::shiftAndPrice(TweakGroup*      tweakGroup,
                                  double           origPrice,
                                  bool             forceReprice){
    static char   routine[] = "SensControl::shiftAndPrice";
    //// performance optimisation if ctrl->skipShifts() is true
    Control* ctrl = getControl(); // does not return null
    if (ctrl->skipShifts()){
        initialValue = 1.0; // to stop eg delta from failing
        initialValSet = true;
        return calcSensPrice(tweakGroup); // still required to 'price' though
    }
    SensMgrOpt    sensMgr(tweakGroup);
    TweakGroupSP  shiftedGroup;
    try{
        shiftedGroup = TweakGroupSP::dynamicCast(sensMgr.shift(this));
    } catch (exception& e){
        throw ModelException(&e, routine);
    }
    double  shiftedPrice;
    try{
        if (!forceReprice && !sensMgr.getShiftStatus()){
            shiftedPrice = origPrice;
        } else {
            shiftedPrice = calcSensPrice(shiftedGroup.get());
        }
    } catch (exception& e){
        sensMgr.restore(); // must restore instrument after successful shift
        throw ModelException(&e, routine);
    }
    sensMgr.restore();  // must restore instrument after successful shift
    return shiftedPrice;
}

/** shifts instrument/market data by this shift and then by given
    shift, reprices, and restores data to original. Returns shifted
    price */
double SensControl::shiftAndPrice(SensControl*     secondShift,
                                  TweakGroup*      tweakGroup,
                                  double           origPrice){
    static char    routine[] = "SensControl::shiftAndPrice";
    //// performance optimisation if ctrl->skipShifts() is true
    Control* ctrl = getControl(); // does not return null
    if (ctrl->skipShifts()){
        initialValue = 1.0; // to stop eg delta from failing
        initialValSet = true;
        return calcSensPrice(tweakGroup); // still required to 'price' though
    }
    SensMgrOpt     sensMgr(tweakGroup);
    TweakGroupSP   shifted1Inst;
    try{
        shifted1Inst = TweakGroupSP::dynamicCast(sensMgr.shift(this));
    } catch (exception& e){
        throw ModelException(&e, routine);
    }
    double  shiftedPrice;
    try{
        bool needReprice = sensMgr.getShiftStatus();
        shiftedPrice = secondShift->shiftAndPrice(shifted1Inst.get(), 
                                                  origPrice,
                                                  needReprice);
    } catch (exception& e){
        sensMgr.restore(); // must restore instrument after successful shift
        throw ModelException(&e, routine);
    }
    sensMgr.restore();    // must restore instrument after successful shift
    return shiftedPrice;
}


/** Applies this shift and then does control->calculate() where the
    control is constructed using sens and outReqs, either or both 
    of which can be null */
ResultsSP SensControl::shiftAndCalculate(
    TweakGroup*                      tweakGroup,
    const SensitivityArrayConstSP&   sens,
    const OutputRequestArrayConstSP& outReqs){
    static char    routine[] = "SensControl::shiftAndCalculate";
    SensMgrOpt    sensMgr(tweakGroup);
    TweakGroupSP  shiftedInst;
    ResultsSP results(new Results);
    bool skipShifts = control && control->skipShifts(); // control might be 0
    try{
        shiftedInst = skipShifts? 
            TweakGroupSP::attachToRef(tweakGroup):
            TweakGroupSP::dynamicCast(sensMgr.shift(this));
    } catch (exception& e){
        throw ModelException(e, routine);
    }
    try{
        SensitivityArrayConstSP mySens = !sens?
            SensitivityArrayConstSP(new SensitivityArray()): sens;
        OutputRequestArrayConstSP myOutReqs = !outReqs?
            OutputRequestArrayConstSP(new OutputRequestArray(0)): outReqs;
        Control myControl(mySens, myOutReqs, false, "", 
                          skipShifts); // pass on flag
        IObjectSP modelSP(shiftedInst->getModel()->clone());
        IModelSP model = IModelSP::dynamicCast(modelSP);
        myControl.calculate(model.get(), shiftedInst->getInstrument(),
                            results.get());
    } catch (exception& e){
        sensMgr.restore(); // must restore instrument after successful shift
        throw ModelException(&e, routine);
    }
    sensMgr.restore();  // must restore instrument after successful shift
    return results;
}

/** Like calculateSens in terms of inputs but calculates a single
    first order one-sided scalar derivative. The name of what to
    be tweaked should be set in 'this' */
double SensControl::calculateOneSidedFirstDeriv(IModel*          algorithm,
                                                CInstrument*     instrument,
                                                Control*         control,
                                                Results*         results){
    static char routine[] = "SensControl::calculateOneSidedFirstDeriv";
    if (!algorithm || !results || !instrument || !control){
        throw ModelException(routine, "Null Inputs");
    }
    this->algorithm = algorithm;
    this->control   = control;
    TweakGroup tweakGroup(InstrumentSP::attachToRef(instrument),
                          IModelSP::attachToRef(algorithm));
    try {
        double deriv = calcOneSidedFirstDeriv(&tweakGroup, results);
        this->algorithm = 0;
        this->control   = 0;
        return deriv;
    } catch (exception& e){
        this->algorithm = 0;
        this->control   = 0;
        throw ModelException(e, routine, "Failed to calculate "+
                             getSensOutputName());
    }
}

/** calculates a first order one-sided scalar derivative */
double SensControl::calcOneSidedFirstDeriv(
    TweakGroup*    tweakGroup,
    CResults*      results){
    // Dummy try/catch to avoid vc6.opt crashes
    try {
        double origPrice = getSensPrice(results, tweakGroup->getInstrument(), tweakGroup->getModel(), getControl());
        double shiftedPrice = shiftAndPrice(tweakGroup, origPrice);
        double divisor = this->divisor();
        return (shiftedPrice - origPrice)/divisor;
    } catch (...) { throw; }
}

/** Calculates 'price' for this sensitivity */
double SensControl::calcSensPrice(TweakGroup*  tweakGroup){
    Results results;
    tweakGroup->getModel()->Price(tweakGroup->getInstrument(), getControl(), 
                                  &results);
    return results.retrievePrice();
    // return getSensPrice(&results, tweakGroup->getInstrument(), tweakGroup->getModel(), getControl());
}

/** resets all internally stored values */
void SensControl::reset() {
    initialValue   = 0.0;
    initialValSet  = false;
    restoreData    = IObjectSP(   );
}

SensControl::SensControl(CClassConstSP clazz, const string& outputName):
    Sensitivity(clazz),
    outputName(outputName),
    initialValue(0.0), initialValSet(false),
    restoreData(0) {
}

class SensControlHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SensControl, clazz);
        SUPERCLASS(Sensitivity);
        FIELD(initialValue, "Initial value");
        FIELD_MAKE_TRANSIENT(initialValue);
        FIELD(initialValSet, "Initial value set");
        FIELD_MAKE_TRANSIENT(initialValSet);
        FIELD(restoreData, "restore date");
        FIELD_MAKE_TRANSIENT(restoreData);
        /* note no fields are stored - here it makes sense for the default
           constructors of inherited classes to correctly build this object */
    }
};

CClassConstSP const SensControl::TYPE = CClass::registerClassLoadMethod(
    "SensControl", typeid(SensControl), SensControlHelper::load);

DEFINE_TEMPLATE_TYPE(SensControlArray);

DRLIB_END_NAMESPACE
