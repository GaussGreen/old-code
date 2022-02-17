//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Sensitivity.cpp
//
//   Description : Base class for sensitivities
//
//   Author      : Mark A Robson
//
//   Date        : 11 Jun 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_SENSITIVITY_CPP
#include "edginc/Addin.hpp"
#include "edginc/Results.hpp"
#include "edginc/NotApplicable.hpp"
#include "edginc/TweakGroup.hpp"

DRLIB_BEGIN_NAMESPACE

Sensitivity::~Sensitivity(){}

/** identifies the packet in which the results are stored. By default this
    returns the same as getSensOutputName() */
const string& Sensitivity::getPacketName() const{
    return getSensOutputName();
}

void Sensitivity::calculateSens(IModel*          algorithm,
                                CInstrument*     instrument,
                                Control*         control,
                                Results*         results){
    static char routine[] = "Sensitivity::calculateSens";
    /* model value is used not only by pricing calls within the tweaking
       layer but also by the actual shifting routines */
    if (!algorithm){
        throw ModelException(routine, "Null Model");
    }
    if (!results){
        throw ModelException(routine, "Null Results");
    }
    if (!instrument){
        throw ModelException(routine, "Null Instrument");
    }
    this->algorithm = algorithm;
    this->control   = control;
    // todo: review SP part
    TweakGroup tweakGroup(InstrumentSP::attachToRef(instrument),
                          IModelSP::attachToRef(algorithm));
    try {
        calculate(&tweakGroup, results);
    } catch (exception& e){
        this->algorithm = 0;
        this->control   = 0;
        throw ModelException(e, routine, "Failed to calculate "+
                             getSensOutputName());
    }
    this->algorithm = 0;
    this->control   = 0;
}

/** returns a reference to the algorithm used for pricing */
IModel* Sensitivity::getModel() const{
    if (!algorithm){
        throw ModelException("Sensitivity::getModel", "Null Model");
    }
    return algorithm;
}

/** returns a reference to the control */
Control* Sensitivity::getControl() const{
    if (!control){
        throw ModelException("Sensitivity::getControl", "Null Control");
    }
    return control;
}

/** scale the results in the Results Object for this sensitivity
    by supplied factor */
void Sensitivity::scaleResult(Results*     results,     // (M)
                              double       scaleFactor) const{
    // can we combine this sensitivity ?
    if (Additive::TYPE->isInstance(this)) {
        const string& packetName = getPacketName();
        const string& sensName = getSensOutputName();
        if (packetName == sensName){
            // scale all results in packetName
            results->scale(packetName, scaleFactor);
        } else {
            // scale result sensName in packet packetName
            results->scale(packetName, sensName, scaleFactor);
        }
    }
}
        

/** Modify results in the Results Object for this sensitivity by
    adding all results in resultsToAdd as indicated by control */
void Sensitivity::addResult(Results*           results,     // (M)
                            const Results*     resultsToAdd,
                            double             scaleFactor) const{
    // can we combine this sensitivity ?
    if (Additive::TYPE->isInstance(this)) {
        const string& packetName = getPacketName();
        const string& sensName = getSensOutputName();
        if (packetName == sensName){
            // add all results in packetName
            results->add(packetName, resultsToAdd, scaleFactor);
        } else {
            // add result sensName in packet packetName
            results->add(packetName, sensName, resultsToAdd, scaleFactor);
        }
    } else {
        results->storeNotApplicable(this, NotApplicable::aggregation);
    }
}

/** Extracts 'price' for this sensitivity  from results set */
double Sensitivity::getSensPrice(Results*      results,
                                 CInstrument*  inst,
                                 IModel*        model,
                                 Control*      control) {
    return results->retrievePrice();
}
 
/** populate from market cache - default provided */
void Sensitivity::getMarket(const IModel* model, const MarketData* market) {
    // do nothing
}

/** does this sensitivity have its own list of underlyings ? */
bool Sensitivity::hasOverrideNames() const {
    if (!toTweak || toTweak->size() == 0) {
        return false;
    }
    return true;
}

/** return its own list of underlyings to tweak */
OutputNameArrayConstSP Sensitivity::overrideNames() const {
    if (!hasOverrideNames()) {
        throw ModelException("Sensitivity::overrideNames",
                             getSensOutputName() + " has no underlying names "
                             "to tweak");
    }

    return toTweak;
}
  
/** store list of underlyings to tweak */
void Sensitivity::storeOverrideNames(OutputNameArraySP names) {
    toTweak = names;
}

/** specialised copy that allows over-ride of algorithm */
Sensitivity * Sensitivity::spawn(IModel* model) const {
    SensitivitySP minime(copy(this));
    minime->algorithm = model;
    return minime.release();
}

/** Removes any names in the array returned by names() that have
    already been calculated for the specified packet. Note that this
    object is modified accordingly */
void Sensitivity::removeOverrideNames(const string&  packetName,
                                      const Results* results){
    if (hasOverrideNames()){
        for (vector<OutputNameSP>::iterator iter = toTweak->begin();
             iter != toTweak->end(); /* ++ in loop*/){
            if (results->exists(packetName, *iter)){
                iter = toTweak->erase(iter);
            } else {
                ++iter;
            }
        }
    }
}
    
Sensitivity::Sensitivity(const CClassConstSP& clazz):
    CObject(clazz), algorithm(0), control(0), toTweak(0) {}

class SensitivityHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spread sheet
        REGISTER(Sensitivity, clazz);
        SUPERCLASS(CObject);
        FIELD(toTweak, "underlyings to tweak");
        FIELD_MAKE_OPTIONAL(toTweak);
    }
};

CClassConstSP const Sensitivity::TYPE = CClass::registerClassLoadMethod(
    "Sensitivity", typeid(Sensitivity), SensitivityHelper::load);

DEFINE_TEMPLATE_TYPE(SensitivityArray);


class ToTweakAddin: public CObject{
    static CClassConstSP const TYPE;

    SensitivitySP sens;
    CStringArray  names;
    CStringArray  names2;

    /** the 'addin function' */
    static IObjectSP addNames(ToTweakAddin* params){
        static const string method = "ToTweakAddin::addNames";
        try {
            bool pairs = false;

            OutputNameArraySP names(new OutputNameArray(params->names.size()));
            
            if (!params->names2.empty()) {
                pairs = true;
                // if creating name pairs check lengths the same
                if (params->names.size() != params->names2.size()) {
                    throw ModelException(method,
                                         "building name pairs but name lists "
                                         "are different lengths");
                }
            }

            for (int i = 0; i < params->names.size(); i++) {
                if (!pairs) {
                    (*names)[i] = OutputNameSP(new OutputName(params->names[i]));
                }
                else {
                    (*names)[i] = OutputNameSP(new OutputName(params->names[i],
                                                              params->names2[i]));
                }                    
            }

            IObjectSP sens(params->sens->clone());

            SensitivitySP sensSP = SensitivitySP::dynamicCast(sens);
            sensSP->storeOverrideNames(names);

            return sens;
        } 
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    ToTweakAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ToTweakAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultToTweakAddin);
        FIELD(sens, "sensitivity");
        FIELD(names, "underlying names");
        FIELD(names2, "second underlying names");
        FIELD_MAKE_OPTIONAL(names2);

        Addin::registerClassObjectMethod("SENS_TWEAK_NAMES",
                                         Addin::RISK,
                                         "add tweak names to sensitivity",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addNames);
    }

    static IObject* defaultToTweakAddin(){
        return new ToTweakAddin();
    }   
};

CClassConstSP const ToTweakAddin::TYPE = CClass::registerClassLoadMethod(
    "ToTweakAddin", typeid(ToTweakAddin), load);


void Sensitivity::setControl(CControl* c) {
    control = c;
}

void Sensitivity::setAlgorithm(IModel* m) {
    algorithm = m;
}

DRLIB_END_NAMESPACE
