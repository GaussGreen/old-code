//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Theta.cpp
//
//   Description : Theta shift
//
//   Author      : Andrew J Swain
//
//   Date        : 16 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/SensControlPerName.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE
Theta::Shift::Shift(){} // empty
Theta::Shift::~Shift(){} // empty
Theta::RestorableShift::~RestorableShift(){} // empty
Theta::INextDaySensitivity::~INextDaySensitivity(){}

/** Sens Control for theta */
const string Theta::NAME = "THETA";
const int Theta::DEFAULT_SHIFT = 1;

/** Simple implementation of above method where you can specify
    a SensControl */
SensitivitySP Theta::INextDaySensitivity::Util::requiredSensitivity(
    const Sensitivity*   nextDaySens,
    SensControlPerNameSP sensCtrl, /* eg Delta when doing DeltaNextDay */
    TweakGroup*          tweakGroup){
    OutputNameArrayConstSP names;
    if (nextDaySens->hasOverrideNames()) {
        names = nextDaySens->overrideNames();
    } else {
        names = sensCtrl->names(tweakGroup);
    }
    sensCtrl->storeOverrideNames(OutputNameArraySP(names.clone()));    
    return SensitivitySP(sensCtrl);
}

/** Simple implementation of above method where
    sensCtrl is the 'this' passed to writeResults. The copyWholePacket 
    determines whether the names within reqdSens is ignored and the whatever
    results are in the packet are just copied over. Otherwise only the names
    listed in reqdSens are copied over */
void Theta::INextDaySensitivity::Util::writeResults(
    const Sensitivity*            nextDaySens,
    const SensitivitySP&          reqdSens,
    bool                          copyWholePacket,
    Results*                      dest,
    const UntweakableSP&          untweakable,
    const Results*                src){
    const string& packetName = nextDaySens->getPacketName();
    if (!untweakable){
        // we've got some results perhaps
        if (!reqdSens->hasOverrideNames() || !src){
            // there was nothing to tweak
            dest->storeNotApplicable(packetName);
        } else {
            // copy them into dest
            if (copyWholePacket){
                vector<pair<OutputNameConstSP, IObjectConstSP> > outputs(
                    src->listPacketResults(reqdSens->getPacketName()));
                for (unsigned int j=0; j < outputs.size(); j++) {
                    IObjectSP output(outputs[j].second->clone());
                    dest->storeCcylessResult(output, packetName,
                                             outputs[j].first);
                }
            } else {
                OutputNameArrayConstSP names(reqdSens->overrideNames());
                const string& srcPacketName = reqdSens->getPacketName();
                for (int i = 0; i < names->size(); i++){
                    const OutputNameSP& name = (*names)[i];
                    IObjectSP myResult(
                        src->retrieveGreek(srcPacketName, name).clone());
                    dest->storeGreek(myResult, packetName, name);
                }
            }
        }
    } else {
        dest->storeGreek(IObjectSP(untweakable), packetName, 
                         OutputNameSP(new OutputName("")));
    }
}


/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns true */
bool Theta::discreteShift() const{
    return true;
}

/** constructors */
Theta::Theta(int offset, HolidaySP hols): SensControlAllNames(TYPE, NAME), 
    offset(offset), hols(hols){
    if (!hols){
        throw ModelException("Theta::Theta", "Holidays are null");
    }
}

Theta::Theta(int offset, HolidaySP hols,
             const CClassConstSP& clazz,
             const string&        outputName): 
            SensControlAllNames(clazz, outputName), offset(offset), hols(hols){
    if (!hols){
        throw ModelException("Theta::Theta", "Holidays are null");
    }
}

Theta::~Theta(){}

Theta::Theta(const CClassConstSP& clazz,
             const string&        outputName): 
    SensControlAllNames(clazz, outputName), offset(0){}

/** identifies the packet in which the results are stored. Theta
    results are stored in the instrument packet */
const string& Theta::getPacketName() const{
    return Results::INSTRUMENT_PACKET;
}

/** for reflection */
Theta::Theta(): SensControlAllNames(TYPE, NAME) {}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double Theta::divisor() const{
    return 1.0;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP Theta::shiftInterface() const{
    return Shift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP Theta::restorableShiftInterface() const{
    return RestorableShift::TYPE;
}

bool Theta::shift(IObjectSP obj) {
    // cast obj to Theta::Shift and then invoke shift method
    Shift& thetaObj = dynamic_cast<Shift&>(*obj);
    return thetaObj.sensShift(this);
}

void Theta::restore(IObjectSP obj) {
    // cast obj to Theta::Shift and then invoke restore method
    RestorableShift& thetaObj = dynamic_cast<RestorableShift&>(*obj);
    thetaObj.sensRestore(this);
}

/** returns date + offset business days */
DateTime Theta::rollDate(const DateTime &date) const{
    return hols->addBusinessDays(date, offset);
}

/** Returns true if assets forward value should be used (ie theta
    forward spot) */
bool Theta::useAssetFwds() const{
    return false;
}

/** writes theta value + results for output requests to result object */
void Theta::writeResults(Results*             results, 
                         const UntweakableSP& untweakable,
                         Results*             shiftedResults){
    OutputNameConstSP thetaName(new OutputName(getSensOutputName()));
    double    divisor;
    if (!untweakable){
        // double origPrice = getSensPrice(results, tweakGroup->getInstrument(), tweakGroup->getModel(), getControl());
        double origPrice = results->retrievePrice();
        double shiftedPrice = shiftedResults->retrievePrice();
        divisor = this->divisor();
        double theta = (shiftedPrice-origPrice)/divisor;
        results->storeScalarGreek(theta, Results::INSTRUMENT_PACKET, thetaName);
    } else {
        results->storeGreek(IObjectSP(untweakable),
                            Results::INSTRUMENT_PACKET, thetaName);
    }
    // then sort out the outRequests
    for (int i = 0; outRequests.get() && i < outRequests->size(); i++) {
        OutputRequestSP& request = (*outRequests)[i];
        const string& packetName = request->getPacketName();
        const string& requestName = request->getRequestName();
        OutputNameSP outputName(new OutputName(requestName));
        OutputNameSP newOutputName(new OutputName(requestName + "_" + 
                                                  getSensOutputName()));
        try {
            IObjectSP requestResult;
            if (!untweakable){
                // cope with original calculation not producing a number
                IObjectConstSP origResult;
                if (results->exists(request.get())){
                    origResult = results->retrieveGreek(packetName, outputName);
                }
                if (!origResult || 
                    !CDouble::TYPE->isInstance(origResult)){
                    if (NotApplicable::TYPE->isInstance(origResult)){
                        requestResult = IObjectSP(new NotApplicable());
                    } else {
                        string reason = !origResult?
                            "calculation not requested": "calculation failed";
                        requestResult = IObjectSP(
                            new Untweakable(requestName+": Original "+reason));
                    }
                } else {
                    double origValue = 
                        results->retrieveScalarGreek(packetName, outputName);
                    IObjectConstSP newResult(
                        results->retrieveGreek(packetName, outputName));
                    if (Untweakable::TYPE->isInstance(newResult)){
                        requestResult = IObjectSP(newResult.clone());
                    } else {
                        double newValue = 
                            shiftedResults->retrieveScalarGreek(packetName, 
                                                                outputName);
                        double theta = (newValue - origValue)/divisor;;
                        requestResult = IObjectSP(CDouble::create(theta));
                    }
                }
            } else {
                requestResult = untweakable;
            }
            results->storeGreek(requestResult, Results::INSTRUMENT_PACKET,
                                newOutputName);
        } 
        catch (exception& e){
            results->storeGreek(IObjectSP(new Untweakable(e)),
                                Results::INSTRUMENT_PACKET,
                                newOutputName);
        }
    }
}

/** calculates given sensitivity - invoked by calculateSens */
void Theta::calculate(TweakGroup*  tweakGroup,
                      CResults*    results) {
    /* whilst we're doing this type of theta do any other next day
     * type sens that uses the same type of theta */
    calculateNextDaySens(getClass(), control, tweakGroup, results);
}

/** Potentially calculates 'nexy day' sensitivities corresponding to the
    specified thetaType (eg theta, delta next day, delta proxy next day etc)
    depending on what's been asked for in the control. Optimises things by
    doing everything together. Does not repeat calculation if result has
    already been calculated. */
void Theta::calculateNextDaySens(const CClassConstSP& thetaType,
                                 Control*             control,
                                 TweakGroup*          tweakGroup,
                                 CResults*            results){
    try{
        SensitivityArrayConstSP allSens(control->getSens());
        // find all the NextDay greeks and theta
        vector<const INextDaySensitivity*> nextDayGreek;
        ThetaSP theta;
        UntweakableSP untweakable;
        for (int i = 0; i < allSens->size(); i++){
            const SensitivitySP& sens = (*allSens)[i];
            // pull off first theta with correct type
            if (!theta && sens->getClass() == thetaType){
                theta = ThetaSP::dynamicCast(sens);
            } else if (INextDaySensitivity::TYPE->isInstance(sens.get())){
                nextDayGreek.push_back(
                    &dynamic_cast<const INextDaySensitivity&>(*sens));
            }
        }
        // see if theta already calculated
        if (theta.get()){
            OutputNameConstSP outputName(
                new OutputName(theta->getSensOutputName()));
            if (results->exists(Results::INSTRUMENT_PACKET, outputName)) {
                theta.reset(); // clear it
            }
        }
        while (theta.get() || !nextDayGreek.empty()){
            // need to aggregate those greeks with same offset/hols
            //vector<const INextDaySensitivity*> nextDayGreekToDo;
            SensitivityArraySP nextDaySensToDo(new SensitivityArray());
            vector<const INextDaySensitivity*> nextDayGreekToDo;
            ThetaSP rollData; // how to roll forward to next day
            if (!theta){
                rollData = nextDayGreek[0]->offsetRequired(tweakGroup->getInstrument());
            } else {
                rollData = theta;
            }
            // then find others with matching offset/hols
            for (vector<const INextDaySensitivity*>::iterator iter = 
                     nextDayGreek.begin();
                 iter != nextDayGreek.end(); /* ++ in loop */){
                // get hold of possible alternative roll data
                ThetaSP altRollData((*iter)->offsetRequired(tweakGroup->getInstrument())); 
                if (altRollData->getClass() != thetaType){
                    // skip this one - it uses a different type of theta
                    iter = nextDayGreek.erase(iter);
                } else if (altRollData->equals(rollData.get())){
                    // get hold of sensitivity to calc on next day
                    SensitivitySP sens(
                        (*iter)->requiredSensitivity(tweakGroup).clone());
                    if (!sens->hasOverrideNames()){
                        // get N/A stored
                        (*iter)->writeResults(sens, results, untweakable, 0);
                    } else {
                        // then remove names that have been calculated already
                        sens->removeOverrideNames(
                            (dynamic_cast<const Sensitivity*>
                             (*iter))->getPacketName(), results);
                        // then see if there's any left do
                        if (sens->hasOverrideNames()){
                            nextDaySensToDo->push_back(sens);
                            nextDayGreekToDo.push_back(*iter);
                        }
                    }
                    iter = nextDayGreek.erase(iter);
                } else {
                    ++iter;
                }
            }
            // have we got anything to do?
            if (theta.get() || !nextDaySensToDo->empty()){
                // now set things up to do price + greeks on rolled date
                OutputRequestArraySP outRequest = 
                    !theta? OutputRequestArraySP(): theta->outRequests;
                // the do actual shifting and pricing
                ResultsSP shiftedResults;
                try {
                    shiftedResults = 
                        rollData->shiftAndCalculate(tweakGroup, 
                                                    nextDaySensToDo,
                                                    outRequest);
                } catch (exception& e){
                    untweakable = UntweakableSP(new Untweakable(e));
                }
                if (theta.get()){
                    // write theta results back
                    theta->writeResults(results, untweakable, 
                                        shiftedResults.get());
                    theta.reset(); // don't do theta again
                }
                // then results for extra greeks
                for (unsigned int i = 0; i < nextDayGreekToDo.size(); i++){
                    nextDayGreekToDo[i]->writeResults((*nextDaySensToDo)[i],
                                                      results, untweakable,
                                                      shiftedResults.get());
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e, "Theta::calculate", "Failed when calculating "
                             "'next day' sensitivities");
    }
}

/** Returns true if theta is the same as this */
bool Theta::equals(const Theta* theta) const{
    // default implementation
    return (theta && 
            theta->getClass() == getClass() && theta->offset == offset && 
            theta->hols->equals(hols.get()));
}

/** populate from market cache */
void Theta::getMarket(const IModel* model, const MarketData* market) {
    hols.getData(model, market);
}

/** Returns the value date (as reported by the instrument) before the
    theta shift. This method should only be called by shift(Theta*)
    methods */
const DateTime& Theta::Util::getOriginalValueDate() const{
    return origValueDate;
}

/** Returns the shifted value date. Equivalent to
    rollDate(getOriginalValueDate()). This method should only be
    called by shift(Theta*) methods */
const DateTime& Theta::Util::getNewValueDate() const{
    return newValueDate;
}

/** Returns true if assets forward value should be used (ie theta
    forward spot) */
bool Theta::Util::useAssetFwds() const{
    return theta->useAssetFwds();
}

//// returns original theta object
const Theta* Theta::Util::getTheta() const{
    return theta;
}

/** scale the results in the Results Object for this sensitivity
    by supplied factor */
void Theta::scaleResult(Results*     results,     // (M)
                        double       scaleFactor) const{
    const string& packetName = getPacketName();
    const string& sensName = getSensOutputName();
    if (packetName == sensName){
        // scale all results in packetName
        results->scale(packetName, scaleFactor);
    } else {
        // scale result sensName in packet packetName
        results->scale(packetName, sensName, scaleFactor);
    }
    
    // scale the output requests
    for (int i = 0; outRequests.get() &&  i < outRequests->size();++i) {
        const string& packetName = getPacketName();
        const string& sensName   = (*outRequests)[i]->getRequestName() +
            "_" +getSensOutputName();
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
void  Theta::addResult(Results*           results,     // (M)
                       const Results*     resultsToAdd,
                       double             scaleFactor) const {
    const string& packetName = getPacketName();
    const string& sensName = getSensOutputName();
    if (packetName == sensName){
        // add all results in packetName
        results->add(packetName, resultsToAdd, scaleFactor);
    } else {
        // add result sensName in packet packetName
        results->add(packetName, sensName, resultsToAdd, scaleFactor);
    }
    
    // scale the output requests
    for (int i=0; outRequests.get() && i < outRequests->size();++i) {
        const string& packetName = getPacketName();
        const string& sensName   = (*outRequests)[i]->getRequestName() + "_THETA";
        if (packetName == sensName){
            // add all results in packetName
            results->add(packetName, resultsToAdd, scaleFactor);
        } else {
            // add result sensName in packet packetName
            results->add(packetName, sensName, resultsToAdd, scaleFactor);
        }
    }
}


Theta::Util::Util(const Theta*    theta, 
                              const DateTime& origValueDate):
    theta(theta),
    origValueDate(origValueDate), 
    newValueDate(theta->rollDate(origValueDate)){}
                                               
/** Returns a Util object corresponding to this shift and
    original value date */
Theta::Util Theta::getUtil(const DateTime& origValueDate){
    return Util(this, origValueDate);
}

class ThetaHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default */
    class Factory: public SensitivityFactory::IDefault {
    public:
        virtual Sensitivity* createDefault(){
            HolidaySP hols(Holiday::weekendsOnly());
            return new Theta(Theta::DEFAULT_SHIFT, hols);
        }
    };
public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Theta, clazz);
        SUPERCLASS(SensControlAllNames);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultTheta);
        FIELD(offset, "number of days to roll");
        FIELD(hols, "holidays");
        FIELD(outRequests, "output requests for additional thetas");
        FIELD_MAKE_OPTIONAL(outRequests);
        // register how to build our sensitivity
        SensitivityFactory::addSens(Theta::NAME, 
                                    new Factory(), 
                                    new Theta(),
                                    Theta::Shift::TYPE);
    }

    static IObject* defaultTheta(){
        return new Theta();
    }
};

CClassConstSP const Theta::TYPE = CClass::registerClassLoadMethod(
    "Theta", typeid(Theta), ThetaHelper::load);

CClassConstSP const Theta::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "Theta::Shift", typeid(Theta::Shift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(Theta::RestorableShift, clazz);
    EXTENDS(Theta::Shift);
}
    
CClassConstSP const Theta::RestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "Theta::RestorableShift", typeid(Theta::RestorableShift), 
    restorableShiftLoad);

static void loadINextDaySensitivity(CClassSP& clazz){
    REGISTER_INTERFACE(IObject, clazz);
    EXTENDS(IObject);
}    
CClassConstSP const Theta::INextDaySensitivity::TYPE = 
CClass::registerInterfaceLoadMethod(
    "Theta::INextDaySensitivity", typeid(Theta::INextDaySensitivity), 
    loadINextDaySensitivity);
DRLIB_END_NAMESPACE
