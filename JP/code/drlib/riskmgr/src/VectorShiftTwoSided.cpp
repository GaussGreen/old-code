
#include "edginc/config.hpp"
#include "edginc/VectorShiftTwoSided.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/SensMgr.hpp"

DRLIB_BEGIN_NAMESPACE
VectorShiftTwoSided::~VectorShiftTwoSided(){}

/** Overridden to perform two sided tweak */
void VectorShiftTwoSided::calculate(TweakGroup*  tweakGroup,
                                    Results*     results){
    static const string method("VectorShiftTwoSided::calculate");
    const string& secondOrderName = getSecondOrderSensOutputName();
    // get list of names to calculate result for
    OutputNameArrayConstSP names(this->names(tweakGroup));
        
    if (hasOverrideNames()){
        // manually cope with any explicit names that aren't present.
        // This should be a virtual method on SensControl/Sensitivity - 
        // also see comments in SensMgr::names(SensControl*, Results*)
        OutputNameArrayConstSP allNames(SensMgr(tweakGroup).allNames(this));
        OutputNameArraySP extraNames(OutputName::difference(names,
                                                            allNames));
        for (int i = 0; i < extraNames->size(); i++){
            results->storeGreek(IObjectSP(new NotApplicable()),
                                getPacketName(),
                                (*extraNames)[i]);
            results->storeGreek(IObjectSP(new NotApplicable()),
                                secondOrderName,
                                (*extraNames)[i]);
        }
    }

    if (names->empty()) {
        results->storeNotApplicable(this);
        results->storeNotApplicable(secondOrderName);
        return; // all done
    } 
    Instrument* inst = tweakGroup->getInstrument();
        
    // see if the instrument has a last sens date method
    LastSensDate* lsd = dynamic_cast<LastSensDate*>(inst);
    DateTime      endDate;
    DateTime      valueDate;
    if (lsd) {
        valueDate = inst->getValueDate();
    }
        
    for (int idx = 0; idx < names->size(); idx++){
        OutputNameConstSP name = (*names)[idx];
        // store the name of what we want to shift
        setMarketDataName(name);
        /* skip over where result has been calculated already */
        if (!results->exists(this)){
            double origShiftSize = getShiftSize(); // save this
            try {
                // get end date for this name/greek
                if (lsd) {
                    endDate   = lsd->endDate(this);
                }
        
                // get expiries for which we should tweak
                cachedExpiries = getExpiries(tweakGroup);
                int numExpiries = cachedExpiries->size();
                // get max tweak sizes allowed
                DoubleArrayConstSP tweakShiftSizes = 
                    calculateTweakSizes(tweakGroup, name, *cachedExpiries);
                // create room for storing the results
                ExpiryResultArraySP firstDerivs(
                    new ExpiryResultArray(numExpiries));
                ExpiryResultArraySP secondDerivs(
                    new ExpiryResultArray(numExpiries));
                // then loop over the expiries
                bool expired = false;
                for (int j = 0; j < numExpiries; j++){
                    // store the expiry/shiftSize which we want to tweak
                    setExpiry((*cachedExpiries)[j]);
                    // calculate sens (if not expired)
                    pair<double, double> derivs = expired ?
                        make_pair(0.0, 0.0) :
                        twoSidedDerivs(tweakGroup, results);
                    (*firstDerivs)[j] = ExpiryResult(expiry, derivs.first);
                    (*secondDerivs)[j] = ExpiryResult(expiry, derivs.second);
                    // do we need to tweak anymore ?
                    if (lsd) {
                        expired = expiry->
                            toDate(valueDate).isGreater(endDate);
                    }
                }
                setShiftSize(origShiftSize); // restore value
                // and store it
                results->storeGreek(firstDerivs, this);
                results->storeGreek(secondDerivs, secondOrderName, name);
            } catch (exception& e) {
                setShiftSize(origShiftSize); // restore value
                IObjectSP untweakable(new Untweakable(e));
                results->storeGreek(untweakable, this);
                results->storeGreek(untweakable, secondOrderName, name);
            }
        }
    }
}


/** Overridden to scale 1st and 2nd order derivatives */
void VectorShiftTwoSided::scaleResult(Results*     results,
                                      double       scaleFactor) const{
    VectorShift::scaleResult(results, scaleFactor); // do 1st order derive
    const string& secondOrderName = getSecondOrderSensOutputName();
    if (getPacketName() == getSensOutputName()) {
        results->scale(secondOrderName, scaleFactor);
    } else {
        results->scale(getPacketName(), secondOrderName, scaleFactor);
    }
}

/** Overridden to add 1st and 2nd order derivatives */
void VectorShiftTwoSided::addResult(Results*           results,     // (M)
                                    const Results*     resultsToAdd,
                                    double             scaleFactor) const{
    VectorShift::addResult(results, resultsToAdd,
                           scaleFactor); // do 1st order derive
    const string& secondOrderName = getSecondOrderSensOutputName();
    if (getPacketName() == getSensOutputName()) {
        results->add(secondOrderName, resultsToAdd, scaleFactor);
    } else {
        results->add(getPacketName(), secondOrderName, resultsToAdd,
                     scaleFactor);
    }
}

/** Returns the shift(s) which have been made for the current pricing
    call */
ScalarShiftArray VectorShiftTwoSided::getComponentShifts() const{
    // need array of ScalarShiftConstSP to make this work - but this breaks
    // the array template - would need const array template...
    VectorShiftTwoSided* shift = const_cast<VectorShiftTwoSided*>(this);
    return ScalarShiftArray(1, ScalarShiftSP::attachToRef(shift));
}

/** Note VectorShiftTwoSided is abstract. Create a scalar shift of
    type clazz, which uses outputName (eg VEGA_PARALLEL) to
    identify results and with given shiftSize */
VectorShiftTwoSided::VectorShiftTwoSided(CClassConstSP clazz,
                                         const string& outputName,
                                         double        shiftSize):
    VectorShift(clazz, outputName, shiftSize){}

/** for reflection */
VectorShiftTwoSided::VectorShiftTwoSided(CClassConstSP clazz,
                                         const string& outputName):
    VectorShift(clazz, outputName){}

void VectorShiftTwoSided::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VectorShiftTwoSided, clazz);
    SUPERCLASS(VectorShift);
}

CClassConstSP const VectorShiftTwoSided::TYPE = CClass::registerClassLoadMethod(
    "VectorShiftTwoSided", typeid(VectorShiftTwoSided), load);

DRLIB_END_NAMESPACE

