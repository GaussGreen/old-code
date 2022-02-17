//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : RecoveryTweakWithShift.cpp
//
//   Description : Applies two tweaks to shift/set the Recovery Rate:
//                 First shifts the recovery rate and then calculates the
//                 RecoveryTweak greek
//
//   Author      : Jose Hilera
//
//   Date        : 19 October 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/RecoveryTweak.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Format.hpp"
#include "edginc/Results.hpp"
#include "edginc/RecoveryTweakWithShift.hpp"

DRLIB_BEGIN_NAMESPACE

static const double INITIAL_SHIFT_SIZE = 0.10; /* 10% */

RecoveryTweakWithShift::~RecoveryTweakWithShift()
{}


RecoveryTweakWithShift::RecoveryTweakWithShift():
    RecoveryShiftBase(TYPE, NAME),
    initialShiftSize(INITIAL_SHIFT_SIZE),
    relativeShift(true)
{}


RecoveryTweakWithShift::RecoveryTweakWithShift(double shiftSize,
                                               bool   relativeShift):
    RecoveryShiftBase(TYPE, NAME, shiftSize), 
    initialShiftSize(INITIAL_SHIFT_SIZE),
    relativeShift(relativeShift)
{}


IObject* RecoveryTweakWithShift::defaultConstructor() {
    return new RecoveryTweakWithShift();
}


/** Returns the new recovery level given the original one */
double RecoveryTweakWithShift::applyShift(double unadjRecovery) {
    static const string method = "RecoveryTweakWithShift::applyShift";

    // make adjustment
    double shiftSize = getShiftSize();
    double recovery = relativeShift? unadjRecovery * (1.0 + shiftSize) :
                                     (unadjRecovery + shiftSize);

    if ((recovery < 0.0) || (recovery > 1.0)) {
         throw ModelException(method, 
                              "Shifted recovery is out of [0,1] (Recovery=" +
                              Format::toString(recovery) + ").");
    }
    return recovery;
}

/** Returns the original recovery level given the adjusted one */
double RecoveryTweakWithShift::undoShift(double adjRecovery){
    double shiftSize = getShiftSize();
    return relativeShift? adjRecovery / (1.0 + shiftSize) :
                          (adjRecovery - shiftSize);
}


/** This calculate method applies two tweaks and returns the difference in value
 * between tweaks.
 * The first tweak is "this", ie, RecoveryTweakWithShift, relative or absolute.
 * The second tweak is a RecoveryTweak, which is created in this method, stored
 * in a SensitivityArray and calculated within the shiftAndCalculate method */
void RecoveryTweakWithShift::calculate(TweakGroup* tweakGroup, CResults* results) {
    static const string method("RecoveryTweakWithShift::calculate");
     
    try {
        UntweakableSP untweakable;
        // copy this tweak group to avoid issues - is this really required?
        IObjectSP ioSP(tweakGroup->clone());
        TweakGroupSP tweakGroupSP = TweakGroupSP::dynamicCast(ioSP);

        // get list of names to calculate result for. Remove blanks
        OutputNameArrayConstSP names(this->names(tweakGroupSP.get(), results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        }
        else {
            // Create the RecoveryTweak - its shiftSize is the "shiftSize" 
            // parameter of this greek
            SensitivitySP sensSP(new RecoveryTweak(getShiftSize()));
            SensitivityArraySP sensToDoAfterShift(new SensitivityArray(1, sensSP));

            // The size of this tweak is the initialShiftSize
            setShiftSize(initialShiftSize);
            ResultsSP shiftedResults;
            for (int idx = 0; idx < names->size(); idx++){
                // store what we want to shift:
                // 1. In the initial shift
                setMarketDataName((*names)[idx]);
                
                // 2. In the greek calculation (we want just that same name)
                OutputNameArraySP outName(new OutputNameArray(1, (*names)[idx]));
                sensSP->storeOverrideNames(outName);

                try {
                    // All the magic happens within shiftAndCalculate
                    shiftedResults = 
                        this->shiftAndCalculate(tweakGroupSP.get(), 
                                                sensToDoAfterShift,
                                                OutputRequestArrayConstSP());

                    // We have the results.
                    // However, the results are under the RecoveryTweak packet of the
                    // "shiftedResults" - Need to copy them into this greek's packet 
                    // within "results"
                    const vector<pair<OutputNameConstSP, IObjectConstSP> > packetResults = 
                        shiftedResults->listPacketResults(RecoveryTweak::NAME);

                    for (unsigned int i=0; i < packetResults.size(); i++) {
                        results->storeGreek(IObjectSP::constCast(packetResults[i].second),
                                            NAME, // The packet name is this greek's name
                                            packetResults[i].first);
                    }
                }
                catch (exception& e) {
                    results->storeGreek(IObjectSP(new Untweakable(e)), 
                                        NAME,
                                        (*names)[idx]);
                }
            }
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Invoked when SmallClass is 'loaded' */
void RecoveryTweakWithShift::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(RecoveryTweakWithShift, clazz);
    SUPERCLASS(RecoveryShiftBase);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(relativeShift, "True: the initial shift is multiplicative "
                 "(ie, times [1+shiftSize]), else additive shift");
    FIELD_MAKE_OPTIONAL(relativeShift);
    FIELD(initialShiftSize, "Size of the initial shift to the par recovery.");
    FIELD_MAKE_OPTIONAL(initialShiftSize);
}


CClassConstSP const RecoveryTweakWithShift::TYPE = CClass::registerClassLoadMethod(
    "RecoveryTweakWithShift", typeid(RecoveryTweakWithShift), load);
const string RecoveryTweakWithShift::NAME = "RECOVERY_TWEAK_WITH_SHIFT";

/** Included in RiskMgrLib::linkInClasses() to force link to include this */
bool RecoveryTweakWithShiftLinkIn() {
    return RecoveryTweakWithShift::TYPE != NULL;
}


DRLIB_END_NAMESPACE
