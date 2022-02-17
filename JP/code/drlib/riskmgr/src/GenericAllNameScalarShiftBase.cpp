/**
 * @file GenericAllNameScalarShift.hpp
 *
 *
 */

#include "edginc/config.hpp"
#include "edginc/GenericAllNameScalarShiftBase.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * A tweak applied to all names simultaneously (rather than to each in
 * turn, as is more usual). Note that this is a base non-templated class. It
 * contains the common code needed by GenericAllNameScalarShift
 *
 * This is not a particularly common requirement.  See
 * e.g. CCMAbsoluteBetaSens for an example instantiation.
 */

/** Name of packet for output ("Instrument") */
const string& GenericAllNameScalarShiftBase::getPacketName() const{
    return Results::INSTRUMENT_PACKET;
}
    
/** Units in which sensitivity outputs are expressed. Default
    implementation is 1.0 */
double GenericAllNameScalarShiftBase::getSensitivityUnit() const {
    return 1.;
}

GenericAllNameScalarShiftBase::GenericAllNameScalarShiftBase(
    CClassConstSP clazz,
    const string& outputName,
    double shiftSize):
    SensControlAllNames(clazz, outputName), shiftSize(shiftSize) {}

/**
 * Compute the greek
 *
 * Calculates a one-sided first derivative.
 */
void GenericAllNameScalarShiftBase::calculate(
    TweakGroup* tweakGroup, Results* results) {
    try {
        OutputNameConstSP sensName(new OutputName(getSensOutputName()));
        // if result already computed then skip
        if (results->exists(getPacketName(), sensName)) {
            return;
        }
        if (!SensMgr::dependsUpon(IObjectSP::attachToRef(tweakGroup),
                                  shiftInterface())){
            // if we don't have any applicable objects return N/A
            results->storeNotApplicable(this);
            return;
        }
        
        try {
            // DON'T use storeScalarGreek(..., this)
            results->storeScalarGreek(
                calcOneSidedFirstDeriv(tweakGroup, results),
                getPacketName(), sensName);
        }
        catch (exception& e) {
            // DON'T use storeGreek(..., this)
            results->storeGreek(IObjectSP(new Untweakable(e)),
                                getPacketName(), sensName);
        }
    }
    catch (exception& e) {
        throw ModelException(e, "GenericAllNameScalarShiftBase::calculate");
    }
}

/**
 * The e to use in calculating (f(a')-f(a)) / e
 *
 * Note that this incorporates SENSITIVITY_UNIT so that the output
 * sensitivities end up scaled into the desired units (e.g. per basis point).
 */
double GenericAllNameScalarShiftBase::divisor() const {
    try {
        if (Maths::isZero(shiftSize)) {
            throw ModelException("GenericAllNameScalarShiftBase::divisor",
                                 "Shift size is zero");
        }
        // We use this as a hook to apply getSensitivityUnit()
        return shiftSize / getSensitivityUnit();
    } 
    catch (ModelException& e) {
        throw ModelException(e, "GenericAllNameScalarShift::divisor");
    }
}

bool GenericAllNameScalarShiftBase::discreteShift() const{
    return false;
}

/** Returns size of the shift */
double GenericAllNameScalarShiftBase::getShiftSize() const {
    return shiftSize;
}

void GenericAllNameScalarShiftBase::load(CClassSP& clazz) {
    REGISTER(GenericAllNameScalarShiftBase, clazz);
    SUPERCLASS(SensControlAllNames);
    IMPLEMENTS(Additive);
    FIELD(shiftSize, "Shift size");
}

CClassConstSP const GenericAllNameScalarShiftBase::TYPE =
CClass::registerClassLoadMethod(
    "GenericAllNameScalarShiftBase", 
    typeid(GenericAllNameScalarShiftBase), load);

DRLIB_END_NAMESPACE
