//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : QuasiContractualBaseCorrelation.hpp
//
//   Description : Composite greek: Applies a multiple scenario as defined in
//                 this class, and then computes the sensitivity returned by
//                 derived classes through the "sensitivityToComputeAfterShift"
//                 virtual method (typically BetaSkewParallel/Matrix, used
//                 to compute the base correlation sensitivities of the 
//                 "quasi-contractual" capital structure representation of a
//                 trade).
//
//   Author      : Jose Hilera
//
//   Date        : April 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/QuasiContractualBaseCorrelation.hpp"

DRLIB_BEGIN_NAMESPACE

#define DEFAULT_CREDIT_SPREAD_LEVEL 0.0060 /* 60 bps */
#define DEFAULT_SKEW_LEVEL          0.8    /* 80%*/

QuasiContractualBaseCorrelation::QuasiContractualBaseCorrelation(
        CClassConstSP clazz, 
        const string& name) : 
    SensControlAllNames(clazz, name), 
    creditSpreadsLevel(DEFAULT_CREDIT_SPREAD_LEVEL), 
    compressionRatio(0.0), 
    piecewiseMappingFunctionLevel(0.0), 
    skewLevel(DEFAULT_SKEW_LEVEL)
{}


/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP QuasiContractualBaseCorrelation::shiftInterface() const {
    return IShift::TYPE;
}


/** Shifts the object (which supports being set by this scenario) using given 
    shift. The return value indicates whether or not components of this object 
    need to be tweaked */
bool QuasiContractualBaseCorrelation::shift(IObjectSP obj) {
    return dynamic_cast<IShift &>(*obj).sensShift(this);
}


/** Returns the level to set the credit spreads to */
double QuasiContractualBaseCorrelation::getCreditSpreadsLevel() const {
    return creditSpreadsLevel;
}


/** Returns the compression ratio to use when rescaling the historical betas
 * dispersion. */
double QuasiContractualBaseCorrelation::getCompressionRatio() const {
    return compressionRatio;
}


/** Returns the value to set the yPoints of the PiecewiseMappingFunction to */
double QuasiContractualBaseCorrelation::getPiecewiseMappingFunctionLevel() const {
    return piecewiseMappingFunctionLevel;
}


/** Returns the value to set the skews in the skewSurface to */
double QuasiContractualBaseCorrelation::getSkewLevel() const {
    return skewLevel;
}


/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns true */
bool QuasiContractualBaseCorrelation::discreteShift() const{
    return true;
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double QuasiContractualBaseCorrelation::divisor() const{
    return 1.0;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP QuasiContractualBaseCorrelation::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}

void QuasiContractualBaseCorrelation::restore(IObjectSP obj) {
    // cast obj to IShift and then invoke restore method
    IRestorableShift& quasiObj = dynamic_cast<IRestorableShift&>(*obj);
    quasiObj.sensRestore(this);
}


/** This calculate method applies two tweaks and returns the difference in value
 * between tweaks.
 * The first tweak is "this", ie, RecoveryTweakWithShift, relative or absolute.
 * The second tweak is a RecoveryTweak, which is created in this method, stored
 * in a SensitivityArray and calculated within the shiftAndCalculate method */
void QuasiContractualBaseCorrelation::calculate(TweakGroup* tweakGroup, 
                                                CResults*   results) 
{
    static const string method("RecoveryTweakWithShift::calculate");
     
    try {
        UntweakableSP untweakable;

        // Get the sensitivity to compute after the initial shift. The actual
        // sensitivity depens on the concrete implementations of this class
        SensitivitySP sensSP = sensitivityToComputeAfterShift();
        // Even though the scenario part of this sensitivity is applied to all 
        // names (this is a SensControlAllNames greek), the "toTweak" parameter 
        // can be populated: it indicates the names that the sensitivity 
        // computed after the scenario should be calculated for.
        if (hasOverrideNames()) {
            OutputNameArrayConstSP names = overrideNames();
            sensSP->storeOverrideNames(OutputNameArraySP(names.clone()));    
        }
        SensitivityArraySP sensToDoAfterShift(new SensitivityArray(1, sensSP));

        ResultsSP shiftedResults;
        try {
            // All the magic happens within shiftAndCalculate
            shiftedResults = 
                this->shiftAndCalculate(tweakGroup,
                                        sensToDoAfterShift,
                                        OutputRequestArrayConstSP());

            // We have the results.
            // However, the results are under the internal sensitivity's packet of the
            // "shiftedResults" - Need to copy them into this greek's packet 
            // within "results"
            const vector<pair<OutputNameConstSP, IObjectConstSP> > packetResults = 
                shiftedResults->listPacketResults(sensSP->getSensOutputName());

            for (unsigned int i=0; i < packetResults.size(); i++) {
                results->storeGreek(IObjectSP::constCast(packetResults[i].second),
                                    getSensOutputName(), // The packet name is this greek's name
                                    packetResults[i].first);
            }
        }
        catch (exception& e) {
            OutputNameConstSP myName(new OutputName(getSensOutputName()));

            results->storeGreek(IObjectSP(new Untweakable(e)), 
                                getSensOutputName(),
                                myName);
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}


void QuasiContractualBaseCorrelation::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(QuasiContractualBaseCorrelation, clazz);
    SUPERCLASS(SensControlAllNames);
    IMPLEMENTS(Additive);
    
    FIELD(creditSpreadsLevel,
                 "Level to set the credit spreads to. Default: " +
                 Format::toString(DEFAULT_CREDIT_SPREAD_LEVEL) + " (" +
                 Format::toString(10000 * DEFAULT_CREDIT_SPREAD_LEVEL) + " bps)");
    FIELD_MAKE_OPTIONAL(creditSpreadsLevel);

    FIELD(compressionRatio,
                 "Compression ratio to use when rescaling historical "
                 "betas dispersion. Default: 0");
    FIELD_MAKE_OPTIONAL(compressionRatio);

    FIELD(piecewiseMappingFunctionLevel,
                 "Value to set the yPoints of the PiecewiseMappingFunction "
                 "to. Default: 0");
    FIELD_MAKE_OPTIONAL(piecewiseMappingFunctionLevel);

    FIELD(skewLevel,
                 "Value to set the skews in the skewSurface to. Default: " +
                 Format::toString(100 * DEFAULT_SKEW_LEVEL) + " %");
    FIELD_MAKE_OPTIONAL(skewLevel);
}


CClassConstSP const QuasiContractualBaseCorrelation::TYPE =
CClass::registerClassLoadMethod("QuasiContractualBaseCorrelation",
                                typeid(QuasiContractualBaseCorrelation),
                                QuasiContractualBaseCorrelation::load);

// IShift
CClassConstSP const QuasiContractualBaseCorrelation::IShift::TYPE =
CClass::registerInterfaceLoadMethod("QuasiContractualBaseCorrelation::IShift",
                                    typeid(QuasiContractualBaseCorrelation::IShift),
                                    0);

QuasiContractualBaseCorrelation::IShift::~IShift() 
{}

// IRestorableShift
static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(QuasiContractualBaseCorrelation::IRestorableShift, clazz);
    EXTENDS(QuasiContractualBaseCorrelation::IShift);
}

CClassConstSP const QuasiContractualBaseCorrelation::IRestorableShift::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "QuasiContractualBaseCorrelation::IRestorableShift", 
        typeid(QuasiContractualBaseCorrelation::IRestorableShift), 
        restorableShiftLoad);

QuasiContractualBaseCorrelation::IRestorableShift::~IRestorableShift() 
{}

DRLIB_END_NAMESPACE
