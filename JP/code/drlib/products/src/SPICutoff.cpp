//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPICutoff.cpp
//
//   Description : Cutoff interface for SPI products
//
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPICutoff.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/
// run time post cutoff exposure behaviour
ISPIPostCutoffRT::~ISPIPostCutoffRT() {}

SPIZeroCutoff::SPIZeroCutoff() {}

double SPIZeroCutoff::getPostCutoffExposure(double B, double BF, int iStep) const {
    return 0.0;
}

double SPIZeroCutoff::getPostCutoffCrash() const {
    return 0.0;
}

SPIBMinusBFCutoff::SPIBMinusBFCutoff(double myCrash, int myFinalRebal) {
    crash = myCrash;
    finalRebal = myFinalRebal;
}

double SPIBMinusBFCutoff::getPostCutoffExposure(double B, double BF, int iStep) const {
    return iStep < finalRebal ?
        (Maths::isPositive(B) ? Maths::max(1.0-BF/B, 0.0) : 0.0) :
        0.0;
}

double SPIBMinusBFCutoff::getPostCutoffCrash() const {
    return crash;
}

/*****************************************************************************/
// this is the external interface for cutoff so that users can bolt in 
// any cutoff type they want - note this includes the SPICutoffWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface ICutoffSPI as soon as possible
void ICutoffSPIInterface::load(CClassSP& clazz) {
    REGISTER_INTERFACE(ICutoffSPIInterface, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}
CClassConstSP const ICutoffSPIInterface::TYPE = CClass::registerInterfaceLoadMethod(
    "ICutoffSPIInterface", typeid(ICutoffSPIInterface), ICutoffSPIInterface::load);

/*****************************************************************************/
ICutoffSPI::~ICutoffSPI() {}

/*****************************************************************************/
// types of exposure behaviour post cutoff
const string SPIPostCutoff::ZERO = "Zero";
const string SPIPostCutoff::B_MINUS_BF = "BMinusBF";

SPIPostCutoff::SPIPostCutoff() : SPIInterfaceIMS(TYPE),
                numNonZero(0), crash(0.0) {} // for reflection

SPIPostCutoff::SPIPostCutoff(CClassConstSP clazz) : SPIInterfaceIMS(clazz),
                numNonZero(0), crash(0.0) {}

// foul and nasty becuase of the presence of the 2-factor algorithm
// cutoff either kills exposure entirely or leaves it at B-BF
// If latter and there are 2 assets then either asset can be zeroed
// Interacts heavily with SPIAlgorithmStd::rebalWeights which needs a
// crash to generate weights as follows:
// w[0] = 1 - w[1]
// w[1] = (crash[0] - crash) / (crash[0] - crash[1])
// For one asset we don't need a crash
void SPIPostCutoff::validatePostCutoff(DoubleArray crashSizes) {
    if (!postCutoffExposure.empty() &&
            crashSizes.size() != postCutoffExposure.size()) {
        throw ModelException("SPIPostCutoff::validatePostCutoff",
                "Must have post-cutoff exposure type for each algorithm asset");
    }
    if (!postCutoffExposure.empty()) {
        for (int i = 0; i < postCutoffExposure.size(); i++) {
            if (postCutoffExposure[i] == B_MINUS_BF) {
                numNonZero++;
                crash += crashSizes[i];
            }
        }
        if (numNonZero > 0) {
            crash /= numNonZero;
        }
    }
}

SPIPostCutoffRTSP SPIPostCutoff::getPostCutoff(int finalRebal) {
    if (numNonZero > 0) {
        SPIPostCutoffRTSP postCutSP = SPIPostCutoffRTSP(new SPIBMinusBFCutoff(crash,
                                                                finalRebal));
        return postCutSP;
    } else {
        SPIPostCutoffRTSP postCutSP = SPIPostCutoffRTSP(new SPIZeroCutoff());
        return postCutSP;
    }
}

bool SPIPostCutoff::isCutoff(double  UE,
                        double  TE,
                        double  B,
                        double  BF) const {
    throw ModelException("SPIPostCutoff::isCutoff", "Internal Error:"
                        "Should not end up in base class function");
}

void SPIPostCutoff::crossValidate(double equityExposureMin,
                                  DoubleArray crashSizes) {
    throw ModelException("SPIPostCutoff::crossValidate", "Internal Error:"
                        "Should not end up in base class function");
}

/** Invoked when Class is 'loaded' */
void SPIPostCutoff::load(CClassSP& clazz){
    REGISTER(SPIPostCutoff , clazz);
    SUPERCLASS(SPIInterfaceIMS);
    EMPTY_SHELL_METHOD(defaultSPIPostCutoff);
    FIELD(postCutoffExposure, "after cutoff what exposure remains");
    FIELD(numNonZero, "number of assets with non-zero exposure post cutoff");
    FIELD(crash, "to fix weights of each risky asset post cutoff");
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(postCutoffExposure);
    FIELD_MAKE_TRANSIENT(numNonZero);
    FIELD_MAKE_TRANSIENT(crash);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

IObject* SPIPostCutoff::defaultSPIPostCutoff (){
    return new SPIPostCutoff ();
}

CClassConstSP const SPIPostCutoff::TYPE = CClass::registerClassLoadMethod(
    "SPIPostCutoff", typeid(SPIPostCutoff), SPIPostCutoff::load);

/*****************************************************************************/

SPICutoffStd::SPICutoffStd (): SPIPostCutoff(TYPE),
                        equityExposureCutoff(0.0) {} // for reflection

bool SPICutoffStd::isCutoff(double  UE,
              double  TE,
              double  B,
              double  BF) const {
    bool isCutoffNow = (UE < equityExposureCutoff ||
                        TE < equityExposureCutoff);
    return isCutoffNow;
}

ICutoffSPISP SPICutoffStd::getCutoffSPI() {
    static const string routine = "SPICutoffStd::getCutofffSPI";
    // essentially just returns itself
    ICutoffSPISP theCutoff = ICutoffSPISP(this, NullDeleter()); // FIXME can be a problem
    return theCutoff;
}

void SPICutoffStd::crossValidate(double equityExposureMin,
                                 DoubleArray crashSizes) {
    if (!Maths::isZero(equityExposureMin*equityExposureCutoff)) {
        throw ModelException("SPICutoffStd::crossValidate", "Either equityExposureMin or equityExposureCutoff must be 0.");
    }
    validatePostCutoff(crashSizes);
}

IObject* SPICutoffStd::defaultSPICutoffStd (){
    return new SPICutoffStd ();
}

/** Invoked when Class is 'loaded' */
void SPICutoffStd::load(CClassSP& clazz){
    REGISTER(SPICutoffStd , clazz);
    SUPERCLASS(SPIPostCutoff);
    IMPLEMENTS(ICutoffSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPICutoffStd);
    FIELD(equityExposureCutoff,             "equityExposureCutoff");
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(equityExposureCutoff);
    clazz->setPublic(); // make visible to EAS/spreadsheet
/*        Addin::registerConstructor("SPI_CUTOFF_STD",
                               Addin::RISK,
                               "Creates an SPICutoffStd instance",
                               TYPE);*/
}
CClassConstSP const SPICutoffStd::TYPE = CClass::registerClassLoadMethod(
    "SPICutoffStd", typeid(SPICutoffStd), SPICutoffStd::load);

////////////////////////////////////////

SPICutoffNearBondFloor::SPICutoffNearBondFloor (): SPIPostCutoff(TYPE),
                        bondFloorBuffer(0.0) {} // for reflection

bool SPICutoffNearBondFloor::isCutoff(double  UE,
              double  TE,
              double  B,
              double  BF) const {
    bool isCutoffNow = (B - BF < bondFloorBuffer);
    return isCutoffNow;
}

ICutoffSPISP SPICutoffNearBondFloor::getCutoffSPI() {
    static const string routine = "SPICutoffNearBondFloor::getCutofffSPI";
    // essentially just returns itself
    ICutoffSPISP theCutoff = ICutoffSPISP(this, NullDeleter()); // FIXME can be a problem
    return theCutoff;
}

void SPICutoffNearBondFloor::crossValidate(double equityExposureMin,
                                           DoubleArray crashSizes) {
    if (!Maths::isZero(equityExposureMin*bondFloorBuffer)) {
        throw ModelException("SPICutoffNearBondFloor::crossValidate", "Either equityExposureMin or bondFloorBuffer must be 0.");
    }
    try {
        validatePostCutoff(crashSizes);
    }
    catch (exception& e) {
        // this catch is to try to prevent vc6 from crashing owing to
        // misoptimisation
        throw ModelException(e, "SPICutoffNearBondFloor::crossValidate()");
    }
}

IObject* SPICutoffNearBondFloor::defaultSPICutoffNearBondFloor (){
    return new SPICutoffNearBondFloor ();
}

/** Invoked when Class is 'loaded' */
void SPICutoffNearBondFloor::load(CClassSP& clazz){
    REGISTER(SPICutoffNearBondFloor , clazz);
    SUPERCLASS(SPIPostCutoff);
    IMPLEMENTS(ICutoffSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPICutoffNearBondFloor);
    FIELD(bondFloorBuffer, "B-BF<bondFloorBuffer means cutoff");
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(bondFloorBuffer);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPICutoffNearBondFloor::TYPE = CClass::registerClassLoadMethod(
    "SPICutoffNearBondFloor", typeid(SPICutoffNearBondFloor), SPICutoffNearBondFloor::load);

////////////////////////////////////////

ICutoffSPISP SPICutoffWrapper::getCutoffSPI() {
    static const string routine = "SPICutoffWrapper::getCutoffSPI";

    ICutoffSPISP theCutoff;
    if (SPICutoffType==SPI_CUTOFF_TYPE_STD) {
        theCutoff = ICutoffSPISP(cutoffStd.get(), NullDeleter()); //FIXME
    } else if (SPICutoffType==SPI_CUTOFF_TYPE_BF) {
        theCutoff = ICutoffSPISP(cutoffNearBF.get(), NullDeleter()); //FIXME
    } else {
        throw ModelException(routine, "Unrecognised cutoff type : " +
                             SPICutoffType);
    }

    // need to cast for validity check
    SPIInterfaceIMS* actualCutoff = dynamic_cast<SPIInterfaceIMS*>(theCutoff.get());
    if (!actualCutoff) {
        throw ModelException(routine, "Internal error in SPI cutoff validation");
    }
    if (!actualCutoff->isValid()) {
        throw ModelException(routine, "Invalid SPICutoff (type " + SPICutoffType +
                             ") : " + actualCutoff->errString());
    }
    return theCutoff;
}

// validation
void SPICutoffWrapper::validatePop2Object(){
    static const string routine = "SPICutoffWrapper::validatePop2Object";

    if (SPICutoffType.empty()){
        throw ModelException(routine, "Blank SPICutoffType specified!");
    }
    if (SPICutoffType==SPI_CUTOFF_TYPE_STD) {
        if (!cutoffStd.get()) {
            throw ModelException(routine, "Expected cutoffStd but none supplied!");
        }
    } else if (SPICutoffType==SPI_CUTOFF_TYPE_BF) {
        if (!cutoffNearBF.get()) {
            throw ModelException(routine, "Expected cutoffNearBF but none supplied!");
        }
    } else {
        throw ModelException(routine, "Unrecognised SPICutoffType " + SPICutoffType +
                             ". Expected " + SPI_CUTOFF_TYPE_STD + " or " + SPI_CUTOFF_TYPE_BF);
    }
}

/** Invoked when Class is 'loaded' */
void SPICutoffWrapper::load(CClassSP& clazz){
    REGISTER(SPICutoffWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICutoffSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPICutoffWrapper);
    FIELD(SPICutoffType, "SPICutoffType");
    FIELD(cutoffStd,  "cutoffStd");
    FIELD_MAKE_OPTIONAL(cutoffStd);
    FIELD(cutoffNearBF,  "cutoffNearBF");
    FIELD_MAKE_OPTIONAL(cutoffNearBF);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

// for reflection
SPICutoffWrapper::SPICutoffWrapper(): CObject(TYPE){}

IObject* SPICutoffWrapper::defaultSPICutoffWrapper(){
    return new SPICutoffWrapper();
}
CClassConstSP const SPICutoffWrapper::TYPE = CClass::registerClassLoadMethod(
    "SPICutoffWrapper", typeid(SPICutoffWrapper), SPICutoffWrapper::load);

DRLIB_END_NAMESPACE
