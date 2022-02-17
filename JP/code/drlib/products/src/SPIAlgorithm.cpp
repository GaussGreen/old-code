//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIAlgorithm.cpp
//
//   Description : Algorithm interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPIAlgorithm.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

// this is the external interface for abstraction so that users can bolt in 
// any algorithm type they want - note this includes the SPIAlgorithmWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface IAlgoritmSPI as soon as possible
void IAlgorithmSPIInterface::load(CClassSP& clazz) {
    REGISTER_INTERFACE(IAlgorithmSPIInterface, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}

CClassConstSP const IAlgorithmSPIInterface::TYPE = CClass::registerInterfaceLoadMethod(
    "IAlgorithmSPIInterface", typeid(IAlgorithmSPIInterface), IAlgorithmSPIInterface::load);

/*****************************************************************************/

IAlgorithmSPI::~IAlgorithmSPI() {};

SPIAlgorithmStd::SPIAlgorithmStd(): SPIInterfaceIMS(TYPE),
        overInvestBound(0.0), underInvestBound(0.0),
        equityExposureMin(0.0), equityExposureMax(0.0),
        bondCrashSize(0.0),
        numAlgAssets(0), crashDenom(0.0),
        Cmax(0.0), Cmin(0.0),
        OIB(0.0), UIB(0.0) {} // for reflection

IAlgorithmSPISP SPIAlgorithmStd::getAlgorithmSPI() {
    static const string routine = "SPIAlgorithmStd::getAlgorithmSPI";
    // essentially just returns itself
    IAlgorithmSPISP theAlgorithm = IAlgorithmSPISP(this, NullDeleter()); // FIXME can be a problem
    return theAlgorithm;
}

void SPIAlgorithmStd::validatePop2Object(){
    static const string method = "SPIAlgorithmStd::validatePop2Object";
    int i;

    for(i=0; i<crashSizes.size(); i++) {
        if (!Maths::isPositive(crashSizes[i]) ||
            !Maths::isNegative(crashSizes[i]-1.0)) {
            isOK = false;
            err = "Crash size for equity #" + Format::toString(i+1) +
                " is " + Format::toString(crashSizes[i]) +
                " but must be strictly between 0% and 100%";
            return;
        }
    }

    if (Maths::isNegative(bondCrashSize)) {
        isOK = false;
        err = "Bond crash size (" + Format::toString(bondCrashSize) +
            ") must be non-negative";
        return;
    }

    numAlgAssets = exposureUpperBound.size();
    if (crashSizes.size() != numAlgAssets) {
        isOK = false;
        err = "Require equal number of crashSizes (" +
            Format::toString(crashSizes.size()) + ") and exposureUpperBounds (" +
            Format::toString(exposureUpperBound.size()) + ")";
            return;
    }
    for(i=0; i<numAlgAssets; i++) {
        if (Maths::isNegative(exposureUpperBound[i]) ||
            Maths::isPositive(exposureUpperBound[i]-1.0)) {
            isOK = false;
            err = "ExposureUpperBound for equity #" + Format::toString(i+1) +
                " is " + Format::toString(exposureUpperBound[i]) +
                " but must be between 0% and 100%";
            return;
        }
    }
    if (numAlgAssets==1) {
        // XXX for this to make sense shouldn't there be a check that exposureUpperBound[0]==100%?
        Cmin = crashSizes[0];
        Cmax = Cmin;
    } else if (numAlgAssets==2) {
        double UXA = exposureUpperBound[0];
        double UXB = exposureUpperBound[1];
        Cmin = UXB * crashSizes[1] + (1. - UXB) * crashSizes[0];
        Cmax = UXA * crashSizes[0] + (1. - UXA) * crashSizes[1];
        if (Maths::isPositive(Cmin-Cmax)) {
            err = "Cmax (" + Format::toString(Cmax) +
                ") must not be less than Cmin (" +
                Format::toString(Cmin) + ")";
            return;
        }
    } else {
        // caught elsewhere : only support 1 or 2 assets
    }

    if (Maths::isNegative(overInvestBound) ||
        Maths::isPositive(overInvestBound-1.0)) {
        isOK = false;
        err = "OverInvestBound is " + Format::toString(overInvestBound) +
            " but must be between 0% and 100%";
        return;
    }
    if (Maths::isNegative(underInvestBound) ||
        Maths::isPositive(underInvestBound-1.0)) {
        isOK = false;
        err = "UnderInvestBound is " + Format::toString(underInvestBound) +
            " but must be between 0% and 100%";
        return;
    }
    if (Maths::isNegative(equityExposureMin)) {
        isOK = false;
        err = "EquityExposureMin is " + Format::toString(equityExposureMin) +
            " but must not be negative.";
        return;
    }
    if (Maths::isPositive(equityExposureMin-equityExposureMax)) {
        isOK = false;
        err = "EquityExposureMin (" + Format::toString(equityExposureMin) +
            ") must not be greater than EquityExposureMax (" +
            Format::toString(equityExposureMax) + ")";
        return;
    }

    pE = DoubleArray(numAlgAssets);
    OIB = 1.+overInvestBound;
    UIB = 1.-underInvestBound;
}

void SPIAlgorithmStd::init(const DoubleArray* bondPrices) {
    static const string method = "SPIAlgorithmStd::init";
    bondCrashAdj = DoubleArray(bondPrices->size(), 0.0);
    for (int i=0; i<bondCrashAdj.size(); i++) {
        // if bondPrices[i] = 0 this may cause -#INF as a result from pow
        // this can only happen beyond bond maturity (o/w getBondPrices will parp)
        // we don't care about sc then
        if (Maths::isPositive((*bondPrices)[i])) {
            bondCrashAdj[i] = pow((*bondPrices)[i], -bondCrashSize) - 1.0;
        }
    }
}

bool SPIAlgorithmStd::doRebalance(double UL,
                         double SL,
                         int    iStep) {
    double minLoss = equityExposureMin*(Cmin + bondCrashAdj[iStep]);
    double maxLoss = equityExposureMax*(Cmax + bondCrashAdj[iStep]);
    bool doRebal = (((UL>Maths::min(SL, maxLoss)*OIB) && (UL>minLoss)) ||
                    ((UL<SL*UIB) && (UL<maxLoss)) ||
                    (Maths::isZero(SL) && (UL>minLoss)));
    return doRebal;
}

const DoubleArray& SPIAlgorithmStd::rebalWeights(double crash) {
    // Only defined for up to 2 assets
    if (numAlgAssets==1) {
        pE[0] = 1.0;
    } else {
        pE[1] = (crashSizes[0] - crash) / crashDenom;
        pE[0] = 1.0 - pE[1];
    }
    return pE;
}

double SPIAlgorithmStd::unbalCrash(const DoubleArray& nE,
                          const DoubleArray& E,
                          int                iStep) const {
    double e = 0.;
    double ec = 0.;
    for(int iAsset=0; iAsset<numAlgAssets; iAsset++) {
        double ei = nE[iAsset] * E[iAsset];
        e += ei;
        ec += ei * (crashSizes[iAsset] + bondCrashAdj[iStep]);
    }
    // XXX Note this is not nicely defined when e==0. Could rephrase spec to use "Loss" and not "Crash"
    // XXX numbers. Consider for future....
    double UC = Maths::isZero(e)? (Cmin + bondCrashAdj[iStep]) : (ec / e);
    return UC;
}

double SPIAlgorithmStd::sustCrash(double buffer, int iStep) const {
    return Maths::min(Cmax + bondCrashAdj[iStep],
                      Maths::max(Cmin + bondCrashAdj[iStep], buffer));
}

double SPIAlgorithmStd::targetExp(double SE) const {
    return Maths::min(Maths::max(SE, equityExposureMin), equityExposureMax);
}

double SPIAlgorithmStd::feeAtMin(IFeesSPIConstSP fees) const {
    // a rather high interdependence here ...
    return fees->getFeeAtMin(exposureUpperBound, equityExposureMin);
}

// incremental gap risk for this step
double SPIAlgorithmStd::gapRiskAtStep(double B,
                             double BF,
                             double equityComp,
                             int    iStep) const {
    double minLoss = equityExposureMin*(Cmin+bondCrashAdj[iStep]);
    double Bmin = BF / (1.0 - minLoss);
    double B_Bmin = B - Bmin;
    double G = 0.0;
    if (Maths::isPositive(B_Bmin)) {
        double gapNotional = equityComp-equityExposureMin*Bmin;
        G = gapNotional * gapNotional / (4.0 * B_Bmin);
    }

    return G;
}

void SPIAlgorithmStd::crossValidate(ICutoffSPISP cutoff, IFeesSPISP fees,
                                    bool isRainbowSPI) {
    if (isRainbowSPI && crashSizes.size() > 1) {
        isOK = false;
        err = "Can't have algorithm parameters (crashSize, exposureUpperBound)"
            "for more than one SPI asset in a Rainbow SPI";
        return;
    }

    // this validation has to be here in case it's a rainbow SPI
    // we might have to change this as rainbow SPI changes over time
    if (crashSizes.size()==1) {
        crashDenom = 1.0;
    } else if (crashSizes.size()==2) {
        crashDenom = crashSizes[0] - crashSizes[1];
        if (!Maths::isPositive(crashDenom)) {
            isOK = false;
            err = "Require crash size for equity #1 (" +
                Format::toString(crashSizes[0]) +
                ") to be greater than crash size for equity #2 (" +
                Format::toString(crashSizes[1]) + ")";
            return;
        }
    } else {
        isOK = false;
        err = "SPIAlgorithmStd can be defined only for 1 or 2 assets";
        return;
        /*throw ModelException(method,
          "SPIAlgorithmStd can be defined only for 1 or 2 assets");*/
    }

    // the proper cross validation
    try {
        cutoff->crossValidate(equityExposureMin, crashSizes);
        fees->crossValidate(equityExposureMin);
    }
    catch (...) {
        // dummy catch is to try to prevent vc6 from crashing owing to
        // misoptimisation
        throw;
    }
}

int SPIAlgorithmStd::getNumRiskyAssets() const {
    return numAlgAssets;
}

IObject* SPIAlgorithmStd::defaultSPIAlgorithmStd(){
    return new SPIAlgorithmStd();
}

/** Invoked when Class is 'loaded' */
void SPIAlgorithmStd::load(CClassSP& clazz){
    REGISTER(SPIAlgorithmStd, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(IAlgorithmSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPIAlgorithmStd);
    FIELD(exposureUpperBound,   "exposureUpperBound");
    FIELD(crashSizes,           "crashSizes");
    FIELD(overInvestBound,      "overInvestBound");
    FIELD(underInvestBound,     "underInvestBound");
    FIELD(equityExposureMin,    "equityExposureMin");
    FIELD(equityExposureMax,    "equityExposureMax");
    FIELD(bondCrashSize,        "bondCrashSize");
    FIELD_MAKE_OPTIONAL(bondCrashSize);
    FIELD(crashDenom,           "crashDenom");
    FIELD_MAKE_TRANSIENT(crashDenom);
    FIELD(numAlgAssets,            "numAlgAssets");
    FIELD_MAKE_TRANSIENT(numAlgAssets);
    FIELD(Cmax,                 "Cmax");
    FIELD_MAKE_TRANSIENT(Cmax);
    FIELD(Cmin,                 "Cmin");
    FIELD_MAKE_TRANSIENT(Cmin);
    FIELD(OIB,                  "OIB");
    FIELD_MAKE_TRANSIENT(OIB);
    FIELD(UIB,                  "UIB");
    FIELD_MAKE_TRANSIENT(UIB);
    FIELD(pE,                   "pE");
    FIELD_MAKE_TRANSIENT(pE);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(exposureUpperBound);
    FIELD_MAKE_OPTIONAL(crashSizes);
    FIELD_MAKE_OPTIONAL(overInvestBound);
    FIELD_MAKE_OPTIONAL(underInvestBound);
    FIELD_MAKE_OPTIONAL(equityExposureMin);
    FIELD_MAKE_OPTIONAL(equityExposureMax);
    clazz->setPublic(); // make visible to EAS/spreadsheet
/*        Addin::registerConstructor("SPI_ALGORITHM_STD",
                               Addin::RISK,
                               "Creates an SPIAlgorithmStd instance",
                               TYPE);*/
}
CClassConstSP const SPIAlgorithmStd::TYPE = CClass::registerClassLoadMethod(
    "SPIAlgorithmStd", typeid(SPIAlgorithmStd), SPIAlgorithmStd::load);


IAlgorithmSPISP SPIAlgorithmWrapper::getAlgorithmSPI() {
    static const string routine = "SPIAlgorithmWrapper::getAlgorithmSPI";
    // possibly some choices later
    if (!algorithmStd->isValid()) {
        throw ModelException(routine, "Invalid SPIAlgorithm (type " + SPIAlgorithmType +
                             ") : " + algorithmStd->errString());
    }
    IAlgorithmSPISP theAlgo = IAlgorithmSPISP(algorithmStd.get(), NullDeleter()); // FIXME

    return theAlgo;
}

// validation
void SPIAlgorithmWrapper::validatePop2Object(){
    static const string routine = "SPIAlgorithmWrapper::validatePop2Object";

    if (SPIAlgorithmType.empty()){
        throw ModelException(routine, "Blank SPIAlgorithmType specified!");
    }
    if (SPIAlgorithmType==SPI_ALGORITHM_TYPE_STD) {
        if (!algorithmStd.get()) {
            throw ModelException(routine, "Expected algorithmStd but none supplied!");
        }
    } else {
        throw ModelException(routine, "Unrecognised SPIAlgorithmType " + SPIAlgorithmType +
                             ". Expected " + SPI_ALGORITHM_TYPE_STD);
    }
}

/** Invoked when Class is 'loaded' */
void SPIAlgorithmWrapper::load(CClassSP& clazz){
    REGISTER(SPIAlgorithmWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IAlgorithmSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPIAlgorithmWrapper);
    FIELD(SPIAlgorithmType, "SPIAlgorithmType");
    FIELD(algorithmStd,  "algorithmStd");
    FIELD_MAKE_OPTIONAL(algorithmStd);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

// for reflection
SPIAlgorithmWrapper::SPIAlgorithmWrapper(): CObject(TYPE){}

IObject* SPIAlgorithmWrapper::defaultSPIAlgorithmWrapper(){
    return new SPIAlgorithmWrapper();
}
CClassConstSP const SPIAlgorithmWrapper::TYPE = CClass::registerClassLoadMethod(
    "SPIAlgorithmWrapper", typeid(SPIAlgorithmWrapper), SPIAlgorithmWrapper::load);

DRLIB_END_NAMESPACE
