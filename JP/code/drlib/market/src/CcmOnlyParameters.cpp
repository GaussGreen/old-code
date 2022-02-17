//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CcmOnlyParameters.cpp
//
//   Description : Per name additional parameters needed for CCM model.
//                 Based on (rather, copy of) the now deprecated CCMParameters
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_CCMONLYPARAMETERS_CPP
#include "edginc/CcmOnlyParameters.hpp"

DRLIB_BEGIN_NAMESPACE

#define CHECK_BOUNDS(a, amin, amax,routine,msg) \
    if (!(a>=amin && a<=amax)) \
         throw ModelException(routine,msg+ #a" must be in ["#amin","#amax"]");

/** Destructor */
CcmOnlyParameters::~CcmOnlyParameters() 
{}

/** Build a CcmOnlyParameters object based on an old-style
    CCMParameters object - which must not be null */
CcmOnlyParameters::CcmOnlyParameters(CCMParametersConstSP ccmParam) :
    RationalisedCreditEngineParameters(TYPE),
    name(ccmParam->getName()),
    seniorCurve(ccmParam->seniorCurve),
    qm(ccmParam->getQM()),
    indepFactor(ccmParam->getIndependenceFactor()),
    recDispersion(ccmParam->getRecDispersion()),
    recCataFactor(ccmParam->getCatastrophicRecoveryFactor()),
    betaR(ccmParam->getBetaR()),
    betaTweak(ccmParam->getBetaTweak())
{}

void CcmOnlyParameters::validatePop2Object() {
    const char *routine = "CcmOnlyParameters::validatePop2Object";
    string msg("name: ");
    msg += getName() + "," ;
    CHECK_BOUNDS(qm,   -1., 1.,routine,msg);
    CHECK_BOUNDS(betaR,-1., 1.,routine,msg);
    CHECK_BOUNDS(indepFactor,    0., 1.,routine,msg);
    CHECK_BOUNDS(recDispersion,  0., 1.,routine,msg);
    CHECK_BOUNDS(recCataFactor,  0., 1.,routine,msg);
    CHECK_BOUNDS(betaTweak,  -1., 1.,routine,msg);
}

/** Returns the name of this object. This is the name with which
    it is stored in the market data cache and is the name with
    which results (eg tweaks) should be reported against */
string CcmOnlyParameters::getName() const{
    return name;
}

/** Pull out the component assets & correlations from the market data */
void CcmOnlyParameters::getMarket(const IModel* model, const MarketData* market){
    if (!seniorCurve.isEmpty()){
        ICDSParSpreads::getMarketData(model, market, seniorCurve);
    }
}

/** Returns whether the name is modelling stochastic recovery
 * - true : stochastic recovery
 * - false : fixed recovery */
bool CcmOnlyParameters::hasStochasticRecovery() const {
    return (betaR==0.0 && recDispersion==0.0) ? false : true;
}

/** Computes dependence survival proba between startDate and endDate */ 
double CcmOnlyParameters::dependenceSurvivalProba(
    const DateTime& startDate, const DateTime& endDate) const
{
    if (seniorCurve.get() != 0) {
        if (startDate < endDate) {
            // Improve performance here ?
            DefaultRatesSP defaultRates(seniorCurve->defaultRates());
            return defaultRates->calcDefaultPV(
                startDate,
                endDate);
        }
        else {
            return 1.0;
        }
    }
    else {
        return 1.0;
    }
}

/** Returns independence factor */
double CcmOnlyParameters::getIndependenceFactor() const {
    return indepFactor;
}

/**
 * Returns catastrophic recovery factor f such that
 * catastrophic loss = notional * (1 - f * marketRecovery)
 * */
double CcmOnlyParameters::getCatastrophicRecoveryFactor() const {
    return recCataFactor;
}

// [TO BE RETIRED]
double CcmOnlyParameters::getBetaTweak() const {
    return betaTweak;
}

// [TO BE RETIRED]
double CcmOnlyParameters::getQM() const {
    return qm;
}

// [TO BE RETIRED]
double CcmOnlyParameters::getRecDispersion() const {
    return recDispersion;
}

// [TO BE RETIRED]
double CcmOnlyParameters::getBetaR() const {
    return betaR;
}

/**
 * Builds a default CcmOnlyParameters object with all fields defaulted to 0
 * (corresponds to Credit Metrics model)
 * */
CcmOnlyParameters* CcmOnlyParameters::buildDefault() {
    return new CcmOnlyParameters();
}

CcmOnlyParameters::CcmOnlyParameters():
    RationalisedCreditEngineParameters(TYPE),
    seniorCurve(0),
    qm(0.0),
    indepFactor(0.0),
    recDispersion(0.0),
    recCataFactor(0.0),
    betaR(0.0),
    betaTweak(0.0){}

IObject* CcmOnlyParameters::defaultConstructor(){
    return new CcmOnlyParameters();
}

void CcmOnlyParameters::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CcmOnlyParameters, clazz);
    SUPERCLASS(RationalisedCreditEngineParameters);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name, "Identifier");
    FIELD_MAKE_OPTIONAL(name); /* will need to be mandatory when it goes in the
                                  market data cache */
    FIELD(seniorCurve   ,"Senior Curve");
    FIELD_MAKE_OPTIONAL(seniorCurve); // to do: review this
    FIELD(qm            ,"skewed gaussian copula skew");
    FIELD(indepFactor   ,"non cata spread assigned to independence");
    FIELD(recDispersion ,"recovery dispersion (impact variance");
    FIELD(recCataFactor ,"catastrophic recovery factor");
    FIELD(betaR         ,"recovery correlation");
    FIELD(betaTweak     ,"beta relative tweak from historical beta");
    FIELD_MAKE_OPTIONAL(betaTweak);

    // Register "indepFactor" as a field that can be calibrated
    Calibrator::IAdjustable::registerField(
        clazz,
        "indepFactor",
        new Range(ClosedBoundary(0.0), ClosedBoundary(1.0)));

    // Register "betaTweak" as a field that can be calibrated
    Calibrator::IAdjustable::registerField(
        clazz,
        "betaTweak",
        new Range(OpenBoundary(-1.0), OpenBoundary(1.0)));
}

CClassConstSP const CcmOnlyParameters::TYPE = CClass::registerClassLoadMethod(
    "CcmOnlyParameters", typeid(CcmOnlyParameters), load);
    
/** Definition of TYPE for MarketWrapper template class */
DEFINE_TEMPLATE_TYPE(CcmOnlyParametersWrapper);

DRLIB_END_NAMESPACE
