#include "edginc/config.hpp"
#include "edginc/CompositeCopulaDefaultsModel.hpp"
#include "edginc/CreditMetricsDefaultsModel.hpp"
#include "edginc/RFLDefaultsModel.hpp"
#include "edginc/CmCcmRflParameters.hpp"
#include "edginc/CmCcmParameters.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Gaussian2DMarketFactorModel.hpp"
#include "edginc/CondLossDistributionsGenRisklessKey.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"
#include <set>
#include "edginc/FunctionOperations.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(CompositeCopulaDefaultsModelArray);

START_PUBLIC_ENUM_DEFINITION(CompositeCopulaDefaultsModel::CopulaTypes,
    "Types of copula models");
ENUM_VALUE(CompositeCopulaDefaultsModel::CREDIT_METRICS, "Credit Metrics model");
ENUM_VALUE(CompositeCopulaDefaultsModel::RFL, "RFL model");
END_ENUM_DEFINITION(CompositeCopulaDefaultsModel::CopulaTypes);

/** Composite Copula implementation of ICondLossDistributionsGenKey */
class CompositeCopulaKey: 
    public CObject,
    public virtual ICondLossDistributionsGenKey
{
public:
    
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~CompositeCopulaKey(){}
    
    /** Constructor */
    CompositeCopulaKey(
        ICondLossDistributionsGenKeySP mixtureCopulaKey,
        double cataThreshold,
        double catastrophicLoss,
        double nonCatastrophicLoss,
        double independenceSP):
            CObject(TYPE),
            mixtureCopulaKey(mixtureCopulaKey),
            cataThreshold(cataThreshold),
            catastrophicLoss(catastrophicLoss),
            nonCatastrophicLoss(nonCatastrophicLoss),
            independenceSP(independenceSP)
    {
        // Builds catastrophicLossDistribution once for all
        DoubleArraySP values(new DoubleArray(1));
        DoubleArraySP probas(new DoubleArray(1));
        (*values)[0] = catastrophicLoss;
        (*probas)[0] = 1.0;
        catastrophicLossDistribution.reset(
            new DiscreteDistribution(values, probas));
    }
    
    /**
     * Computes a survival probability conditional on
     * a "market factor" value
     * */
    virtual IDistribution1DConstSP conditionalLossDistribution(
        IMarketFactorValueConstSP marketFactorValue) const
    {
        // Expects marketFactor to be a double array
        // NB: use static_cast for performance
        const DoubleArrayMarketFactor* mfValue =
            static_cast<const DoubleArrayMarketFactor*>(marketFactorValue.get());

        // Tests if we are in the "catastrophic scenario"
        if (mfValue->getValue()[1] < cataThreshold)
        {
            // Catastrophic scenario, just return the catastrophicLossDistribution !
            return catastrophicLossDistribution;
        }
        else
        {
            // Non catastrophic scenario
            
            // Retrieves "mixture copula" loss distribution
            // NB: may be non-binary if using e.g. a local recovery model
            DiscreteDistributionConstSP mixtureLD = 
                mixtureCopulaKey->conditionalLossDistribution(marketFactorValue)->discretise();

            // Merges "mixture copula" loss distribution with "independence copula"
            // using the following rule:
            // - use nonCatastrophicLoss if
            //   no catastrophic default (always true here because sp1(z)=1) AND
            //   no "mixture default" [this case has a probability sp1(z)*sp2(z)*(1-sp3)]
            // - use "mixture" loss (distribution) if
            //   no catastrophic default (always true here because sp1(z)=1) AND
            //   at LEAST one "mixture default" [this case has a probability sp1(z)*(1-sp2(z))]

            if (independenceSP == 1.0)
            {
                // degenerated case
                return mixtureLD;
            }

            DoubleArrayConstSP mixtureValues = mixtureLD->getValues();
            DoubleArrayConstSP mixtureProbas = mixtureLD->getProbabilities();

            DoubleArraySP mergedValues, mergedProbas;

            // Checks if mixtureLD has a value = 0.0 ("zero loss")
            int indexZeroLoss;
            bool hasZeroLoss = mixtureLD->search(0.0, indexZeroLoss, 3e-15);

            // Checks if mixtureLD has a value = nonCatastrophicLoss
            int indexNCLoss;
            bool hasNCLoss = mixtureLD->search(nonCatastrophicLoss, indexNCLoss, 3e-15);

            if (hasZeroLoss)
            {
                double mixtureSP = (*mixtureProbas)[indexZeroLoss];
                
                if (mixtureSP == 0.0)
                {
                    // degenerated case where we can ignore independenceSP
                    return mixtureLD;
                }

                if (hasNCLoss)
                {
                    // mixtureValues has already all relevant elements
                    mergedValues.reset(
                        new DoubleArray(mixtureValues->begin(), mixtureValues->end()));
                    mergedProbas.reset(
                        new DoubleArray(mixtureProbas->begin(), mixtureProbas->end()));
                    
                    // update nonCatastrophicLoss (NB: use addition)
                    (*mergedProbas)[indexNCLoss] += mixtureSP * (1.0 - independenceSP);
                }
                else
                {
                    // need to add nonCatastrophicLoss to mergedValues
                    int nbMergedValues = mixtureValues->size() + 1;
                    mergedValues.reset(new DoubleArray(nbMergedValues));
                    mergedProbas.reset(new DoubleArray(nbMergedValues));

                    // merge
                    int i;
                    for (i = 0; i < indexNCLoss && i < mixtureValues->size(); ++i) {
                        (*mergedValues)[i] = (*mixtureValues)[i];
                        (*mergedProbas)[i] = (*mixtureProbas)[i];
					}
                    (*mergedValues)[indexNCLoss] = nonCatastrophicLoss;
                    (*mergedProbas)[indexNCLoss] = mixtureSP * (1.0 - independenceSP);
                    for (i = indexNCLoss + 1; i < nbMergedValues; ++i) {
                        (*mergedValues)[i] = (*mixtureValues)[i-1];
                        (*mergedProbas)[i] = (*mixtureProbas)[i-1];
                    }
                }
                // update "zero loss" (NB: use multiplication)
                (*mergedProbas)[indexZeroLoss] *= independenceSP;
            }
            else
            {
                // degenerated case (mixtureSP == 0.0) where we can ignore independenceSP
                return mixtureLD;
            }

            return DiscreteDistributionSP(new DiscreteDistribution(
                mergedValues, mergedProbas));
        }
    }

    /**
     * Computes a survival probability conditional on
     * a "market factor" value
     * */
    virtual double conditionalSurvProb(
        IMarketFactorValueConstSP marketFactorValue) const
    {
        throw ModelException ("Internal error: method "
                              "CompositeCopulaKey::conditionalSurvProb "
                              "not yet supported");
    }

    /** Access to catastrophic threshold */
    double getCatastrophicSurvivalProba() const
    {
        return cataThreshold;
    }
    
    /** Access to "mixture copula" key */ 
    ICondLossDistributionsGenKeyConstSP getMixtureCopulaKey() const
    {
        return mixtureCopulaKey;
    }

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    // input fields
    ICondLossDistributionsGenKeyConstSP mixtureCopulaKey;
    double cataThreshold;
    double catastrophicLoss;
    double nonCatastrophicLoss;
    double independenceSP; // independence survival probability "s3"

    // internally computed fields
    DiscreteDistributionSP catastrophicLossDistribution;
};

void CompositeCopulaKey::load(CClassSP& clazz)
{
    clazz->setPrivate();
    REGISTER(CompositeCopulaKey, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICondLossDistributionsGenKey);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(mixtureCopulaKey);
    FIELD_NO_DESC(cataThreshold);
    FIELD_NO_DESC(catastrophicLoss);
    FIELD_NO_DESC(nonCatastrophicLoss);
    FIELD_NO_DESC(independenceSP);
    FIELD_NO_DESC(catastrophicLossDistribution);
    FIELD_MAKE_TRANSIENT(catastrophicLossDistribution);
}

/** TYPE (for reflection) */        
CClassConstSP const CompositeCopulaKey::TYPE =
CClass::registerClassLoadMethod(
    "CompositeCopulaKey",
    typeid(CompositeCopulaKey),
    CompositeCopulaKey::load);

IObject* CompositeCopulaKey::defaultConstructor()
{
    return new CompositeCopulaKey(
        ICondLossDistributionsGenKeySP(0),
        0.0,
        0.0,
        0.0,
        0.0);
}

/** TYPE (for reflection) */        
CClassConstSP const CompositeCopulaDefaultsModel::TYPE =
CClass::registerClassLoadMethod(
    "CompositeCopulaDefaultsModel",
    typeid(CompositeCopulaDefaultsModel),
    CompositeCopulaDefaultsModel::load);

/** Virtual destructor */
CompositeCopulaDefaultsModel::~CompositeCopulaDefaultsModel()
{
}

/** Constructor */
CompositeCopulaDefaultsModel::CompositeCopulaDefaultsModel():
    CObject(TYPE),
    integrator(0),
    marketFactorModel(0),
    mixtureCopulaType(CREDIT_METRICS){}

/** Called immediately after object constructed */
void CompositeCopulaDefaultsModel::validatePop2Object()
{
    static const string method("CompositeCopulaDefaultsModel::validatePop2Object");        
    try
    {
        // builds mixtureCopula model using consistent integrator and marketFactorModel
        switch (mixtureCopulaType) {
            case CREDIT_METRICS:
                mixtureCopula.reset(
                    new CreditMetricsDefaultsModel(integrator, marketFactorModel));
                break;
            case RFL:
                mixtureCopula.reset(
                    new RFLDefaultsModel(integrator, marketFactorModel));
                break;
            default:
                throw ModelException(
                    method,
                    "Type of mixture copula not supported: " +
                    BoxedEnum<CopulaTypes>::toString(mixtureCopulaType) +
                    ".");
                break;
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

IObject* CompositeCopulaDefaultsModel::defaultConstructor()
{
    return new CompositeCopulaDefaultsModel();
}

void CompositeCopulaDefaultsModel::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CompositeCopulaDefaultsModel, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IConditionalDefaultsModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(integrator, "Market factor integrator - used for internal thresholds calibration");
    FIELD(marketFactorModel, "Market factor model - defines distribution of the market factor");
    FIELD(mixtureCopula, "Internally built 'mixture copula' model");
    FIELD_MAKE_TRANSIENT(mixtureCopula);
    FIELD(mixtureCopulaType,
        "Define the type of the 'mixture copula' (eg: Credit Metrics or RFL)");
}

/** [Implements IConditionalDefaultsModel] */
ICondLossDistributionsGenKeySP CompositeCopulaDefaultsModel::initialise(
    double defaultProba,
    double expectedLoss,
    double notional,
    CreditEngineParametersConstSP modelParameters,
    const DateTime& startDate,
    const DateTime& endDate) const
{
    static const string method("CompositeCopulaDefaultsModel::initialise");        
    
    if (expectedLoss == 0.0 || defaultProba == 0.0)
    {
        return ICondLossDistributionsGenKeySP(
            new CondLossDistributionsGenRisklessKey());
    }
    
    // -----------------------------------
    // Retrieves relevant model parameters
    // -----------------------------------
    
    // Retrieves CcmOnlyParameters within modelParameters
    const CcmOnlyParameters* ccmParams = DYNAMIC_CONST_CAST(
        CcmOnlyParameters, modelParameters->getEngineParams(CcmOnlyParameters::TYPE).get());

    // Retrieves "mixture copula" parameters within modelParameters
    CreditEngineParametersConstSP mixtureCopulaParameters;
    switch (mixtureCopulaType) {
        case CREDIT_METRICS:
            // Retrieves CmOnlyParameters within modelParameters
            mixtureCopulaParameters = modelParameters->getEngineParams(CmOnlyParameters::TYPE);
            break;
        case RFL:
        {
            // Retrieves CmOnlyParameters within modelParameters
            CmOnlyParametersConstSP ccmParam = CmOnlyParametersConstSP::attachToRef(
                DYNAMIC_CONST_CAST(CmOnlyParameters, modelParameters->getEngineParams(CmOnlyParameters::TYPE).get()));

            // Retrieves RflOnlyParameters within modelParameters
            RflOnlyParametersConstSP rflParam = RflOnlyParametersConstSP::attachToRef(
                DYNAMIC_CONST_CAST(RflOnlyParameters, modelParameters->getEngineParams(RflOnlyParameters::TYPE).get()));
            
            mixtureCopulaParameters.reset(new CmRflParameters(ccmParam, rflParam));
            break;
        }
        default:
            throw ModelException(method,
                "Type of mixture copula not supported: " +
                BoxedEnum<CopulaTypes>::toString(mixtureCopulaType) +
                ".");
    }
    
    // --------------------------------------------------------
    // Splits survival proba sp into
    // sp1 (dependence), sp2 ("mixture") and sp3 (independence)
    // --------------------------------------------------------

    double sp = 1.0 - defaultProba;

    double sp1 = ccmParams->dependenceSurvivalProba(startDate, endDate);
    
    if (sp1 < sp)
    {
        throw ModelException(method,
            "Dependence survival probability (" +
            Format::toString(sp1) +
            ") is smaller than survival probability (" +
            Format::toString(sp) +
            ") for period " +
            startDate.toString() +
            " to " +
            endDate.toString() +
            ".");
    }

    double independenceFactor = ccmParams->getIndependenceFactor();

    double sp3 = pow(sp / sp1, independenceFactor);
    double sp2 = sp / (sp1 * sp3); // to ensure sp = sp1*sp2*sp3

    // ------------------------------------
    // Calibrates non-catastrophic recovery
    // ------------------------------------
    //
    // Here we have a closed form:
    // - if at least 1 "catastrophic" default, use cataRecovery
    //   (adds notional * (1-cataRecovery) * (1-sp1(Z)) to the conditional expected loss)
    // - if at least 1 default but no "catastrophic" default, use nonCataRecovery
    //   (adds notional            *
    //         (1-nonCataRecovery) *
    //         sp1(Z)              * [no "catastrophic" default]
    //         (1-sp2(Z)*sp3(Z))   * [but at least 1 default]
    //    to the conditional expected loss)
    //
    // Then we integrate the conditional expected loss over the market
    // factor Z and get:
    // EL = notional * [(1-cataRecovery) * (1 - sp1) + (1-nonCataRecovery) * (sp1 - sp)]
    //
    // NB: would require numerical integration over market factor for
    //     more generic composite copula models
    
    if (notional == 0.0 || defaultProba == 0.0)
    {
        throw ModelException(method,
            "Asset is riskless (notional = " +
            Format::toString(notional) +
            ", default probability = " +
            Format::toString(defaultProba) +
            ")");
    }
    double cataRecovery = ccmParams->getCatastrophicRecoveryFactor() *
        (1.0 - expectedLoss / (notional * defaultProba));
    double nonCataRecovery = 1.0 - 
        ((expectedLoss / notional) - (1.0 - cataRecovery) * (1.0 - sp1)) /
        (sp1 - sp);
    
    // ---------------------
    // Calibrates thresholds
    // ---------------------

    // "catastrophic" threshold, known closed form !
    double cataThreshold = 1.0 - sp1;
    
    // "mixture copula" threshold
    ICondLossDistributionsGenKeySP mixtureCopulaKey = mixtureCopula->initialise(
        1.0 - sp2, // default proba
        notional * (1.0 - sp2) * (1.0 - nonCataRecovery), // expected loss
        notional, // notional
        mixtureCopulaParameters, // parameters
        startDate,
        endDate);

    return ICondLossDistributionsGenKeySP(
        new CompositeCopulaKey(
            mixtureCopulaKey,
            cataThreshold,
            notional * (1.0 - cataRecovery),
            notional * (1.0 - nonCataRecovery),
            sp3));
}

/** [Implements IConditionalDefaultsModel] */
double CompositeCopulaDefaultsModel::integrateCondFunction(
    const MFunctionND* condFunction,
    ICondLossDistributionsGenKeyArrayConstSP condKeys,
    const DateTime& time) const
{
    static const string method("CompositeCopulaDefaultsModel::integrateCondFunction"); 
    
    // Retrieves all catastrophic survival probabilities and
    // "mixture keys" inside CompositeCopulaKey
    set<double> allCataSP;
    allCataSP.insert(0.0);
    allCataSP.insert(1.0);

    if (!condKeys) {
        throw ModelException(method, "Cannot use this model (no keys provided)");
    }
    ICondLossDistributionsGenKeyArraySP mixtureCopulaKeys(
        new ICondLossDistributionsGenKeyArray(condKeys->size()));
    for (int i = 0; i < (int) condKeys->size(); ++i) {
        const CompositeCopulaKey* key =
            dynamic_cast<const CompositeCopulaKey*>((*condKeys)[i].get());
        if (key == 0)
        {
            throw ModelException(method,
                "Internal error: unable to cast to CompositeCopulaKey");
        }                

        allCataSP.insert(key->getCatastrophicSurvivalProba());

        // const_cast is not very nice, ideally mixtureCopulaKeys should be
        // an array<ICondLossDistributionsGenKeyConstSP, ICondLossDistributionsGenKey>.
        // Note that mixtureCopulaKeys is then passed below as a
        // ICondLossDistributionsGenKeyArrayConstSP in mixtureCopula->integrateCondEL
        (*mixtureCopulaKeys)[i] =
            ICondLossDistributionsGenKeySP::constCast(key->getMixtureCopulaKey());
	}
    
    // Integrates
    double integral = 0.0;
    double b, weight, y, partialIntegral;
    set<double>::const_iterator iter = allCataSP.begin();
    double a = (*iter);
    ++iter;
    for (; iter != allCataSP.end(); ++iter) {
        // Builds a "partial" 1D function for each catastrophic survival probabilities interval
        b = *(iter);
        weight = b-a;
        y = (a+b)/2;
        MFunctionNDSP partialFunction = FunctionOperations::partialFunction(*condFunction, 1, y);
        a = b;

        // Computes "partial" 1D integral
        partialIntegral = 
            mixtureCopula->integrateCondFunction(partialFunction.get(), 
                                                 mixtureCopulaKeys, 
                                                 time);
        
        // Updates (global) 2D integral
        integral += weight * partialIntegral;
	}
    return integral;
}

/** [Implements IConditionalDefaultsModel] */
bool CompositeCopulaDefaultsModel::isCondFunctionIntegrationTimeDependent() const {
    return true;    // integrateCondFunction uses the (time-specific) condKeys
}


/** [Implements IConditionalDefaultsModel] */
int CompositeCopulaDefaultsModel::marketFactorDimension() const
{
    // NB: "1" corresponds to the "catastrophic" dimension where we
    //     know a simple closed form to perform the integration
    return mixtureCopula->marketFactorDimension() + 1;
}

/**
 * Returns the type of CreditEngineParameters expected
 * by this IConditionalDefaultsModel
 * */
CClassConstSP CompositeCopulaDefaultsModel::engineParamsType() const
{
    switch (mixtureCopulaType) {
        case CREDIT_METRICS:
            return CmCcmParameters::TYPE;
        case RFL:
            return CmCcmRflParameters::TYPE;
        default:
            throw ModelException(
                "CompositeCopulaDefaultsModel::engineParamsType",
                "Type of mixture copula not supported: " +
                BoxedEnum<CopulaTypes>::toString(mixtureCopulaType) +
                ".");
    }
}

/* external symbol to allow class to be forced to be linked in */
bool CompositeCopulaDefaultsModelLoad(){
    return (CompositeCopulaDefaultsModel::TYPE != 0);
}

DRLIB_END_NAMESPACE
