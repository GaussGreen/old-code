#include "edginc/config.hpp"
#include "edginc/CreditMetricsDefaultsModel.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CCMSkew.hpp"
#include "edginc/GaussianMarketFactorModel.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/CondLossDistributionsGenRisklessKey.hpp"
#include "edginc/FunctionOperations.hpp"
#include "edginc/FunctionOperations.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(CreditMetricsDefaultsModelArray);

/** Credit Metrics implementation of ICondLossDistributionsGenKey */
class CMKey:
    public CObject,
    public virtual ICondLossDistributionsGenKey
{
public:
    
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~CMKey(){}
    
    /** Constructor */
    CMKey(double threshold,
        double beta,
        DiscreteDistributionSP condLossDistribution):
            CObject(TYPE),
            threshold(threshold),
            beta(beta),
            condLossDistribution(condLossDistribution) {}
    
    /**
     * Computes a survival probability conditional on
     * a "market factor" value
     * */
    virtual IDistribution1DConstSP conditionalLossDistribution(
        IMarketFactorValueConstSP marketFactorValue) const
    {
        double condSurvProb = conditionalSurvProb(marketFactorValue);

        // Update DiscreteDistribution
        condLossDistribution->setProbability(0, condSurvProb);
        condLossDistribution->setProbability(1, 1.0 - condSurvProb);
            
        return condLossDistribution;            
    }

    /**
     * Computes a survival probability conditional on
     * a "market factor" value
     * */
    virtual double conditionalSurvProb(
        IMarketFactorValueConstSP marketFactorValue) const
    {
        // Expects marketFactor to be a double array
        // NB: use static_cast for performance
        const DoubleArrayMarketFactor* mfValue =
            static_cast<const DoubleArrayMarketFactor*>(marketFactorValue.get());

        // NB: CCMSkew::condSurvivalProba is a very bad name because we
        //     expect this method to return the conditional DEFAULT probability
        //     when we use the convention "negative market factor values = bad state
        //     of the economy".
        double condDefaultProba = CCMSkew::condSurvivalProba(
            1.0,                   // indep survival proba
            threshold,             // gauss survival threshold
            beta,                  // beta
            0.0,                   // skew of M
            mfValue->getValue()[0]); // M

        return 1.0 - condDefaultProba; // ie, cond survival probability
    }

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    double threshold;
    double beta;
    DiscreteDistributionSP condLossDistribution;
};

void CMKey::load(CClassSP& clazz)
{
    clazz->setPrivate();
    REGISTER(CMKey, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICondLossDistributionsGenKey);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(threshold);
    FIELD_NO_DESC(beta);
    FIELD_NO_DESC(condLossDistribution);
}

/** TYPE (for reflection) */        
CClassConstSP const CMKey::TYPE =
CClass::registerClassLoadMethod(
    "CMKey",
    typeid(CMKey),
    CMKey::load);

IObject* CMKey::defaultConstructor()
{
    return new CMKey(
        0.0,
        0.0,
        DiscreteDistributionSP(0));
}

/** TYPE (for reflection) */        
CClassConstSP const CreditMetricsDefaultsModel::TYPE =
CClass::registerClassLoadMethod(
    "CreditMetricsDefaultsModel",
    typeid(CreditMetricsDefaultsModel),
    CreditMetricsDefaultsModel::load);

/** Virtual destructor */
CreditMetricsDefaultsModel::~CreditMetricsDefaultsModel()
{
}

/** Constructor */
CreditMetricsDefaultsModel::CreditMetricsDefaultsModel():
    CObject(TYPE),
    integrator(0),
    marketFactorModel(0){}

/** Public constructor */
CreditMetricsDefaultsModel::CreditMetricsDefaultsModel(
    IIntegratorSP integrator,
    IMarketFactorModelSP marketFactorModel):
        CObject(TYPE),
        integrator(integrator),
        marketFactorModel(marketFactorModel)
{
    validatePop2Object();
}

/** Called immediately after object constructed */
void CreditMetricsDefaultsModel::validatePop2Object()
{
    if (!GaussianMarketFactorModel::TYPE->isInstance(marketFactorModel.get()))
    {
        throw ModelException(
            "CreditMetricsDefaultsModel::validatePop2Object",
            "Type of market factor not supported: " +
            marketFactorModel->getClass()->getName());
    }
}

IObject* CreditMetricsDefaultsModel::defaultConstructor()
{
    return new CreditMetricsDefaultsModel();
}

void CreditMetricsDefaultsModel::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CreditMetricsDefaultsModel, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IConditionalDefaultsModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(integrator, "Market factor integrator - used for internal thresholds calibration");
    FIELD(marketFactorModel, "Market factor model - defines distribution of the market factor");
}

/** [Implements IConditionalDefaultsModel] */
ICondLossDistributionsGenKeySP CreditMetricsDefaultsModel::initialise(
    double defaultProba,
    double expectedLoss,
    double notional,
    CreditEngineParametersConstSP modelParameters,
    const DateTime& startDate,
    const DateTime& endDate) const
{
    if (expectedLoss == 0.0 || defaultProba == 0.0)
    {
        return ICondLossDistributionsGenKeySP(
            new CondLossDistributionsGenRisklessKey());
    }
    
    // Expects modelParameters to be a CmOnlyParameters
    double beta =
        DYNAMIC_CONST_CAST(CmOnlyParameters, modelParameters.get())->getBeta();
    
    // Known closed form for GaussianMarketFactorModel
    // NB: Type of marketFactorModel (expected: GaussianMarketFactorModel)
    //     has been checked in "validatePop2Object()"
    double threshold = N1InverseBetter(defaultProba);

    // Initialise condLossDistribution
    double lossLevel = expectedLoss / defaultProba;
    DiscreteDistributionSP condLossDistribution(
        new DiscreteDistribution(lossLevel, 0.0));

    return ICondLossDistributionsGenKeySP(
        new CMKey(threshold, beta, condLossDistribution));
}

/** [Implements IConditionalDefaultsModel] */
double CreditMetricsDefaultsModel::integrateCondFunction(
    const MFunctionND* condFunction,
    ICondLossDistributionsGenKeyArrayConstSP condKeys,
    const DateTime& time) const
{
    FunctionNDDoubleConstSP density = marketFactorModel->density(time);
    
    MFunctionNDSP integrand = FunctionOperations::multiply(
        *condFunction, *density, density->getIntervals());
    
    return DYNAMIC_CAST(CDouble,
        integrator->integrate(integrand.get()).get())->doubleValue();
}

/** [Implements IConditionalDefaultsModel] */
bool CreditMetricsDefaultsModel::isCondFunctionIntegrationTimeDependent() const {
    return marketFactorModel->isDensityTimeDependent();
}

/** [Implements IConditionalDefaultsModel] */
int CreditMetricsDefaultsModel::marketFactorDimension() const
{
    return marketFactorModel->dimension();
}

/** [Implements IConditionalDefaultsModel] */
CClassConstSP CreditMetricsDefaultsModel::engineParamsType() const
{
    return CmOnlyParameters::TYPE;
}

/* external symbol to allow class to be forced to be linked in */
bool CreditMetricsDefaultsModelLoad(){
    return (CreditMetricsDefaultsModel::TYPE != 0);
}

DRLIB_END_NAMESPACE
