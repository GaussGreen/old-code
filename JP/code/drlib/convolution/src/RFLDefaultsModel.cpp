#include "edginc/config.hpp"
#include "edginc/RFLDefaultsModel.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CCMSkew.hpp"
#include "edginc/GaussianMarketFactorModel.hpp"
#include "edginc/CmRflParameters.hpp"
#include "edginc/CondLossDistributionsGenRisklessKey.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/FunctionOperations.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(RFLDefaultsModelArray);

/** RFL implementation of ICondLossDistributionsGenKey */
class RFLKey:
    public CObject,
    public virtual ICondLossDistributionsGenKey
{
public:
    
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~RFLKey(){}
    
    /** Constructor */
    RFLKey(double threshold,
        double betaHist,
        RflOnlyParametersConstSP rflParams,
        DiscreteDistributionSP condLossDistribution):
            CObject(TYPE),
            threshold(threshold),
            betaHist(betaHist),
            rflParams(rflParams),
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
        const DoubleArrayMarketFactor* mfValue =
            static_cast<const DoubleArrayMarketFactor*>(marketFactorValue.get());

        double condDefaultProba;
    
        // NB: RflOnlyParameters::condSurvivalProba is a very bad name because we
        //     expect this method to return the conditional DEFAULT probability
        //     when we use the convention "negative market factor values = bad state
        //     of the economy".
        (*rflParams).condSurvivalProba(
            1.0,                   // indep survival proba
            threshold,             // gauss survival threshold
            betaHist,              // beta
            mfValue->getValue()[0],// M
            condDefaultProba);

        return 1.0 - condDefaultProba; // ie, cond survival probability
    }

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    double threshold;
    double betaHist;
    RflOnlyParametersConstSP rflParams;
    DiscreteDistributionSP condLossDistribution;
};

void RFLKey::load(CClassSP& clazz)
{
    clazz->setPrivate();
    REGISTER(RFLKey, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICondLossDistributionsGenKey);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(threshold);
    FIELD_NO_DESC(betaHist);
    FIELD_NO_DESC(rflParams);
    FIELD_NO_DESC(condLossDistribution);
}

/** TYPE (for reflection) */        
CClassConstSP const RFLKey::TYPE =
CClass::registerClassLoadMethod(
    "RFLKey",
    typeid(RFLKey),
    RFLKey::load);

IObject* RFLKey::defaultConstructor()
{
    return new RFLKey(
        0.0,
        0.0,
        RflOnlyParametersConstSP(0),
        DiscreteDistributionSP(0));
}

/** TYPE (for reflection) */        
CClassConstSP const RFLDefaultsModel::TYPE =
CClass::registerClassLoadMethod(
    "RFLDefaultsModel",
    typeid(RFLDefaultsModel),
    RFLDefaultsModel::load);

/** Virtual destructor */
RFLDefaultsModel::~RFLDefaultsModel()
{
}

/** Constructor */
RFLDefaultsModel::RFLDefaultsModel():
    CObject(TYPE),
    integrator(0),
    marketFactorModel(0) {}

/** Public constructor */
RFLDefaultsModel::RFLDefaultsModel(
    IIntegratorSP integrator,
    IMarketFactorModelSP marketFactorModel):
        CObject(TYPE),
        integrator(integrator),
        marketFactorModel(marketFactorModel)
{
    validatePop2Object();
}

/** Called immediately after object constructed */
void RFLDefaultsModel::validatePop2Object()
{
    if (!GaussianMarketFactorModel::TYPE->isInstance(marketFactorModel.get()))
    {
        throw ModelException(
            "RFLDefaultsModel::validatePop2Object",
            "Type of market factor not supported: " +
            marketFactorModel->getClass()->getName());
    }
}

IObject* RFLDefaultsModel::defaultConstructor()
{
    return new RFLDefaultsModel();
}

void RFLDefaultsModel::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(RFLDefaultsModel, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IConditionalDefaultsModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(integrator, "Market factor integrator - used for internal thresholds calibration");
    FIELD(marketFactorModel, "Market factor model - defines distribution of the market factor");
}

/** [Implements IConditionalDefaultsModel] */
ICondLossDistributionsGenKeySP RFLDefaultsModel::initialise(
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
    
    // Expects modelParameters to be a CmRflParameters
    const CmRflParameters* cmRflParameters = DYNAMIC_CONST_CAST(CmRflParameters, modelParameters.get());

    // get the historical beta from the cm parameters
    double betaHist = (cmRflParameters->getCmParameters())->getBeta();

    // get the rfl parameters
    RflOnlyParametersConstSP rflParams = cmRflParameters->getRflParameters();

    double threshold;

    Key key(defaultProba, betaHist, rflParams);

    ThresholdsCache::const_iterator iter = thresholdsCache.find(key);
    if (iter != thresholdsCache.end())
    {
        // already computed, just return the cached value
        threshold = iter->second;
    }
    else
    {
        // compute...

        // This is the time consuming operation !
        // NB1: Type of marketFactorModel (expected: GaussianMarketFactorModel)
        //      has been checked in "validatePop2Object()"
        // NB2: We pass a DEFAULT probability as first input in threshold(...) method
        //      because we use the convention "negative market factor values = bad state
        //      of the economy".
        (*rflParams).threshold(defaultProba, threshold, betaHist);

        // ... and store in the cache
        
        // Need to clone rflParams to ensure nobody can modify them
        Key newKey(defaultProba, betaHist, IObjectSP(rflParams->clone()));
        thresholdsCache.insert(MapEntry(newKey, threshold));
    }
    
    // Initialise condLossDistribution
    double lossLevel = expectedLoss / defaultProba;
    DiscreteDistributionSP condLossDistribution(
        new DiscreteDistribution(lossLevel, betaHist));
    
    return ICondLossDistributionsGenKeySP(
        new RFLKey(threshold, betaHist, rflParams, condLossDistribution));
}

/** [Implements IConditionalDefaultsModel] */
double RFLDefaultsModel::integrateCondFunction(
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
bool RFLDefaultsModel::isCondFunctionIntegrationTimeDependent() const {
    return marketFactorModel->isDensityTimeDependent();
}

/** [Implements IConditionalDefaultsModel] */
int RFLDefaultsModel::marketFactorDimension() const
{
    return marketFactorModel->dimension();
}

/** [Implements IConditionalDefaultsModel] */
CClassConstSP RFLDefaultsModel::engineParamsType() const
{
    return CmRflParameters::TYPE;
}
/** Override clone() to copy the local cache */
IObject* RFLDefaultsModel::clone() const
{
    IObject* myCopy = CObject::clone();
    RFLDefaultsModel& rflDefaultsModel = 
        dynamic_cast<RFLDefaultsModel&>(*myCopy);
    // copy the cache !
    rflDefaultsModel.thresholdsCache = thresholdsCache;
    return myCopy;
}

/* external symbol to allow class to be forced to be linked in */
bool RFLDefaultsModelLoad(){
    return (RFLDefaultsModel::TYPE != 0);
}

/** Key constructor (for the thresholds cache) */
RFLDefaultsModel::Key::Key(double defaultProba,
    double betaHist,
    IObjectConstSP rflParams):
        defaultProba(defaultProba),
        betaHist(betaHist),
        rflParams(rflParams) {}

/** Key "==" operator (for the thresholds cache) */
bool RFLDefaultsModel::Key::operator==(const Key& key) const
{
    return (defaultProba == key.defaultProba &&
        betaHist == key.betaHist &&
        rflParams->equalTo(key.rflParams.get()));
}

/** KeyHash "()" operator (for the thresholds cache) */
size_t RFLDefaultsModel::KeyHash::operator()(const Key& key) const {
    return (CDouble::hashCode(key.defaultProba) ^
        CDouble::hashCode(key.betaHist) ^
        key.rflParams->hashCode());
}

DRLIB_END_NAMESPACE
