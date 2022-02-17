#include "edginc/config.hpp"
#include "edginc/RFLMixtureDefaultsModel.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CCMSkew.hpp"
#include "edginc/GaussianMarketFactorModel.hpp"
#include "edginc/CmRflParameters.hpp"
#include "edginc/CondLossDistributionsGenRisklessKey.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/FunctionOperations.hpp"
#include "edginc/GridIntegrator.hpp"
#include "edginc/mathlib.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(RFLMixtureDefaultsModelArray);

/** RFL implementation of ICondLossDistributionsGenKey */
class RFLMixtureKey:
    public CObject,
    public virtual ICondLossDistributionsGenKey
{
public:
    static CClassConstSP const TYPE;

    /** Constructor */
    RFLMixtureKey(double barrier,
           RflMixtureParametersConstSP rflParams,
           DiscreteDistributionSP condLossDistribution):
    CObject(TYPE),
    m_barrier(barrier),
    m_rflParams(rflParams),
    m_condLossDistribution(condLossDistribution) {}
    
    /**
     * Computes a survival probability conditional on
     * a "market factor" value
     * */
    IDistribution1DConstSP conditionalLossDistribution(
        IMarketFactorValueConstSP marketFactorValue) const
    {
        // Expects marketFactor to be a double array
        const DoubleArrayMarketFactor* mfValue =
            static_cast<const DoubleArrayMarketFactor*>(marketFactorValue.get());

        double condSP;

        m_rflParams->condSurvivalProba(m_barrier,
                                    mfValue->getValue()[0],
                                    condSP);

        // Update DiscreteDistribution
        m_condLossDistribution->setProbability(0, condSP);
        m_condLossDistribution->setProbability(1, 1.0 - condSP);
            
        return m_condLossDistribution;            
    }

    // Write thresholds to provided container
    void Thresholds(set<double>& thresholds) const { m_rflParams->Thresholds(thresholds); }

private:
    double m_barrier;
    RflMixtureParametersConstSP m_rflParams;
    DiscreteDistributionSP m_condLossDistribution;

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

};

void RFLMixtureKey::load(CClassSP& clazz)
{
    clazz->setPrivate();
    REGISTER(RFLMixtureKey, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICondLossDistributionsGenKey);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(m_barrier);
    FIELD_NO_DESC(m_rflParams);
    FIELD_NO_DESC(m_condLossDistribution);
}

/** TYPE (for reflection) */        
CClassConstSP const RFLMixtureKey::TYPE =
CClass::registerClassLoadMethod(
    "RFLMixtureKey",
    typeid(RFLMixtureKey),
    RFLMixtureKey::load);

IObject* RFLMixtureKey::defaultConstructor()
{
    return new RFLMixtureKey(
        0.0,
        RflMixtureParametersConstSP(0),
        DiscreteDistributionSP(0));
}

/** TYPE (for reflection) */        
CClassConstSP const RFLMixtureDefaultsModel::TYPE =
CClass::registerClassLoadMethod(
    "RFLMixtureDefaultsModel",
    typeid(RFLMixtureDefaultsModel),
    RFLMixtureDefaultsModel::load);

/** Constructor */
RFLMixtureDefaultsModel::
RFLMixtureDefaultsModel()
: CObject(TYPE),
  m_c(0.0), m_w(0.0), m_factorShift(0.0),
  m_lowerBound(-7.5), m_upperBound(7.5), m_numNodes(101)
{
  // empty
}

/** Public constructor */
RFLMixtureDefaultsModel::
RFLMixtureDefaultsModel(double c,
                 double w,
                 double factorShift,
                 double lowerBound,
                 double upperBound,
                 unsigned int numNodes)
: CObject(TYPE), m_c(c), m_w(w), m_factorShift(factorShift), m_lowerBound(lowerBound), m_upperBound(upperBound), m_numNodes(numNodes)
{
  // empty
}

/** Called immediately after object constructed */
void RFLMixtureDefaultsModel::validatePop2Object()
{
/*    if (!GaussianMarketFactorModel::TYPE->isInstance(marketFactorModel.get()))
    {
        throw ModelException(
            "RFLMixtureDefaultsModel::validatePop2Object",
            "Type of market factor not supported: " +
            marketFactorModel->getClass()->getName());
    }*/
}

IObject* RFLMixtureDefaultsModel::defaultConstructor()
{
    return new RFLMixtureDefaultsModel();
}

void RFLMixtureDefaultsModel::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(RFLMixtureDefaultsModel, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IConditionalDefaultsModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(m_c, "Mixing shift");
    FIELD(m_w, "Mixing weight");
    FIELD(m_factorShift, "Factor shift");
    FIELD_MAKE_OPTIONAL(m_factorShift);
    FIELD(m_lowerBound, "Lower bound for factor quadrature");
    FIELD_MAKE_OPTIONAL(m_lowerBound);
    FIELD(m_upperBound, "Upper bound for factor quadrature");
    FIELD_MAKE_OPTIONAL(m_upperBound);
    FIELD(m_numNodes, "Number of nodes in factor quadrature");
    FIELD_MAKE_OPTIONAL(m_numNodes);
}

/** [Implements IConditionalDefaultsModel] */
ICondLossDistributionsGenKeySP RFLMixtureDefaultsModel::initialise(
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
    
    // Expects modelParameters to be a RFLMixtureParameters
    const RflMixtureParameters* const params = DYNAMIC_CONST_CAST(RflMixtureParameters, modelParameters.get());

    RflMixtureParametersConstSP rflParams = RflMixtureParametersConstSP::attachToRef(params);

    Key key(defaultProba, rflParams);

    double barrier;

    ThresholdsCache::const_iterator iter = barriersCache.find(key);
    if (iter != barriersCache.end())
    {
        // already computed, just return the cached value
        barrier = iter->second;
    }
    else
    {
        // compute...

        // This is the time consuming operation !
        // NB: Type of marketFactorModel (expected: GaussianMarketFactorModel)
        //     has been checked in "validatePop2Object()"
        rflParams->defaultBarrier(1.0 - defaultProba, barrier, m_c, m_w);

        // ... and store in the cache
        
        // Need to clone rflParams to ensure nobody can modify them
        Key newKey(defaultProba, IObjectSP(rflParams->clone()));
        barriersCache.insert(MapEntry(newKey, barrier));
    }
    
    // Initialise condLossDistribution
    double LGD = expectedLoss / defaultProba;
    DiscreteDistributionSP condLossDistribution(
        new DiscreteDistribution(LGD, defaultProba));
    
    return ICondLossDistributionsGenKeySP(
        new RFLMixtureKey(barrier, rflParams, condLossDistribution));
}

// Implement function class to support shifts of the factor distribution
// Only used in function below
class ShiftedGaussian {
public:
  ShiftedGaussian(double shift) : m_shift(shift) { }
  double operator()(double x) const { return N1Density(x-m_shift); }
private:
  double m_shift;
};

/** [Implements IConditionalDefaultsModel] */
double RFLMixtureDefaultsModel::integrateCondFunction(
    const MFunctionND* convolutedCondEL,
    ICondLossDistributionsGenKeyArrayConstSP condELKeys,
    const DateTime& time) const
{
  // Check if we have a real mixture
  const bool mixture = (fabs(m_c) <= 1.0e-6 || fabs(m_w) <= 1.0e-10) ? false : true;
   
  // Find all discontinuities in factor-dependence
  set<double> discSet0, discSet1;
  set<double> thresholds;
  ICondLossDistributionsGenKeyArraySP RFLKeys(
        new ICondLossDistributionsGenKeyArray(condELKeys->size()));
  for (int i = 0; i < (int) condELKeys->size(); ++i) 
  {
    const RFLMixtureKey* key = dynamic_cast<const RFLMixtureKey*>((*condELKeys)[i].get());
    if (key == 0)
      throw ModelException("RFLMixtureDefaultsModel::integrateCondEL", "Internal error: unable to cast to RFLMixtureKey");
    // Get thresholds for this name
    key->Thresholds(thresholds);
    discSet0.insert(thresholds.begin(),thresholds.end());
  }
  if (mixture)
    for (set<double>::const_iterator iD = discSet0.begin(); iD != discSet0.end(); ++iD)
      discSet1.insert(*iD - m_c);

  // Construct first integrator
  GridIntegrator integrator0 = GridIntegrator::Simpson(m_numNodes, m_lowerBound + m_factorShift, m_upperBound + m_factorShift, discSet0, 1.0e-10);
  integrator0.ScaleWeightsByWeightingFunction(ShiftedGaussian(m_factorShift));

  // Declare outputs
  vector<double> el0(1,0.0), el1(1,0.0);

  // Do first integral
  integrator0(*convolutedCondEL, el0);

  // Do second integral if required
  if (mixture)
  {
    GridIntegrator integrator1 = GridIntegrator::Simpson(m_numNodes, m_lowerBound + m_factorShift, m_upperBound + m_factorShift, discSet1, 1.0e-10);
    integrator1.ScaleWeightsByWeightingFunction(ShiftedGaussian(m_factorShift));

    CDoubleArray shift(1, m_c);
    integrator1(*FunctionOperations::shift(shift, *convolutedCondEL), el1);
  }

  // Return weighted sum of integrals
  const double res = mixture? (1.0 - m_w) * el0[0] + m_w * el1[0] : el0[0];
  return res;
}

/** [Implements IConditionalDefaultsModel] */
int RFLMixtureDefaultsModel::marketFactorDimension() const
{
    return 1;
}

/** [Implements IConditionalDefaultsModel] */
CClassConstSP RFLMixtureDefaultsModel::engineParamsType() const
{
    return RflMixtureParameters::TYPE;
}

/** Override clone() to copy the local cache */
IObject* RFLMixtureDefaultsModel::clone() const
{
    IObject* myCopy = CObject::clone();
    RFLMixtureDefaultsModel& rflDefaultsModel = 
        dynamic_cast<RFLMixtureDefaultsModel&>(*myCopy);
    // copy the cache !
    rflDefaultsModel.barriersCache = barriersCache;
    return myCopy;
}

/* external symbol to allow class to be forced to be linked in */
bool RFLMixtureDefaultsModelLoad(){
    return (RFLMixtureDefaultsModel::TYPE != 0);
}

/** Key constructor (for the barriers cache) */
RFLMixtureDefaultsModel::Key::Key(double defaultProba,
    IObjectConstSP rflParams):
        defaultProba(defaultProba),
        rflParams(rflParams) {}

/** Key "==" operator (for the barriers cache) */
bool RFLMixtureDefaultsModel::Key::operator==(const Key& key) const
{
    return (defaultProba == key.defaultProba &&
        rflParams->equalTo(key.rflParams.get()));
}

/** KeyHash "()" operator (for the barriers cache) */
size_t RFLMixtureDefaultsModel::KeyHash::operator()(const Key& key) const 
{
    return (CDouble::hashCode(key.defaultProba) ^
        key.rflParams->hashCode());
}

DRLIB_END_NAMESPACE
