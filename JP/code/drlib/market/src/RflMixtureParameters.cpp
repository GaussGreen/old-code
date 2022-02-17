//----------------------------------------------------------------------------
//
//   Group       : QR Credit Hybrids
//
//   Description : Container class for RFL-specific per name model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Format.hpp"
#include "edginc/MappingFunction.hpp"
#include "edginc/PiecewiseFlatMappingFunction.hpp"
#include "edginc/PiecewiseFlatIncrementalMappingFunction.hpp"
#include "edginc/RFLParameters.hpp"
#include "edginc/PiecewiseLinearMappingFunction.hpp"
#include "edginc/Nrfns.hpp" //for zbrentUseful root finder
#include "edginc/PortfolioName.hpp" // for the betaTweak function
#define QLIB_RFLMIXTUREPARAMETERS_CPP
#include "edginc/RflMixtureParameters.hpp"

DRLIB_BEGIN_NAMESPACE


////////////////////////// RflMixtureParameters members //////////////////////////

// Return sum of absolute values of discontinuities
double RflMixtureParameters::
ParameterBadness() const
{
  if (fabs(m_badnessWeight) <= 1.0e-10 || m_betas.empty())
    return 0.0;
  double badness = 0.0;
  CDoubleArray::const_iterator iB = m_betas.begin();
  CDoubleArray::const_iterator iB_ = iB++;
  for (; iB != m_betas.end(); ++iB)    
  {
    badness += fabs(*iB - *iB_);
    iB_ = iB;
  }
  return m_badnessWeight * badness;
}

// Return systemic term for given factor value
double RflMixtureParameters::
f(double factorValue) const
{
  // Due to our representation of thresholds as offsets we must contend with linear search
  // This should not be a performance issue since the number of thresholds is low (<5)
  CDoubleArray::const_iterator iT = m_thresholds.begin();
  CDoubleArray::const_iterator iB = m_betas.begin();
  CDoubleArray::const_iterator iA = m_alphas.begin();
  for (; iT != m_thresholds.end(); ++iT, ++iB, ++iA)
    if (factorValue <= *iT)
      return *iB * factorValue + *iA;
  return *iB * factorValue + *iA;
}

double RflMixtureParameters::
survivalProb(double barrier,   // Value of default barrier
             double c,   // displacement of mixture components (usually < 0)
             double w) const // weight of displaced component (usually w << 1)
{
  // Check if we have a real mixture
  const bool mixture = (fabs(c) <= 1.0e-6 || fabs(w) <= 1.0e-10) ? false : true;

  const double cc = mixture? c : 0.0;
  const double ww = mixture? w : 0.0;

  double denom = sqrt(1.0 + m_betas[0] * m_betas[0]);
  double t_ = m_thresholds[0] - cc;
  double a = (barrier - m_betas[0] * cc - m_alphas[0]) / denom;
  double r = m_betas[0] / denom;
  
  double p = N2Std(a, t_, r);

  const size_t Nt = m_thresholds.size();
  size_t i = 1;
  for (; i < Nt; ++i)
  {
    denom = sqrt(1.0 + m_betas[i] * m_betas[i]);
    double t = m_thresholds[i] - cc;
    a = (barrier - m_betas[i] * cc - m_alphas[i]) / denom;
    r = m_betas[i] / denom;

    p += N2Std(a, t, r) - N2Std(a, t_, r);

    t_ = t;
  }

  denom = sqrt(1.0 + m_betas[i] * m_betas[i]);
  a = (barrier - m_betas[i] * cc - m_alphas[i]) / denom;
  r = m_betas[i] / denom;

  p += N1(a) - N2Std(a, t_, r);

  double q = 1.0 - p;

  if (!mixture)
    return q;
  else
    return w * q + (1.0 - w) * survivalProb(barrier);
}


/** compute the conditional survival probability for a name */
void RflMixtureParameters::condSurvivalProba(double barrier,
                                          double factorValue,
                                          double& pm) const
{
    pm = 1.0 - N1(barrier - f(factorValue));
}

// Write thresholds to provided container
void
RflMixtureParameters::
Thresholds(set<double>& thresholds) const
{
  thresholds.clear();
  double t_ = 0;
  for (size_t i = 0; i < m_thresholds.size(); ++i)
  {
    const double t = t_ + m_thresholds[i];
    thresholds.insert(t);
    t_ = t;
  }
}


// structure used only in the threshold solving
// contains a pointer to the class itself
// and the value of the probability p for which
// to find the threshold Fx-1(p)
typedef struct barrierHelper_ 
{
    const RflMixtureParameters *t;
    double         survProb;
    double         c, w;
} barrierHelper;

// function used in Brent solving for tgauss
// returns FX(tgauss) - pgauss
static double barrierFunctionToSolve(double x, void  *s) {
    barrierHelper *str = (barrierHelper *) s;
    return str->t->survivalProb(x, str->c, str->w) - str->survProb;
}

/** compute the threshold associated to the survival probability for a name */
void RflMixtureParameters::defaultBarrier(double survProb, 
                                       double &barrier, 
                                       double c,
                                       double w) const
{
    static char method[] = "RflMixtureParameters::threshold";
    barrierHelper str;
    str.t = this;
    str.survProb = survProb;
    str.c = c;
    str.w = w;

    // solve for x s.t. XCumProba(x) = pgauss
    double rangeLow  = -10.;
    double rangeHigh =  10.;
    double funcLow;
    double funcHigh;
    ZbracReturn zbrac = zbracUseful(&barrierFunctionToSolve,
                                    &str,
                                    &rangeLow,
                                    &rangeHigh,
                                    &funcLow,
                                    &funcHigh);

    // zbracUseful returns ZBRAC_SUCCESS if everything went well. Else:
    if (zbrac == ZBRAC_FUNC_ERROR) {
        throw NRException(method,
                          "Zbrac error evaluating the function in the "
                          "initial boundary [" + 
                          Format::toString(rangeLow) + ", " +
                          Format::toString(rangeHigh) + "].");
    }
    else if (zbrac == ZBRAC_BRAC_ERROR) {
        // Throw an NRException
        throw NRException(method,
                          "Zbrac error: Failed to identify the range"
                          " (last range tested: [" +
                          Format::toString("%.12f", rangeLow) + ", " + 
                          Format::toString("%.12f", rangeHigh) + "]).");
    }

    // Call zbrent
    try {
        barrier = zbrentUsefulBoundary(
                  &barrierFunctionToSolve, // Function to find the root of
                  &str,       // cds, for
                  rangeLow,  // Low value for x
                  rangeHigh, // High value for x
                  1.e-10,    // Tolerance
                  funcLow,
                  funcHigh);
    }
    catch (NRException&) {
        // This is zbrent failing to solve - hide "zbrent" message,
        // which confuses users, and output a more meaningful error
        throw NRException(
             method,
             "Zbrent error: failed to find tgauss");
    }
}

/** Returns the name of this object */
string RflMixtureParameters::getName() const 
{
    return name;
}

/** Invoked when Class is 'loaded' */
void RflMixtureParameters::load(CClassSP& clazz) 
{
    clazz->setPublic();
    REGISTER(RflMixtureParameters, clazz);
    SUPERCLASS(RationalisedCreditEngineParameters);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(name, "Identifier");

    FIELD(m_thresholds, "First threshold and additional increments for factor loadings (systemic term)");
    FIELD(m_betas, "Factor loading values (slopes for piecewise linear systemic term)");
    FIELD(m_alphas, "Intercepts for piecewise linear systemic term [all 0.0]");
    FIELD_MAKE_OPTIONAL(m_alphas);
    FIELD(m_badnessWeight, "Calibration weight of parameter badness [0.0].");
    FIELD_MAKE_OPTIONAL(m_badnessWeight);

    Calibrator::IAdjustable::registerField(
        clazz,
        "m_thresholds",
        new Range(OpenBoundary(0.0), ClosedBoundary(1.0)));

    Calibrator::IAdjustable::registerField(
        clazz,
        "m_betas",
        new Range(OpenBoundary(0.0), ClosedBoundary(10.0)));

    Calibrator::IAdjustable::registerField(
        clazz,
        "m_alphas",
        new Range(ClosedBoundary(0.0), ClosedBoundary(1.0)));
}


void RflMixtureParameters::validatePop2Object()
{
  // Set alphas to have correct size
  m_alphas.clear();
  m_alphas.reserve(m_betas.size());
  for (size_t i = 0; i < m_betas.size(); ++i)
    m_alphas.push_back(0.0);
}  

void RflMixtureParameters::fieldsUpdated(const CFieldArray& fields)
{
    bool betasUpdated = false;
    
    string betasFieldName = getClass()->getDeclaredField("m_betas")->getName();
    for (unsigned int i = 0; i < fields.size() && !betasUpdated; ++i) 
    {
        string updatedFieldName = fields[i]->getName();
        betasUpdated = (updatedFieldName == betasFieldName);
	  }
    
    if (betasUpdated ) 
    {
        // update alphas to have correct size
        m_alphas.clear();
        m_alphas.reserve(m_betas.size());
        for (size_t i = 0; i < m_betas.size(); ++i)
        m_alphas.push_back(0.0);
    }
}

RflMixtureParameters::RflMixtureParameters() :
    RationalisedCreditEngineParameters(TYPE),
      m_badnessWeight(0.0)
{}
    
/** Default constructor */
IObject* RflMixtureParameters::defaultConstructor() 
{
    return new RflMixtureParameters();
}

/** TYPE (for reflection) */        
CClassConstSP const RflMixtureParameters::TYPE =
    CClass::registerClassLoadMethod("RflMixtureParameters",
                                    typeid(RflMixtureParameters),
                                    RflMixtureParameters::load);

DEFINE_TEMPLATE_TYPE(RflMixtureParametersWrapper);


DRLIB_END_NAMESPACE
