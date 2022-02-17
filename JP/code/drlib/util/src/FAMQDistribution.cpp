//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : FAMQDistribution.cpp
//
//   Description : Forward measure MQ probability distribution
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/FAMQDistribution.hpp"
#include "edginc/Range.hpp"
#include "edginc/Integrator.hpp"
#include "edginc/MultiQDistribution.hpp"

DRLIB_BEGIN_NAMESPACE


//-------------------------------
// Inner classes
//-------------------------------
FAMQDistribution::myFuncX::myFuncX(double (FAMQDistribution::*psi) (double) const,
                                   const FAMQDistribution *fa,
                                   Range r)
    :
    Function1DDouble(r), psi(psi), fa(fa), r(r) 
{}

double FAMQDistribution::myFuncX::operator() (double x) const
{
    return (fa->*psi) (x);
}

FAMQDistribution::myFuncFx::myFuncFx(double (FAMQDistribution::*psi) (const Function1DDouble&, double) const,
                                     const Function1DDouble& f,
                                     const FAMQDistribution *fa,
                                     Range r)
    :
    Function1DDouble(r), psi(psi), fa(fa), f(f), r(r) 
{}

double FAMQDistribution::myFuncFx::operator() (double x) const
{
    return (fa->*psi) (f, x);
}

//-------------------------------
// FAMQDistribution
//-------------------------------

FAMQDistribution::FAMQDistribution(CClassConstSP clazz)
{}


FAMQDistribution::FAMQDistribution(MultiQDistributionSP mq,
                                   double convexA, 
                                   double convexP, 
                                   double convexT, 
                                   double convexFreq,
                                   double delayA, 
                                   double delayP, 
                                   double delayT, 
                                   double delayFreq)
    :
    mq(mq), convexA(convexA), convexP(convexP), convexT(convexT), convexFreq(convexFreq),
    delayA(delayA), delayP(delayP), delayT(delayT), delayFreq(delayFreq)
{}


FAMQDistribution::~FAMQDistribution()
{}


double FAMQDistribution::xpdf(double x) const
{
    return x * pdf(x);
}

double FAMQDistribution::fpdf(const Function1DDouble& payoff, double x) const
{
    return (payoff)(x) * pdf(x);
}

double FAMQDistribution::cdf(double x) const       //todo
{
    static const string method = "FAMQDistribution::cdf";
    double sum;

    OpenRomberg1D integ(FAMQ_TINY, FAMQ_MAX_STEPS, RichardsonPolyOrder() );
    OpenBoundary upperB(x);
    OpenBoundary lowerB(FAMQ_TINY);
    Range r(lowerB, upperB);

    myFuncX integrand(&FAMQDistribution::pdf, this, r);

    try
    {

        sum = integ.integrate(integrand);
    }
    catch(exception &e) 
    {
        throw ModelException(&e, method);
    }

    return sum / norm();
}



double FAMQDistribution::expectation() const
{
    static const string method = "FAMQDistribution::expection";
    double sum;

    OpenRomberg1D integ(FAMQ_TINY, FAMQ_MAX_STEPS, RichardsonPolyOrder() );
    OpenBoundary inf(FAMQ_LARGE_YIELD);
    OpenBoundary lowerB(FAMQ_TINY);
    Range r(lowerB, inf);

    myFuncX integrand(&FAMQDistribution::xpdf, this, r);

    try
    {
        sum = integ.integrate(integrand);
    }
    catch(exception &e) 
    {
        throw ModelException(&e, method);
    }

    return sum / norm();
}


double FAMQDistribution::expectation(const Function1DDouble& payoff) const
{
    static const string method = "FAMQDistribution::expection";
    double sum;

    OpenRomberg1D integ(FAMQ_TINY, FAMQ_MAX_STEPS, RichardsonPolyOrder() );
    OpenBoundary inf(FAMQ_LARGE_YIELD);
    OpenBoundary lowerB(FAMQ_TINY);
    Range r(lowerB, inf);

    myFuncFx integrand(&FAMQDistribution::fpdf, payoff, this, r);

    try
    {
        sum = integ.integrate(integrand);
    }
    catch(exception &e) 
    {
        throw ModelException(&e, method);
    }

    return sum / norm();
}

void FAMQDistribution::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(FAMQDistribution, clazz);
    IMPLEMENTS(IDistribution1D);
    FIELD(mq, "MultiQ");
    FIELD(convexA, "");
    FIELD(convexP,"");
    FIELD(convexT,"");
    FIELD(convexFreq,"");
    FIELD(delayA,"");
    FIELD(delayP,"");
    FIELD(delayT,"");
    FIELD(delayFreq,"");
}

//static const
const double FAMQDistribution::FAMQ_TINY        = 1.0e-6;
const double FAMQDistribution::FAMQ_LARGE_YIELD = 0.5;
const int    FAMQDistribution::FAMQ_MAX_STEPS   = 500;

CClassConstSP const FAMQDistribution::TYPE = CClass::registerClassLoadMethod
    ("FAMQDistribution", 
     typeid(FAMQDistribution), 
     load);

DRLIB_END_NAMESPACE
