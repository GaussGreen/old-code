//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : FAMQDistribution.hpp
//
//   Description : Base Class for forward measure probability distribution
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#ifndef _FAMQDISTRIBUTION_HPP
#define _FAMQDISTRIBUTION_HPP

#include "edginc/IDistribution1D.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(FAMQDistribution);
FORWARD_DECLARE(MultiQDistribution);


//-------------------------------
// FAMQDistribution
//-------------------------------
class UTIL_DLL FAMQDistribution : public IDistribution1D
{
public:
    static CClassConstSP const TYPE;

    FAMQDistribution(MultiQDistributionSP mq,
                     double convexA, double convexP, double convexT, double convexFreq,
                     double delayA, double delayP, double delayT, double delayFreq);

    virtual ~FAMQDistribution();

    //-----------------------------
    //IDistribution1D method
    //-----------------------------
    virtual double cdf(double x) const;

    virtual double expectation() const;

    virtual double expectation(const Function1DDouble& payoff) const;

    virtual double variance() const = 0;

    //-------------- new methods -------------
    // x * pdf(x)
    virtual double xpdf(double x) const;
    
    virtual double fpdf(const Function1DDouble& payoff, double x) const;

protected:
    //helper classes
    // f(x)
    class UTIL_DLL myFuncX : public Function1DDouble
    {
    public:
        myFuncX(double (FAMQDistribution::*psi) (double) const,
                const FAMQDistribution *fa,
                Range r);
        virtual double operator()(double x) const;
        virtual ~myFuncX() {};
        
    private:
        double (FAMQDistribution::*psi)(double) const;
        const FAMQDistribution *fa;
        Range r;
    };

    // f(g, x)
    class UTIL_DLL myFuncFx : public Function1DDouble
    {
    public:
        myFuncFx(double (FAMQDistribution::*psi) (const Function1DDouble&, double) const,
                 const Function1DDouble& f,
                 const FAMQDistribution *fa,
                 Range r);

    virtual double operator()(double x) const;
    virtual ~myFuncFx() {};

    private:
        double (FAMQDistribution::*psi)(const Function1DDouble&, double) const;
        const FAMQDistribution *fa;
        const Function1DDouble& f;
        Range r;
    };

protected:
    FAMQDistribution(CClassConstSP = TYPE);

    virtual int RichardsonPolyOrder() const = 0;
    virtual double norm() const = 0;

    static const double FAMQ_TINY;
    static const double FAMQ_LARGE_YIELD;
    static const int    FAMQ_MAX_STEPS;

    MultiQDistributionSP mq;
    double convexA;
    double convexP;
    double convexT;
    double convexFreq;
    double delayA;
    double delayP;
    double delayT;
    double delayFreq;


private:
    static void load(CClassSP& clazz);

};

DECLARE(FAMQDistribution);

DRLIB_END_NAMESPACE
#endif
