//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : FAMQDistributionCMS.cpp
//
//   Description : Forward measure probability distribution for CMS
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/FAMQDistributionCMS.hpp"
#include "edginc/MultiQDistribution.hpp"
#include "edginc/Integrator.hpp"
#include "edginc/mathlib.hpp"
#include <math.h>

DRLIB_BEGIN_NAMESPACE

FAMQDistributionCMS::FAMQDistributionCMS(CClassConstSP clazz)
    :CObject(clazz)
{}


FAMQDistributionCMS::~FAMQDistributionCMS()
{}

FAMQDistributionCMS::FAMQDistributionCMS(MultiQDistributionSP mq,
                                         double convexA, 
                                         double convexP, 
                                         double convexT, 
                                         double convexFreq,
                                         double delayA, 
                                         double delayP, 
                                         double delayT, 
                                         double delayFreq,
                                         int RichardsonOrder):
    CObject(TYPE),
    FAMQDistribution(mq, convexA,  convexP, convexT, convexFreq,
                     delayA,  delayP, delayT, delayFreq), 
    RichardsonOrder(RichardsonOrder)
{
    calcNorm();
}

void FAMQDistributionCMS::calcNorm()
{
    static const string method = "FAMQDistributionCMS::calcNorm";
    double sum = 0;

    OpenRomberg1D integ(FAMQ_TINY, FAMQ_MAX_STEPS, RichardsonOrder);
    OpenBoundary inf(FAMQ_LARGE_YIELD);
    OpenBoundary lowerB(FAMQ_TINY);
    Range r(lowerB, inf);

    
    
    myFuncX integrand(&FAMQDistribution::pdf, this, r);

    try
    {
        sum = integ.integrate(integrand);
    }
    catch(exception &e) 
    {
        throw ModelException(&e, method);
    }
    
    normNumber = sum;
} 

int FAMQDistributionCMS::RichardsonPolyOrder() const
{
    return RichardsonOrder;
}

double FAMQDistributionCMS::norm() const
{
    return normNumber;
}

//---------------------------------
// IDistribution1D method
//---------------------------------

double FAMQDistributionCMS::pdf(double y) const
{
    static const char method[] = "FAMQDistributionCMS::pdf";
    static const double Q3_FA_BDRY = 1.0e-2;
    static const double TINY = 1.0e-10;

    double yCutoff = Q3_FA_BDRY * mq->forward();
    double yCut = Maths::max(y, yCutoff);
    double p = Maths::min(y/yCutoff, 1.0);

    double z, x, d, ibpv, zero;

    try
    {
        //annuity
        if(y > TINY)
        {
            z = pow(yCut, convexP);
            x = 1 + convexA * z / convexFreq;
            d = pow(x, -convexFreq * convexT);
            ibpv = yCut / (1 - d);
            ibpv = p * ibpv + (1-p) / convexT;
        } else {
            ibpv = 1 / convexT;
        }

        //delay
        if(y > TINY)
        {
            yCut = max(y, yCutoff);
            z = pow(yCut, delayP);
            x = 1 + delayA * z / delayFreq;
            zero = 1 / pow(x, delayFreq * delayT);
            zero = p * zero + (1-p);
        } else {
            zero = 1;
        }

        double dens = mq->pdf(y) * ibpv * zero;
        //double dens = N1Density(y) * ibpv * zero;
        return dens;

    }
    catch(exception &e) 
    {
        throw ModelException(&e, method);
    }

}

double FAMQDistributionCMS::mgf(double z, int n) const
{
    throw ModelException("FAMQDistributionCMS::mgf", "not implemented yet");
    return 0;
}

double FAMQDistributionCMS::variance() const
{
    throw ModelException("FAMQDistributionCMS::variance", "not implemented yet");
    return 0;
}

DiscreteDistributionConstSP FAMQDistributionCMS::discretise() const
{
    throw ModelException("FAMQDistributionCMS::discretise", "not implemented yet");
    //return NULL;
}


//---------------------------------------------------------------------
//Invoked when class is loaded
//---------------------------------------------------------------------
void FAMQDistributionCMS::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(FAMQDistributionCMS, clazz);
    SUPERCLASS(FAMQDistribution);
    EMPTY_SHELL_METHOD(defaultFAMQCMS);
    FIELD(RichardsonOrder, "");
    FIELD(normNumber, "");
    FIELD_MAKE_TRANSIENT(normNumber);
}

IObject* FAMQDistributionCMS::defaultFAMQCMS() {
    return new FAMQDistributionCMS();
}


CClassConstSP const FAMQDistributionCMS::TYPE = CClass::registerClassLoadMethod
    ("FAMQDistributionCMS", 
     typeid(FAMQDistributionCMS), 
     load);


DRLIB_END_NAMESPACE
