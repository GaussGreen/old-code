//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : FAMQDistributionCMS.hpp
//
//   Description : Class for forward measure probability distribution for CMS
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#ifndef _FAMQDISTRIBUTIONCMS_HPP
#define _FAMQDISTRIBUTIONCMS_HPP

#include "edginc/FAMQDistribution.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(MultiQDistribution);


//-------------------------------
// FAMQDistributionCMS 
//-------------------------------
class UTIL_DLL FAMQDistributionCMS : public CObject,
                                     public FAMQDistribution
{
public:
    static CClassConstSP const TYPE;

    FAMQDistributionCMS(MultiQDistributionSP mq,
                        double convexA, 
                        double convexP, 
                        double convexT, 
                        double convexFreq, 
                        double delayA, 
                        double delayP, 
                        double delayT, 
                        double delayFreq,
                        int RichardsonOrder);

    virtual ~FAMQDistributionCMS();
    //--------------------
    // IDistribution1D method
    //--------------------
    virtual double pdf(double x) const;
    virtual double mgf(double z, int n) const;
    virtual double variance() const;
    virtual DiscreteDistributionConstSP discretise() const;

protected:
    virtual int RichardsonPolyOrder() const;
    virtual double norm() const;
    virtual void calcNorm();

private:
    int RichardsonOrder;
    double normNumber;

    FAMQDistributionCMS(CClassConstSP clazz = TYPE);
    static void load(CClassSP& clazz);
    static IObject* defaultFAMQCMS();

};

DECLARE(FAMQDistributionCMS);

DRLIB_END_NAMESPACE
#endif
