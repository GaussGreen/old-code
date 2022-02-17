//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : FAMQDistributionCMCDS.hpp
//
//   Description : Class for forward measure probability distribution for CMCDS
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#ifndef _FAMQDISTRIBUTIONCMCDS_HPP
#define _FAMQDISTRIBUTIONCMCDS_HPP

#include "edginc/FAMQDistribution.hpp"
#include "edginc/MultiQDistribution.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

//-------------------------------
// FAMQDistributionCMCDS 
//-------------------------------
class UTIL_DLL FAMQDistributionCMCDS : public virtual FAMQDistribution,
                                       public CObject
{
public:
    static CClassConstSP const TYPE;

    FAMQDistributionCMCDS(MultiQDistributionSP mq,
                          double convexA, 
                          double convexP, 
                          double convexT,
                          double convexRecRate, 
                          double convexRiskFreeRate, 
                          double delayA, 
                          double delayP, 
                          double delayT, 
                          int RichardsonOrder);

    virtual ~FAMQDistributionCMCDS();

    //--------------------
    // IDistribution method
    //--------------------
    virtual double pdf(double x) const;
    
protected:
    virtual int RichardsonPolyOrder();
    virtual double norm();
    virtual void calcNorm();
    
private:
    int RichardsonOrder;
    double normNumber;
    double convexRecoveryRate;
    double convexRiskFreeRate;
    
    FAMQDistributionCMCDS();
    static void load(CClassSP& clazz);
    static IObject* defaultFAMQCMCDS();
    
};

typedef smartConstPtr<FAMQDistributionCMCDS> FAMQDistributionCMCDSConstSP;
typedef smartPtr<FAMQDistributionCMCDS>      FAMQDistributionCMCDSSP;

DRLIB_END_NAMESPACE
#endif
