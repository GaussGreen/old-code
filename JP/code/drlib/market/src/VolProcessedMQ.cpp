//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : VolProcessedMQ.cpp
//
//   Description : Base class for MultiQ processed vol
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/VolProcessedMQ.hpp"
#include "edginc/Function.hpp"
#include "edginc/FAMQDistribution.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/* Make sure the class links */
bool  VolProcessedMQLoad() {
    return (VolProcessedMQ::TYPE != 0);
}



/********************************************************
 ** class MyFunc
 ********************************************************
 */
class MyFunc : public Function1DDouble
{
public:
    MyFunc(double capValue, 
           double floorValue, 
           double spread, 
           bool isAdditive)
        :
        Function1DDouble(),
        capValue(capValue), 
        floorValue(floorValue), 
        spread(spread),
        isAdditive(isAdditive)
    {};

    virtual double operator() (double x) const
    {
        if(isAdditive) 
        {
            x += spread;
        } else {
            x *= spread;
        }

        double y = Maths::min(Maths::max(x, floorValue), capValue);

        return y;
    };

private:
    double capValue;
    double floorValue;
    double spread;
    bool   isAdditive;
};


/********************************************************
 ** class IVolProcessedMQ  
 **
 ********************************************************
 */

VolProcessedMQ::VolProcessedMQ(CClassConstSP clazz)
    :CObject(clazz)
{}

VolProcessedMQ::~VolProcessedMQ() 
{}


VolProcessedMQ::VolProcessedMQ(string name,
                               TimeMetricSP metric,
                               IDistribution1D *mq, 
                               IDistribution1D *famq)
    :CObject(TYPE), name(name), metric(metric), mq(mq), famq(famq)
{}


//------------------------------------------
// fwds
//------------------------------------------
double VolProcessedMQ::fwd()
{
    return famq->expectation();
}

double VolProcessedMQ::fwd(double capValue, 
                           double floorValue, 
                           double spread, 
                           bool isAdditive) 
{
    //make a Function1DDoubleSP with take formula "min(max(x, floorValue),capValue)+spread"
    MyFunc payoff(capValue, floorValue, spread, isAdditive);
    return famq->expectation(payoff);
}

double VolProcessedMQ::fwd(const Function1DDouble& payoff) 
{
    return famq->expectation(payoff);
}


//------------------------------------------
// multiQ related methods
//------------------------------------------
/**
   void VolProcessedMQ::MQ2FAMQ()
   {
   FAData convexFA = Q3VNFMZero2Swap(zeroMat, zeroRate);
   FAData delayFA = Q3VNFMZero2Swap(payMat, payZeroRate);
   
   famq = FAMQDistributionCMSSP(new FAMQDistributionCMS(mq, convexA, convexP, convexT, 
   convexFreq, delayA, delayP, delayT, delayFreq));
   }
*/

//-----------------------
// IVolProcessed methods
//-----------------------
string VolProcessedMQ::getName() const
{
    return name;
}

double VolProcessedMQ::calcTradingTime(const DateTime &date1, const DateTime &date2) const
{
    return metric->yearFrac(date1, date2);
}

TimeMetricConstSP VolProcessedMQ::GetTimeMetric() const
{
    return TimeMetricConstSP(metric);
}



//------------------------------------------
// load
//------------------------------------------

void VolProcessedMQ::load(CClassSP& clazz) {
    REGISTER(VolProcessedMQ, clazz);
    IMPLEMENTS(IVolProcessed);
    SUPERCLASS(CObject);
    FIELD(name, "name");
    FIELD(metric, "time metric");
    FIELD(mq, "MultiQDistribution");
    FIELD(famq, "Forward measure MultiQDistribution");
}

CClassConstSP const VolProcessedMQ::TYPE = 
    CClass::registerClassLoadMethod("VolProcessedMQ", typeid(VolProcessedMQ), load);


DRLIB_END_NAMESPACE
