//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : VolProcessedMQ.hpp
//
//   Description : Base class for MultiQ processed vol
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#ifndef _VOLPROCESSEDMQ_HPP
#define _VOLPROCESSEDMQ_HPP

#include "edginc/VolProcessed.hpp"
#include "edginc/Function.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IDistribution1D.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IDistribution1D);


/********************************************************
 ** class IVolProcessedMQ  
 **
 ********************************************************
 */
class MARKET_DLL VolProcessedMQ: public virtual IVolProcessed,
                                 public CObject {
public:
    static CClassConstSP const TYPE;
    virtual ~VolProcessedMQ();

    VolProcessedMQ(string name,
                   TimeMetricSP metric,
                   IDistribution1D *mq, 
                   IDistribution1D *famq);

    //-----------------------
    // IVolProcessed methods
    //-----------------------
    virtual string getName() const;
    virtual double calcTradingTime(const DateTime &date1, const DateTime &date2) const;
    virtual TimeMetricConstSP GetTimeMetric() const;

    //-----------------------
    // VolProcessedMQ methods
    //-----------------------
    /** simply return adjusted fwd */
    virtual double fwd();

    /** return fwd that can be capped, floored, or a spread added */
    virtual double fwd(double capValue, 
                       double floorValue, 
                       double spraed, 
                       bool isAdditive);

    /** return fwd that can be manipulated in any way according "payoff" */
    virtual double fwd(const Function1DDouble& payoff);

    /** return vanilla option price */
    //virtual double vanillaOptionPrice(bool isCall,
    //                                  double strike) const;
    //
    //virtual double binaryOptionPrice(bool isCall,
    //                                 double strike) const;

protected:
    VolProcessedMQ(CClassConstSP clazz = TYPE);

protected:
    string       name;
    TimeMetricSP metric;
    IDistribution1DSP mq;
    IDistribution1DSP famq;

private:
    static void load(CClassSP& clazz);
};

DECLARE(VolProcessedMQ);

DRLIB_END_NAMESPACE

#endif
