//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : VolProcessedMQCMS.hpp
//
//   Description : MultiQ processed vol for CMS
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#ifndef _VOLPROCESSEDMQCMS_HPP
#define _VOLPROCESSEDMQCMS_HPP

#include "edginc/VolProcessedMQ.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IDistribution1D);


class MARKET_DLL VolProcessedMQCMS : public VolProcessedMQ
{
public:
    static CClassConstSP const TYPE;
    virtual ~VolProcessedMQCMS();

    VolProcessedMQCMS(string name,
                      TimeMetricSP metric,
                      IDistribution1D *mq, 
                      IDistribution1D *famq);

protected:
    VolProcessedMQCMS(CClassConstSP clazz = TYPE);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultProcessed();
};

DECLARE(VolProcessedMQCMS);

DRLIB_END_NAMESPACE

#endif
