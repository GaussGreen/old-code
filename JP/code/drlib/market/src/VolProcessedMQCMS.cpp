//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : VolProcessedMQCMS.cpp
//
//   Description : class for MultiQ processed vol for CMS
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/VolProcessedMQCMS.hpp"
#include "edginc/IDistribution1D.hpp"

DRLIB_BEGIN_NAMESPACE

/* Make sure the class links */
bool  VolProcessedMQCMSLoad() {
    return (VolProcessedMQCMS::TYPE != 0);
}


/********************************************************
 ** class VolProcessedMQCMS  
 ********************************************************
 */

VolProcessedMQCMS::VolProcessedMQCMS(CClassConstSP clazz)
    : VolProcessedMQ(clazz)
{}

VolProcessedMQCMS::~VolProcessedMQCMS() 
{}


VolProcessedMQCMS::VolProcessedMQCMS(string name,
                                     TimeMetricSP metric,
                                     IDistribution1D *mq, 
                                     IDistribution1D *famq)
    : VolProcessedMQ(name, metric, mq, famq)
{}

//------------------------------------------
// load
//------------------------------------------

void VolProcessedMQCMS::load(CClassSP& clazz) {
    REGISTER(VolProcessedMQCMS, clazz);
    SUPERCLASS(VolProcessedMQ);
    EMPTY_SHELL_METHOD(defaultProcessed);
}

IObject* VolProcessedMQCMS::defaultProcessed()
{
    return new VolProcessedMQCMS();
}

CClassConstSP const VolProcessedMQCMS::TYPE = 
    CClass::registerClassLoadMethod("VolProcessedMQCMS", typeid(VolProcessedMQCMS), load);


DRLIB_END_NAMESPACE
