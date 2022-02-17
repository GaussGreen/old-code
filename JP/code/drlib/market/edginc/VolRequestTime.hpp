//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestTime.hpp
//
//   Description : Vol request for dealing with [trading] time
//
//   Author      : Mark A Robson
//
//   Date        : 1 Nov 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOLREQUESTTIME_HPP
#define EDR_VOLREQUESTTIME_HPP

#include "edginc/VolProcessed.hpp"
#include "edginc/VolRequest.hpp"

DRLIB_BEGIN_NAMESPACE
/** This is a request for a processed vol that supports the bare minimum 
    functions that the ProcessedVol class gives. It should be supported by
    any vol.*/
class MARKET_DLL VolRequestTime: public CVolRequest {
public:
    static CClassConstSP const TYPE;

    virtual ~VolRequestTime();
    VolRequestTime();

    /** Creates a VolProcessed which supports only those methods that 
        VolProcessed does (ie suitable for dealing with [trading] time) */
    static IVolProcessed* createVolProcessed(const string&     name, 
                                             TimeMetricConstSP metric);
    /** Creates a VolProcessed which supports only those methods that 
        VolProcessed does (ie suitable for dealing with [trading] time).
        Here, simple DateTime::yearFrac is used */
    static IVolProcessed* createVolProcessed(const string&     name);
private:
    class Processed;
    VolRequestTime(const VolRequestTime &rhs);
    VolRequestTime& operator=(const VolRequestTime& rhs);
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
};

// smart pointers for VolRequestTime
typedef smartConstPtr<VolRequestTime> VolRequestTimeConstSP;
typedef smartPtr<VolRequestTime> VolRequestTimeSP;
#ifndef QLIB_VOLREQUESTTIME_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<VolRequestTime>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<VolRequestTime>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<VolRequestTime>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<VolRequestTime>);
#endif

DRLIB_END_NAMESPACE

#endif
