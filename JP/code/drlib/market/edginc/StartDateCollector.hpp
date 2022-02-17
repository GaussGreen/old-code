//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StartDateCollector.hpp
//
//   Description : Start date collector class
//
//   Author      : Andre Segger
//
//   Date        : 27 Apr 2001
//
//
//----------------------------------------------------------------------------

#ifndef STARTDATE_COLLECT_HPP
#define STARTDATE_COLLECT_HPP
#include "edginc/Collector.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** A class to validate cross validate start dates */
class MARKET_DLL StartDateCollector: public CObject, public virtual ICollector {
public:
    friend class StartDateCollHelper;
    static CClassConstSP const TYPE;

    void startDateValidate(const DateTime& dateToValidate,
                           const string&   text,
                           bool            checkStrictlyEquals);

    StartDateCollector(const DateTime& startDate,
                       const string&   source,
                       bool            fwdStarting);
private:
    static void load(CClassSP& clazz);
    StartDateCollector(const StartDateCollector &rhs);
    StartDateCollector& operator=(const StartDateCollector& rhs);

    DateTime    startDate; // $unregistered
    string      source; // $unregistered
    bool        fwdStarting; // $unregistered
};

typedef smartPtr<StartDateCollector> StartDateCollectorSP;
typedef smartPtr<const StartDateCollector> StartDateCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
