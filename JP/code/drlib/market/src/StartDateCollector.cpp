//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StartDateCollector.cpp
//
//   Description : Start date collector class
//
//   Author      : Andre Segger
//
//   Date        : 27 Apr 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/StartDateCollector.hpp"

DRLIB_BEGIN_NAMESPACE
/** Invoked when TestCollect is 'loaded' */
void StartDateCollector::load(CClassSP& clazz){
    REGISTER(StartDateCollector, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICollector);
}

CClassConstSP const StartDateCollector::TYPE = CClass::registerClassLoadMethod(
    "StartDateCollector", typeid(StartDateCollector), load);

StartDateCollector::StartDateCollector(
    const DateTime& startDate,
    const string&   source,
    bool            fwdStarting): 
    CObject(TYPE), startDate(startDate),
    source(source), fwdStarting(fwdStarting) {}

void StartDateCollector::startDateValidate(const DateTime& dateToValidate,
                                           const string&   text,
                                           bool            checkStrictlyEquals)
{
    if (fwdStarting)
    {
        /* ensure option does not start before basket */
        if ( checkStrictlyEquals ) {
            if ( !dateToValidate.equals(startDate) ) {
                throw ModelException("StartDateCollector::startDateValidate", 
                                     text + " start date (" + dateToValidate.toString() + 
                                     ") must be the same as " + source + " start date (" +
                                     startDate.toString() +")");
            }
        } else {
            if ( dateToValidate.isGreater(startDate) ) {
                throw ModelException("StartDateCollector::startDateValidate",
                                     text + " start date (" + dateToValidate.toString() +
                                     ") can't be after " + source + " start date (" +
                                     startDate.toString() +")");
            }
        }
    }
    else
    {
        /* ensure basket starts on or before start date */
        if (dateToValidate.isGreater(startDate))
        {
            throw ModelException("StartDateCollector::startDateValidate",
                                 text + " start date (" + dateToValidate.toString() +
                                 ") can't be after " + source + " start date (" +
                                 startDate.toString() +")");
        }
    }
}

DRLIB_END_NAMESPACE

