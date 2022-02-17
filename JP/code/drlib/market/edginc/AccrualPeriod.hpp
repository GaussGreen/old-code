//----------------------------------------------------------------------------------------
//
//   Filename    : AccrualPeriod.hpp
//
//   Description : A simple container for an accrual period (start & end dates)
//
//-----------------------------------------------------------------------------------------

#ifndef QLIB_ACCRUALPERIOD_HPP
#define QLIB_ACCRUALPERIOD_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL AccrualPeriod : public CObject
{
public:
    static CClassConstSP const TYPE;

    AccrualPeriod(); //empty dates
    AccrualPeriod(DateTime start,
                  DateTime end);

    /** Retrieve start date */
    DateTime startDate();

    /** Retrieve end date */
    DateTime endDate();

    ~AccrualPeriod();

private:
    //infrastructure support
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    //-------
    // Fields
    //-------
    DateTime start;
    DateTime end;
};

DECLARE(AccrualPeriod)

DRLIB_END_NAMESPACE

#endif
