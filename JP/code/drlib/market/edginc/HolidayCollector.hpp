//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HolidayCollector.hpp
//
//   Description : holiday collector class
//
//   Author      : André Segger
//
//   Date        : 17 Jul 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_HOLIDAY_COLLECT_HPP
#define EDR_HOLIDAY_COLLECT_HPP
#include "edginc/smartPtr.hpp"
#include "edginc/Collector.hpp"
#include "edginc/Holiday.hpp"


DRLIB_BEGIN_NAMESPACE

/** A class to retrieve the market holiday from an asset */
class MARKET_DLL HolidayCollector: public CObject, 
                        public virtual ICollector {
public:
    friend class HolidayCollHelper;
    static CClassConstSP const TYPE;
    void      setHoliday(const HolidayConstSP& holiday);
    HolidayConstSP getHoliday();
    HolidayCollector();

private:
    static void load(CClassSP& clazz);
    HolidayCollector(const HolidayCollector &rhs);
    HolidayCollector& operator=(const HolidayCollector& rhs);

    HolidayConstSP  holiday; // $unregistered
};

typedef smartPtr<HolidayCollector> HolidayCollectorSP;
typedef smartPtr<const HolidayCollector> HolidayCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
