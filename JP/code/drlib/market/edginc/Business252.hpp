//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Business252.hpp
//
//   Description : Bus/252 day count convention
//
//   Author      : Xiaolan Zhang
//
//   Date        : 7 November 2005
//
//----------------------------------------------------------------------------

#ifndef BUSINESS252_HPP
#define BUSINESS252_HPP

#include "edginc/Holiday.hpp"
#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

/** Bus/252 day count convention */
class MARKET_DLL Business252 : public DayCountConvention,
                               public IGetMarket
{
public:
    static CClassConstSP const TYPE;
    friend class Business252Helper;

    Business252();
    Business252(const string& name);
    Business252(const HolidayConstSP& hols);

    virtual ~Business252();
    
    /** how many days between two dates */
    virtual int days(const DateTime &lodate, const DateTime &hidate) const;
    
    /** year fraction between two dates */
    virtual double years(const DateTime &lodate, const DateTime &hidate) const;
    
    /** returns a string description e.g. Bus/252 */
    virtual string toString() const;

    /** populate with supplied holidays */
    void setHoliday(const HolidayConstSP& hols) const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Called after reflection is used to fill the internal fields */
    virtual void validatePop2Object();

private:
    mutable bool holsSpecified;
    mutable MarketWrapper<Holiday> hols;
};

DRLIB_END_NAMESPACE

#endif
