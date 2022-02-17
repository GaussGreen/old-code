//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : Interface to an object capable of performing 
//                 bad day adjustmens
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#ifndef IBADDAYADJUSTER_HPP
#define IBADDAYADJUSTER_HPP

#include "edginc/DateTime.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/** Interface for objects capable of performing bad day adjustmens */
class MARKET_DLL IBadDayAdjuster : public virtual IObject {
public:
    static CClassConstSP const TYPE;
    virtual ~IBadDayAdjuster();

    /** Returns "date" bad day adjusted using the bad day convention
     * and holidays in this object */
    virtual DateTime badDayAdjust(const DateTime& date) const = 0;

    /** Add a number of business days to a date */
    virtual DateTime addBusinessDays(const DateTime& from, int busDays) const = 0;

protected:
    IBadDayAdjuster();

private:
    static void load(CClassSP& clazz);
};

DECLARE(IBadDayAdjuster);

DRLIB_END_NAMESPACE

#endif
