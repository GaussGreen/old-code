//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DateBuilder.hpp
//
//   Description : Parameterized Date Builder
//
//   Author      : Manos Venardos
//
//   Date        : 23 February 2006
//
//
//----------------------------------------------------------------------------

#ifndef DATE_BUILDER_HPP
#define DATE_BUILDER_HPP

#include "edginc/DateTime.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/Model.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for parameterized date sequences */
class MARKET_DLL IDateBuilder: virtual public IObject,
                               virtual public IGetMarket {
public:
    static CClassConstSP const TYPE;

    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the array of dates */
    virtual DateTimeArraySP dates() const = 0;

    /** Returns the first date in the list */
    virtual const DateTime start() const = 0;

    /** Returns the first date in the list */
    virtual const DateTime end() const = 0;

    /** Returns the number of dates */
    virtual int size() const = 0;

    /** Returns a particular date from the final array of dates */
    virtual const DateTime date(int index) const = 0;

    /** Returns a particular date from the final array of dates 
        just calls date() */
    const DateTime operator[](int index) const;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IDateBuilder> IDateBuilderSP;

/*****************************************************************************/
/** Interface for parameterized date sequences 
    a special one for SPI cos we need to interact with assets 
    for holidays - all rather fiddly */
class MARKET_DLL ISPIDateBuilder: virtual public IDateBuilder {
public:
    static CClassConstSP const TYPE;

    // note the ccyHols and the exclude flag are NOT fields in the concrete 
    // implementation classes because then interface would have them in 
    // every date builder, i.e. the one for rebalance dates, lock in dates,
    // coupon dates etc
    // Note the exclude flag is here because some existing deals have rebal
    // on asset history and then ISDA.
    virtual void constructDates(const IMultiFactors*          assets,
                                const ObservationSourceArray& sources,
                                const Holiday*                ccyHols,
                                bool                          excludeAssetHols) = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<ISPIDateBuilder> ISPIDateBuilderSP;

DRLIB_END_NAMESPACE

#endif
